library(DescTools)
library(tidyverse)
library(data.table)
library(taxonomizr)
library(dtplyr)
library(rlang)
library(ggforce)
library(openxlsx)
library(rioja)
library(patchwork)
library(wesanderson)

source("butterfly/dmg.R")
source("butterfly/get_calculate_plot_grid.R")
source("butterfly/perk.R")
source("butterfly/damage_est_function.R")
source("butterfly/perk_wrapper.R")
source("butterfly/perk_wrapper_function.R")
source("butterfly/get_dmg_decay_fit.R")
source("butterfly/median.R")
source("butterfly/filter.R")


# ----------------- ARGS ----------------- #

metadata_file <- "new_metadata.tsv"
metadmg_data <- "reps.Tjornen.metadmg_data.tsv"
bamfilter_data <- "reps.Tjornen.bam-filter_stats.tsv"

taxonomy_out <- "results/taxonomy/holi-fb-taxids.tsv.gz"
taxonomy_in <- "/projects/lundbeck/scratch/for_abigail/holi-fb-taxids.tsv.gz"

names <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
nodes <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
acc2taxid_file <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"

used_saved_taxid = TRUE

plot_directory <- "replicateplots/"

# -----------------  METADATA ----------------- #

metadata <- read.table(metadata_file,he=T,as.is=T)
metadata$cgg <- gsub("_", "", metadata$CGG_ID)
names(metadata)[2] <- "date" # !!!!!!!!!!!!!!!! this must be changed
metadata <- metadata[, c("cgg","date")]
metadata$sample_id <- paste0(metadata$cgg,"-","merge")
metadata$specific_feature <- "lake"
metadata$label_fig <- paste0(metadata$cgg,"-",metadata$date)


# -----------------  METADMG  ----------------- #

holi_data <- fread(metadmg_data, header=T, sep="\t", fill=T, nThread=20)
names(holi_data)[1] <- "file"
holi_data$cgg <- sapply(holi_data$file, function(x) strsplit(strsplit(x,"/")[[1]][2],"_")[[1]][1])
holi_data$label <- sapply(holi_data$file, function(x) gsub("^.*/([^/]+_[^_]+)_.+$", "\\1", x))

# Merge with metadata
holi_data <- inner_join(holi_data, metadata, by ="cgg")


filter_metadmg <- function(df, samples){
  holi_data_sp_euk <- df |>
    filter(rank == "species") |>
    filter(grepl("Eukaryota", taxa_path)) |>
    mutate(PlantAnimal = case_when(
      grepl("Viridiplantae", taxa_path) ~ "plant",
      grepl("Metazoa", taxa_path) ~ "animal",
    )) |>
    rename(tax_name = taxid, n_reads = nreads)

  # Let's get the damage fits using CCC

  dat_all <- dmg_fwd_CCC(holi_data_sp_euk, samples, ci = "asymptotic", nperm = 100, nproc = 24)

  # Define good and bad hits
  dat_filt <- inner_join(dat_all, holi_data_sp_euk) |>
    mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))
}

dat_filt <- filter_metadmg(holi_data, holi_data$label  |> unique())

# ----------------- BAMFILTER & TAXONOMY ----------------- #

read.names.sql(names, sqlFile = "nameNode.sqlite", overwrite=TRUE)
read.nodes.sql(nodes, sqlFile = "nameNode.sqlite", overwrite=TRUE)

# Bamfilter uses accessions but we need taxids
# Re-implementing this in c++ so we never have to do this part
if (used_saved_taxid == TRUE){
  tax_ids <- read_table(taxonomy_in)
} else {
  accs <- fb_data |>
    select(accession.version = reference) |>
    distinct() |>
    pull(accession.version)

  acc2taxid <- fread(acc2taxid_file, # sloooow
    tmpdir = "~/",
    nThread = 24,
    showProgress = TRUE)

  tax_ids <- acc2taxid %>%
    filter(accession.version %in% accs)
  write_tsv(tax_ids, taxonomy_out)
}

# Gets the unique taxids
tax_ids_lst <- tax_ids |>
  select(taxid) |>
  distinct() |>
  pull(taxid)

# Gets the taxonomic info for the taxids in the data
tax_data <- getTaxonomy(tax_ids_lst, "nameNode.sqlite") |>
  as_tibble() |>
  mutate(taxid = tax_ids_lst)


fb_data <- fread(bamfilter_data, header=T, sep="\t", fill=T, nThread=20)
names(fb_data)[1] <- "file"
fb_data$cgg <- sapply(fb_data$file, function(x) strsplit(x, "/")[[1]][9])
fb_data$label <- sapply(fb_data$file, function(x) gsub("^.*/([^/]+_[^_]+)_.+$", "\\1", x))

seqrun_df <- unique(fb_data[,c("file","cgg","label")])
seqrun_df$seq <- sapply(seqrun_df$file, function(x) strsplit(strsplit(x,"/")[[1]][10],"_")[[1]][2])
seqrun_df$seqdate <- sapply(seqrun_df$file, function(x) strsplit(strsplit(x,"/")[[1]][10],"_")[[1]][1])


# ----------------- THE BIG FILTER & STATS MACHINE ----------------- #


agg_stats <- get_stats(tax_ids, tax_data, fb_data, dat_filt, metadata=metadata, mode="replicates")
agg_stats <- inner_join(agg_stats,seqrun_df, by="label") # this is saved so we dont have to do all the above all the time
# read.csv("repagg.csv",h=T)
# ----------------- BASIC PLOTS  ----------------- #

# ----------------- ranks classified ----------------- #

ranks_to_plot <- holi_data %>%
  group_by(rank) %>%
  summarize(totalreads = sum(nreads)) %>%
  arrange(desc(totalreads)) %>% slice_head(n=10) %>% select(rank)

df <- holi_data %>%
  filter(grepl("Eukaryota", taxa_path)) %>%
  filter(rank %in% ranks_to_plot$rank) %>%
  group_by(label, rank) %>%
  summarize(numreads = sum(nreads)) %>%
  mutate(proportion = numreads / sum(numreads)) %>%
  inner_join(seqrun_df) %>% inner_join(metadata)

plot_total <- ggplot(df, aes(y = as.factor(date), x = numreads, fill = rank)) +
  geom_col() +
  ggh4x::facet_nested(~seqdate*seq) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_proportion <- ggplot(df, aes(y = as.factor(date), x = proportion, fill = rank)) +
  geom_col() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggh4x::facet_nested(~seqdate*seq) +
  scale_x_continuous(labels = scales::percent)


plot_total / plot_proportion +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(plot_directory, "/euk.ranks.png"))


# ----------------- good vs bad fit ----------------- #


goodbad <- agg_stats %>%
  group_by(seq, fit, date, seqdate,PlantAnimal) %>%
  summarise(total_nreads=sum(n_reads))

ggplot(goodbad, aes(y=as.factor(date), x=total_nreads, fill=fit)) +
  geom_col(position="stack") +
  ggh4x::facet_nested(PlantAnimal~seqdate*seq) +
  labs(y = "Date (CE)", x = "Total number of reads")
ggsave(paste0(plot_directory, "/good_vs_bad.png"))


goodbad <- agg_stats %>%
  group_by(seq, seqdate, date,PlantAnimal) %>%
  summarise(
    prop_good = sum(n_reads[fit == 'good']) / sum(n_reads),
    prop_bad = sum(n_reads[fit == 'bad']) / sum(n_reads),
    .groups = 'drop'
  ) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "fit",
    values_to = "proportion"
  )

ggplot(goodbad, aes(y=as.factor(date), x=proportion, fill=fit)) +
  geom_col(position="stack") +
  ggh4x::facet_nested(PlantAnimal~seqdate*seq) +
  labs(y = "Date (CE)", x = "Proportion of reads")
ggsave(paste0(plot_directory, "/good_vs_bad.proportions.png"))


# ----------------- damage model ----------------- #


qdata <- read.table("qdata.tsv",h=TRUE)
agg_stats <- inner_join(agg_stats, qdata) %>% mutate(status = ifelse(median_A_b < lwr, "fail", "pass"),
  status = ifelse(!is_first_dates & fit == "bad", "fail", status)) %>%
  mutate(alpha = ifelse(status=="pass",1,0.2))


goodbad <- agg_stats %>%
  group_by(seq, status, date, seqdate,PlantAnimal) %>%
  summarise(total_nreads=sum(n_reads))

ggplot(goodbad, aes(y=as.factor(date), x=total_nreads, fill=status)) +
  geom_col(position="stack") +
  ggh4x::facet_nested(PlantAnimal~seqdate*seq) +
  labs(y = "Date (CE)", x = "Total number of reads")
ggsave(paste0(plot_directory, "/pass_vs_fail.png"))


goodbad <- agg_stats %>%
  group_by(seq, seqdate, date, PlantAnimal) %>%
  summarise(
    prop_pass = sum(n_reads[status == 'pass']) / sum(n_reads),
    prop_fail = sum(n_reads[status == 'fail']) / sum(n_reads),
    .groups = 'drop'
  ) %>%
  pivot_longer(
    cols = starts_with("prop_"),
    names_to = "status",
    values_to = "proportion"
  )

ggplot(goodbad, aes(y=as.factor(date), x=proportion, fill=status)) +
  geom_col(position="stack") +
  ggh4x::facet_nested(PlantAnimal~seqdate*seq) +
  labs(y = "Date (CE)", x = "Proportion of reads")
ggsave(paste0(plot_directory, "/pass_vs_fail.proportions.png"))


# ----------------- plotting filtered data ----------------- #
get_top <- function(df, n, planimal){
  tmp <- df %>%
  filter(PlantAnimal == planimal) %>%
  group_by(genus) %>%
  summarize(total_reads = sum(n_reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = n)
  return(df %>% filter(genus %in% tmp$genus))
}

doplot <- function(df, g){
  d <- df %>% filter(genus == g)
  p <- ggplot(d) +
    geom_point(aes(y = as.factor(date), x = median_A_b, color = status, size = n_reads, alpha = ifelse(status == "pass", 1, 0.4)), show.legend = c(color = TRUE, alpha = FALSE)) +
    scale_alpha_identity() +
    ggh4x::facet_nested(~seqdate*seq) +
    labs(
      title = g,
      x = "Median Damage (A_b)",
      y = "Date (CE)"
    ) +
    xlim(0,0.5)
    print(p)
}

wrap_plot <- function(planimal, outname){
  df <- get_top(agg_stats,20,planimal)
  pdf(file = outname)
  for (g in unique(df$genus)){
    doplot(df,g)
  }
  dev.off()
}

wrap_plot("animal", paste0(plot_directory, "/damage.animals.20.pdf"))
wrap_plot("plant", paste0(plot_directory, "/damage.plants.20.pdf"))


# ----------------- heatmaps ----------------- #

# compare replicates with heatmaps, one box per genus, heat is nreads
# one with status pass, one with good fit only, one with rm==remove true
tech_replicates <- agg_stats %>%
  group_by(seqdate, cgg) %>%
  summarise(seq_combinations = list(seq)) %>%
  ungroup() %>%
  mutate(seq_combinations = map(seq_combinations, ~ sort(unique(.)))) %>%
  unnest_wider(seq_combinations, names_sep = "")

create_plot <- function(row) {
  print(row)
  seq1 <- row[["seq_combinations1"]]
  seq2 <- row[["seq_combinations2"]]
  this_cgg <- row[["cgg"]]
  seqdate <- row[["seqdate"]]


  
  d <- agg_stats %>% filter(cgg == this_cgg & seq %in% c(seq1, seq2) & status == "pass")
  p <- ggplot(d, aes(x = seq, y = genus, fill = n_reads)) +
    geom_tile() +
    geom_text(aes(label = n_reads), color = "white") +
    ggtitle(paste("CGG:", cgg, "SeqDate:", seqdate))
  
  return(p)
}

plots <- lapply(1:nrow(tech_replicates), function(i) create_plot(tech_replicates[i, ]))

pdf(paste0(plot_directory, "tech_replicates_heatmap.pdf"))
lapply(plots, print)
dev.off()
# ----------------- venn diagrams ----------------- #

# shared and unique genera
