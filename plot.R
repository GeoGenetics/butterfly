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
metadmg_data <- "Tjornen.metadmg_data.tsv"
bamfilter_data <- "Tjornen.bam-filter_stats.tsv"

taxonomy_out <- "results/taxonomy/holi-fb-taxids.tsv.gz"
taxonomy_in <- "/projects/lundbeck/scratch/for_abigail/holi-fb-taxids.tsv.gz"

names <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
nodes <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
acc2taxid_file <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"

used_saved_taxid = TRUE

neg_metadmg_data <- "negs/NEGS.Tjornen.metadmg_data.tsv"
neg_bamfilter_data <- "negs/NEGS.Tjornen.bam-filter_stats.tsv"

plot_directory <- "finalplots/"

# -----------------  METADATA ----------------- #

metadata <- read.table(metadata_file,he=T,as.is=T)
metadata$cgg <- gsub("_", "", metadata$CGG_ID)
names(metadata)[2] <- "date" # !!!!!!!!!!!!!!!! this must be changed
metadata <- metadata[, c("cgg","date")]
metadata$sample_id <- paste0(metadata$cgg,"-","merge")
metadata$specific_feature <- "lake"
metadata$label_fig <- paste0(metadata$cgg,"-",metadata$date)
metadata$label <- metadata$cgg


# -----------------  METADMG  ----------------- #

holi_data <- fread(metadmg_data, header=T, sep="\t", fill=T, nThread=20)
names(holi_data)[1] <- "sample_id"
holi_data$cgg <- sub(".*(CGG[0-9]+)\\..*", "\\1", holi_data$sample_id)
holi_data$label <- holi_data$cgg

# Merge with metadata
holi_data <- inner_join(holi_data, metadata, by ="label")

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

dat_filt <- filter_metadmg(holi_data, metadata$cgg |> unique())

# negatives 
neg_holi_data <- fread(neg_metadmg_data, header=T, sep="\t", fill=T, nThread=20)
names(neg_holi_data)[1] <- "sample_id"
neg_holi_data$label <- sapply(neg_holi_data$sample_id, function(x) strsplit(x,"_")[[1]][3])

neg_dat_filt <- filter_metadmg(neg_holi_data, neg_holi_data$label |> unique())

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


fb_data <- read_delim(bamfilter_data, delim = "\t", escape_double = FALSE, col_names = TRUE, trim_ws = TRUE)
names(fb_data)[1] <- "sample_id"
fb_data$label <- sub(".*(CGG[0-9]+)\\..*", "\\1", fb_data$sample_id)

neg_fb_data <- fread(neg_bamfilter_data, header=T, sep="\t", fill=T, nThread=20)
names(neg_fb_data)[1] <- "sample_id"
neg_fb_data$label <- sapply(neg_fb_data$sample_id, function(x) strsplit(x,"_")[[1]][4])


# ----------------- THE BIG FILTER & STATS MACHINE ----------------- #


agg_stats <- get_stats(tax_ids, tax_data, fb_data, dat_filt, metadata=metadata, mode="library")


# ----------------- NEGATIVE CALCULATOR ----------------- #

# Get genera in blanks 
neg_agg_stats <- get_stats(tax_ids, tax_data, neg_fb_data, neg_dat_filt, mode="negatives") %>%
  mutate(genus = factor(genus, levels = sort(unique(genus),decreasing=TRUE))) %>% 
  group_by(label) %>%
  ungroup()

ggplot(neg_agg_stats, aes(x = label, y = genus, fill = n_reads)) +
  geom_tile() +
  geom_text(aes(label = n_reads), color = "black", size = 3) +
  scale_fill_gradientn(colors = wes_palette("Zissou1"),values=scales::rescale(c(min(neg_agg_stats$n_reads),median(neg_agg_stats$n_reads),max(neg_agg_stats$n_reads)))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(paste0(plot_directory, "/negatives_heatmap.png"))

contam_genera <- neg_agg_stats %>%
    group_by(genus) %>%
    summarise(
      reads_blanks = sum(n_reads) 
    ) %>% arrange(desc(reads_blanks))



# ----------------- BASIC PLOTS  ----------------- #

# ----------- A) ranks classified  ----------- #
ranks_to_plot <- holi_data %>%
  group_by(rank) %>%
  summarize(totalreads = sum(nreads)) %>%
  arrange(desc(totalreads)) %>% slice_head(n=10) %>% select(rank)

df <- holi_data %>%
  filter(grepl("Eukaryota", taxa_path)) %>%
  filter(rank %in% ranks_to_plot$rank) %>%
  group_by(date, rank) %>%
  summarize(numreads = sum(nreads)) %>%
  mutate(proportion = numreads / sum(numreads))

plot_total <- ggplot(df, aes(x = as.factor(date), y = numreads, fill = rank)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL)

# Create the second plot showing the proportions of reads
plot_proportion <- ggplot(df, aes(x = as.factor(date), y = proportion, fill = rank)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  scale_y_continuous(labels = scales::percent)

combined_plot <- plot_total / plot_proportion +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(plot_directory, "/A.euk.ranks.png"))

# ----------------- CONSRTUCTING THE DAMAGE MODEL  ----------------- #

# ----------- B) good and bad example fits ----------- #
plot_example_fits <- function(dat, goodbad, plotname="", nreads = 1, howmany){
  tax <- dat |>
    ungroup() |>
    filter(fit == goodbad & n_reads > nreads) |>
    group_by(label) |>
    slice_sample(n = howmany) |>
    ungroup()

  samples <- tax$label |> unique()

  plots <- purrr::map(.x = samples, dat = tax, .f = function(x, dat, orient = "fwd", pos = 25, p_breaks = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)) {
    data <- dat |> filter(label == x)
    grid_size <- calculate_plot_grid(length(data$tax_name))
    l <- lapply(data$name, function(X) {
      df1 <- data |>
        filter(name == X)
      p <- get_dmg_decay_fit(df1, orient = orient, pos = pos, p_breaks = p_breaks)
      p <- p + ggtitle(X)
      return(p)
    })
    plot <- ggpubr::ggarrange(plotlist = l, ncol = grid_size$cols, nrow = grid_size$rows, align = "hv")
    tit <- paste0(x, " -- ", data$label_fig |> unique())
    ggpubr::annotate_figure(plot, top = ggpubr::text_grob(tit,
                                                          color = "black", face = "bold", size = 12
    ))
  }, .progress = TRUE)

  names(plots) <- samples
  pdf(plotname, width = 20, height = 20)
  print(plots)
  dev.off()
}


plot_example_fits(dat_filt, "good", paste0(plot_directory, "/B.good_examples.100reads.pdf"), 100, 10)
plot_example_fits(dat_filt, "bad", paste0(plot_directory, "/B.bad_examples.100reads.pdf"), 100, 10)

# ----------- B) good vs bad over time ----------- #
calculate_proportion_good <- function(data, n_reads_threshold) {
  data %>%
    filter(n_reads >= n_reads_threshold) %>%
    group_by(date, PlantAnimal) %>%
    summarise(total = n(),
              good_count = sum(fit == "good")) %>%
    mutate(proportion_good = good_count / total,
           threshold = paste0("n_reads >= ", n_reads_threshold)) %>%
    ungroup()
}

thresholds <- c(0, 20, 100, 500)
df_list <- lapply(thresholds, calculate_proportion_good, data = agg_stats)
df_combined <- bind_rows(df_list)

ggplot(df_combined, aes(x = date, y = proportion_good, color = PlantAnimal)) +
  geom_line() +
  facet_wrap(~ threshold) +
  labs(x = "Date",
       y = "Proportion of Good Fits",
       fill = "Plant/Animal")

ggsave(paste0(plot_directory, "/B.good_vs_bad.png"))


# ----------- C) damage all data ----------- #

ggplot(agg_stats[agg_stats$n_reads>100,], aes(x = median_A_b, y = date, color = fit, size=n_reads)) +
  geom_point(alpha=0.6) + facet_wrap(~PlantAnimal)

ggsave(paste0(plot_directory, "/C.dmg.alldata.100reads.png"))


# ----------- D) damage in "good" plants 500  ----------- #

get_top <- function(df, n, planimal){
  tmp <- df %>%
  filter(PlantAnimal == planimal) %>%
  group_by(genus) %>%
  summarize(total_reads = sum(n_reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = n)
  return(df %>% filter(genus %in% tmp$genus))
}

df <- get_top(agg_stats, 20, "plant") %>% filter(n_reads >= 500)
ggplot(df, aes(x = median_A_b, y = date, color = fit, size=n_reads)) +
  geom_point(alpha=0.8) +
  facet_wrap(~genus)
ggsave(paste0(plot_directory, "/D.dmg.plants.500reads.png"))


# ----------- E) damage model  ----------- #


get_conditional_quantile <- function(df) {
  # get latest 10 dates
  latest_dates <- df %>%
    distinct(date) %>%
    arrange(date) %>%
    slice_tail(n=10)

  qdata <- df %>%
    mutate(is_first_dates = date %in% latest_dates$date) %>%
    group_by(date, is_first_dates) %>%
    filter(n() >= 10) %>%

    summarise(
      QFILT = if (first(is_first_dates)) {
        quantile(median_A_b, probs = 0.05, na.rm = TRUE)
      } else {
        quantile(median_A_b[fit == 'good'], probs = 0.05, na.rm = TRUE)
      }, .groups = 'drop'
    )
  return(qdata)
}

plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & n_reads >= 500)
qdata <- get_conditional_quantile(plants500)
loess_fit <- loess(QFILT ~ date, data = qdata, span = 0.3)
pred <- predict(loess_fit, newdata = qdata, se = TRUE)
qdata <- qdata %>%
  mutate(
    fit_q = pred$fit,
    lwr = pred$fit - 1.96 * pred$se.fit,
    upr = pred$fit + 1.96 * pred$se.fit,
  )

plants500 <- inner_join(plants500, qdata) %>% mutate(alpha = ifelse(fit == "good", 1, 0.2), alpha = ifelse(is_first_dates == TRUE, 1, alpha))

ggplot(plants500) +
  geom_point(aes(y = date, x = median_A_b, color=fit, alpha=alpha)) +
  geom_ribbon(data = qdata, aes(y = date, xmin = lwr, xmax = upr), alpha = 0.2) +
  geom_path(data = qdata, aes(y = date, x = QFILT)) +
  geom_path(data = qdata, color="blue", aes(y = date, x = fit_q)) +
  scale_alpha_identity() +
    labs(y = "Date",
       x = "A_b",
       title = "Damage model")
ggsave(paste0(plot_directory, "/E.dmgmodel.png"))



# ----------------- FILTERING USING THE DAMAGE MODEL ----------------- #


agg_stats <- inner_join(agg_stats, qdata) %>% mutate(status = ifelse(median_A_b < lwr, "fail", "pass"),
  status = ifelse(!is_first_dates & fit == "bad", "fail", status))

filter_occ <- function(df, reads, occ){
  filtered_data <- df %>%
  group_by(genus, label) %>%
  filter(n_reads >= reads & status == "pass") %>%
  group_by(genus) %>%
  filter(n_distinct(label) >= occ) %>%
  ungroup()
  unique(filtered_data$genus)
}

# or we can do agg_stats <- read.table("agg.tsv",h=T)

good_genera <- filter_occ(agg_stats, 100, 2) # at least 2 occurences of 100 reads, pass

agg_stats <- agg_stats %>%
  filter(genus %in% good_genera)


# ----------- F) dmg in filtered data  ----------- #

plot_filtered <- function(qdata, df, plotsave, plotname){

  annotations <- df %>%
    group_by(genus) %>%
    summarise(
      n_reads = sum(n_reads),
      pass_proportion = mean(status == "pass")
    ) %>%
    left_join(contam_genera) %>%
    mutate(reads_blanks = scales::comma(reads_blanks), pass_proportion = scales::percent(pass_proportion), n_reads = scales::comma(n_reads))

  nplots <- length(unique(df$genus))

  rows <- 4
  cols <- 3
  num_pages <- ceiling(nplots / (rows * cols))


  pdf(file = plotsave)

  for (page in 1:num_pages) {
    if (nplots < 0) {
      facet_params <- list(facet_wrap_paginate(~genus, page = page))
    } else {
      facet_params <- list(facet_wrap_paginate(~genus, page = page, nrow = rows, ncol = cols))
    }
    nplots <- nplots - (rows * cols)

    # Creating the plot
    plot <- ggplot(df) +
      geom_ribbon(data = qdata, aes(y = date, xmin = lwr, xmax = upr), alpha = 0.4) +
      geom_path(data = qdata, aes(y = date, x = fit_q)) +
      geom_point(aes(y = date, x = median_A_b, color = status, size = n_reads, alpha = ifelse(status == "pass", 1, 0.4)), show.legend = c(color = TRUE, alpha = FALSE)) +
      scale_alpha_identity() +
      geom_text(
          data = annotations,
              aes(
                  label = paste(
                      "Reads: ", n_reads, 
                      "\nPass: ", pass_proportion
                  ),
                  x = Inf, y = Inf
              ),
              hjust = 1.1, vjust = 1.1,
              inherit.aes = FALSE, size = 3
          ) +
      labs(
        x = "Median Damage (A_b)",
        y = "Date (CE)"
      ) +
       geom_text(
          data = annotations,
              aes(
                  label = paste(

                      ifelse(!is.na(reads_blanks), paste("\nReads in Blanks: ", reads_blanks), "")
                  ),
                  x = Inf, y = min(df$date)
              ),
              hjust = 1.1, vjust = 0,
              inherit.aes = FALSE, size = 3, color="red"
          ) +

      xlim(0,0.5)+  
      facet_params[[1]]

    print(plot)

    }
    dev.off()
}

plant_df <- agg_stats[agg_stats$PlantAnimal == "plant",]
animal_df <- agg_stats[agg_stats$PlantAnimal == "animal",]

plot_filtered(qdata, animal_df, paste0(plot_directory, "/F.damage.filt.animals.pdf"), "Damage in filtered animals")
plot_filtered(qdata, plant_df, paste0(plot_directory, "/F.damage.filt.plants.pdf"), "Damage in filtered plants")

pass_animal_df <- animal_df[animal_df$status == "pass",]
pass_plant_df <- plant_df[plant_df$status == "pass",]


# ----------- G) ani, gini, breadth in filtered data  ----------- #

plot_filtered_var <- function(qdata, df, var, plotsave, plotname){

  nplots <- length(unique(df$genus))

  rows <- 4
  cols <- 3
  num_pages <- ceiling(nplots / (rows * cols))

  pdf(file = plotsave)

  for (page in 1:num_pages) {
    if (nplots < 0) {
      facet_params <- list(facet_wrap_paginate(~genus, page = page))
    } else {
      facet_params <- list(facet_wrap_paginate(~genus, page = page, nrow = rows, ncol = cols))
    }
    nplots <- nplots - (rows * cols)
    # Creating the plot
    plot <- ggplot(df) +
      geom_point(aes(y = date, x = !!sym(var), color = median_A_b, size = n_reads)) +
      labs(
        y = "Date (CE)"
      ) +
      facet_params[[1]]

    print(plot)

    }
    dev.off()
}

plot_filtered_var(qdata, pass_animal_df, "mean_read_ani_median", paste0(plot_directory, "/G.ani.filt.animals.pdf"), "ANI in filtered animals")
plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_breadth_exp_ratio", paste0(plot_directory, "/G.breadth.filt.animals.pdf"), "Breadth in filtered animals")
plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_gini", paste0(plot_directory, "/G.gini.filt.animals.pdf"), "Gini in filtered animals")

plot_filtered_var(qdata, pass_plant_df, "mean_read_ani_median", paste0(plot_directory, "/G.ani.filt.plants.pdf"), "ANI in filtered plants")
plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_breadth_exp_ratio", paste0(plot_directory, "/G.breadth.filt.plants.pdf"), "Breadth in filtered plants")
plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_gini", paste0(plot_directory, "/G.gini.filt.plants.pdf"), "Gini in filtered plants")


# ----------- H) % strat plots  ----------- #

do_strat_percentage <- function(dat){
  aggregated_data <- dat %>%
      group_by(date, genus) %>%
      summarize(n_reads = sum(n_reads, na.rm = TRUE), .groups = 'drop') %>%
      group_by(date) %>%
      mutate(total_reads = sum(n_reads),
             pct_reads = (n_reads / total_reads) * 100) %>%
      ungroup() %>%
      select(date, genus, pct_reads)

  data_wide <- aggregated_data %>%
      pivot_wider(names_from = genus, values_from = pct_reads, values_fill = 0) %>%
      arrange(date)

  data_matrix <- as.matrix(data_wide %>% select(-date))
  rownames(data_matrix) <- data_wide$date

  date_values <- as.numeric(rownames(data_matrix))

  strat.plot(data_matrix, y.rev=FALSE, plot.line=TRUE, plot.poly=FALSE, plot.bar=TRUE,
             lwd.bar=4, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01,
             x.pc.lab=TRUE, x.pc.omit0=TRUE, srt.xlabel=45, las=2,
             exag=TRUE, exag.mult=5, ylabel = "Date (CE)", yvar = date_values)
}


# Save the plot to a PDF
pdf(file = paste0(plot_directory, "/H.strat_plants_top20_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_plant_df, 20, "plant"))
dev.off()

pdf(file = paste0(plot_directory, "/H.strat_plants_top50_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_plant_df, 50, "plant"))
dev.off()


pdf(file = paste0(plot_directory, "/H.strat_animals_top20_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_animal_df, 20, "animal"))
dev.off()


pdf(file = paste0(plot_directory, "/H.strat_animals_top50_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_animal_df, 50, "animal"))
dev.off()


# ----------------- EXCEL OUT ----------------- #

# Function to transform the tibble
transform_tibble <- function(tibble, calculation = c("sum", "proportion")) {
  calculation <- match.arg(calculation)  # Ensure the argument is either "sum" or "proportion"

  tibble <- tibble %>%
    select(date, genus, n_reads) %>%  # Select only the required columns
    group_by(date, genus) %>%         # Group by date and genus
    summarize(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop")  # Sum n_reads

  if (calculation == "proportion") {
    tibble <- tibble %>%
      group_by(date) %>%
      mutate(total_reads = sum(n_reads)) %>%  # Calculate the total reads for each date
      mutate(proportion = n_reads / total_reads) %>%  # Calculate the proportion
      select(-total_reads, -n_reads) %>%  # Remove the total_reads and n_reads columns
      pivot_wider(names_from = genus, values_from = proportion)
  } else {
    tibble <- tibble %>%
      pivot_wider(names_from = genus, values_from = n_reads)
  }

  tibble <- tibble %>%
    arrange(desc(date)) %>%             # Arrange by date in descending order
    select(date, tidyselect::peek_vars()) %>%  # Sort genus columns alphabetically
    select(date, order(names(.)[-1]) + 1)  # Move date column to the front

  return(tibble)
}

# Transform the tibbles
out_a <- transform_tibble(pass_animal_df, "sum")
out_p <- transform_tibble(pass_plant_df, "sum")
p_out_a <- transform_tibble(pass_animal_df, "proportion")
p_out_p <- transform_tibble(pass_plant_df, "proportion")

# Create a new workbook
wb <- createWorkbook()

# Add worksheets for each transformed tibble
addWorksheet(wb, "Animals")
writeData(wb, "Animals", out_a)

addWorksheet(wb, "Plants")
writeData(wb, "Plants", out_p)

addWorksheet(wb, "ProportionAnimals")
writeData(wb, "ProportionAnimals", p_out_a)

addWorksheet(wb, "ProportionPlants")
writeData(wb, "ProportionPlants", p_out_p)

addWorksheet(wb, "GeneraInBlanks")
writeData(wb, "GeneraInBlanks", contam_genera)

# Save the workbook to a file
saveWorkbook(wb, paste0(plot_directory, "/Tjornin.xlsx"), overwrite = TRUE)
