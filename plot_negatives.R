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

metadmg_data <- "negs/NEGS.Tjornen.metadmg_data.tsv"
bamfilter_data <- "negs/NEGS.Tjornen.bam-filter_stats.tsv"

taxonomy_out <- "results/taxonomy/holi-fb-taxids.tsv.gz"
taxonomy_in <- "/projects/lundbeck/scratch/for_abigail/holi-fb-taxids.tsv.gz"

names <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
nodes <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
acc2taxid_file <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"

used_saved_taxid = TRUE



# -----------------  METADMG  ----------------- #

holi_data <- fread(metadmg_data, header=T, sep="\t", fill=T, nThread=20)
names(holi_data)[1] <- "sample_id"
holi_data$label <- sapply(holi_data$sample_id, function(x) strsplit(x,"_")[[1]][3])

holi_data_sp_euk <- holi_data |>
  filter(rank == "species") |>
  filter(grepl("Eukaryota", taxa_path)) |>
  mutate(PlantAnimal = case_when(
    grepl("Viridiplantae", taxa_path) ~ "plant",
    grepl("Metazoa", taxa_path) ~ "animal",
  )) |>
  rename(tax_name = taxid, n_reads = nreads)

# Let's get the damage fits using CCC
samples <- holi_data$label |> unique()
dat_all <- dmg_fwd_CCC(holi_data_sp_euk, samples, ci = "asymptotic", nperm = 100, nproc = 24)

# Define good and bad hits
dat_filt <- inner_join(dat_all, holi_data_sp_euk) |>
  mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
  mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))



# ----------------- BAMFILTER & TAXONOMY ----------------- #

read.names.sql(names, sqlFile = "neg.nameNode.sqlite", overwrite=TRUE)
read.nodes.sql(nodes, sqlFile = "neg.nameNode.sqlite", overwrite=TRUE)

# Load in parsed bamfilter data
fb_data <- fread(bamfilter_data, header=T, sep="\t", fill=T, nThread=20)
names(fb_data)[1] <- "sample_id"
fb_data$label <- sapply(fb_data$sample_id, function(x) strsplit(x,"_")[[1]][3])

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
tax_data <- getTaxonomy(tax_ids_lst, "neg.nameNode.sqlite") |>
  as_tibble() |>
  mutate(taxid = tax_ids_lst)


# ----------------- THE BIG FILTER & STATS MACHINE ----------------- #

# ANI must be 95+

# Adds a column "rm"
# Keep if penalized_weighted_median_breadth_exp_ratio > 0.8 and either penalized_weighted_median_gini < 0.6 or penalized_weighted_median_entropy > 0.75
# Remove taxa without a genus

# Join bamfilter data and metadmg data together with taxonomy info
# Get this at genus level



agg_stats <- get_stats()


