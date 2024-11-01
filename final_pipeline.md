
---
title: "Analysis of Tjornen Data"
author: "Abigail Ramsoe"
output: html_document
---

## Dependencies 

Load necessary internal and external libraries 
```{r}
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

# internal functions
source("dmg.R")
source("get_calculate_plot_grid.R")
source("perk.R")
source("damage_est_function.R")
source("perk_wrapper.R")
source("perk_wrapper_function.R")
source("get_dmg_decay_fit.R")
source("median.R")
source("filter.R")
```


## Filepaths 

We need names, nodes, and mapping of accession to taxid, these don't change unless the reference collection changes 
```{r names-nodes}
# probably don't need to change 
NAMES <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
NODES <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
ACC2TAXID <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"
```

Path for result directories (metadmg aggregate and bamfilter stats)
```{r result-dirs}
# change per project
METADMG_DIR <- "/maps/projects/caeg/scratch/iceland/raw_results/stats/metadmg/aggregate/"
BAMFILTER_DIR <- "/maps/projects/caeg/scratch/iceland/raw_results/results/bamfilter/"

# negative data 
NEG_METADMG_DIR <- "/maps/projects/caeg/scratch/iceland/raw_results/stats/metadmg/aggregate/lib/"
NEG_BAMFILTER_DIR <- "/maps/projects/caeg/scratch/iceland/raw_results/results/bamfilter/lib"
```

Metadata tsv file, with sample IDs and dates 
```{r metadata-file}
METADATA <- "/maps/projects/caeg/scratch/iceland/new_metadata.tsv"
```

Places to save output 
```{r output-paths}
PLOT_DIRECTORY <- "/maps/projects/caeg/scratch/iceland/finalplotsV2/"

# have we already run this? if so we can save a lot of time by using save taxonomy info
TAXONOMY_TSV <- "/maps/projects/caeg/scratch/iceland/saved_taxids.tsv" 

# save agg stats here so we can skip this if needed
AGG_TSV <- "/maps/projects/caeg/scratch/iceland/agg_stats.tsv"
NEGAGG_TSV <- "/maps/projects/caeg/scratch/iceland/neg-agg_stats.tsv"

# save damage model here 
GLOBAL_DMG_MODEL <- "/maps/projects/caeg/scratch/iceland/qdata.tsv"

# final excel 
EXCEL_OUT <- "TJORNIN.XLSX" # this is just a basename, it is saved in PLOT_DIRECTORY
```

## Read in all data 

Read in the metadata table, remove the "_" from the CGG number, and save just the "date" and "label" columns, this is all we care about in this case. 

```{r load-metadata}
metadata <- read.table(METADATA,he=T,as.is=T)
metadata$label <- gsub("_", "", metadata$CGG_ID)
names(metadata)[2] <- "date" # we know this column should be called date, this must be changed per project 
metadata <- metadata[, c("label","date")]

head(metadata)
```

We need a small helper function to read large files and add the file name itself to a column, we will use this a lot 
```{r fun_read-file}
read_file <- function(f) {
  df <- fread(f, header=T, sep="\t", fill=T, nThread=20)
  df$sample_id <- fpac
  return(df)
}
```

### Metadmg 

Load in metadmg data to a big dataframe. Firstly we get a list of all the metadmg files, and then read them all in to one big dataframe. Lastly, we extract the CGG number from the filepath. 
```{r load-metadmg}
file_list <- list.files(path = METADMG_DIR, pattern = "*gz", recursive = FALSE, full.name=TRUE)
holi_data <- do.call(rbind, lapply(file_list, read_file))
holi_data$label <- sub(".*(CGG[0-9]+)\\..*", "\\1", holi_data$sample_id) # this will have to change when we have non-CGG ids 

head(holi_data)
```

Now we merge the metadmg data with the metadata 
```{r merge-holi-metadata}
holi_data <- inner_join(holi_data, metadata, by = "label")
```

We need a little function to filter the metadmg data. We are only interested in taxa that are classified at species level, and that are Eukaryotes. We make a new column that shows if a taxa is a plant or an animal (based on the taxon path). We then find if the damage profile for each taxa is "good" or "bad" fit to the expected damage model. 
```{r merge-holi-metadata}
filter_metadmg <- function(df, samples){
  holi_data_sp_euk <- df |>
    filter(rank == "species") |> # must be classified at species level 
    filter(grepl("Eukaryota", taxa_path)) |> # must be a Eukaryote 
    mutate(PlantAnimal = case_when( # define plant/animal based on taxpath 
      grepl("Viridiplantae", taxa_path) ~ "plant",
      grepl("Metazoa", taxa_path) ~ "animal",
    )) |>
    rename(tax_name = taxid, n_reads = nreads) # rename these for ease later

  # get the damage fits using CCC
  dat_all <- dmg_fwd_CCC(holi_data_sp_euk, samples, ci = "asymptotic", nperm = 100, nproc = 24)

  # define good and bad hits
  # bad if confidence interval is above 1 or below 0 (despite other params)
  # good if rho_c > 0.85 and C_b > 0.9 and rho_c p-value < 0.1
  dat_filt <- inner_join(dat_all, holi_data_sp_euk) |>
    mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) |>
    mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))
}
```

Apply this filter on our data 
```{r make-dat_filt}
dat_filt <- filter_metadmg(holi_data, metadata$label |> unique())
```