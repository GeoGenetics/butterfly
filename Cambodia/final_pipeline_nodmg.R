
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
library(optparse)
library(yaml)

# internal functions
source("../dmg.R")
source("../get_calculate_plot_grid.R")
source("../perk.R")
source("../damage_est_function.R")
source("../perk_wrapper.R")
source("../perk_wrapper_function.R")
source("../get_dmg_decay_fit.R")
source("../median.R")
source("../filter.R")

# ----------------- ARGS  ----------------- #

# example run
#    Rscript final_pipeline.R --metadata /projects/caeg/data/pmp/Cambodia/CAM2311M/data.csv --outdir /projects/caeg/data/pmp/Cambodia/CAM2311M/ \





metadata <- read.csv("cambodiadb.csv")
metadata <- metadata %>% select(field_sample_parent_id, field_sample_id, archive_sample_id, library_id, "Master.Depth..cm.", "ArchiveSampleDepthCalTape")
metadata <- metadata %>% distinct()
metadata$time <- metadata$ArchiveSampleDepthCalTape # change this 

finished <- read.csv("../files/241129.files.csv", header=FALSE)
names(finished) <- c("basedir", "basename")
finished$library_id <- sapply(finished$basename, function(x) strsplit(x,"_")[[1]][2])


big <- inner_join(metadata, finished)

for (core in unique(big$field_sample_parent_id)) {
  print(core)
  df <- big %>% filter(field_sample_parent_id == core)
  outdir <- paste0("output/nodmgfilter/", core)
  print(names(df))
  df <- df %>%
    rowwise() %>%
    mutate(
      aggregate_stat_files = 
        as.character(list.files(
          path = paste0(basedir, "/stats/metadmg/aggregate/"),
          pattern = "*_collapsed.stat.gz$",
          full.names = TRUE
        )
      )
  )

  dir.create(outdir, showWarnings = FALSE)

  write_tsv(df, paste0(outdir, "/data.csv"))
  # when this is made, us it as input to get_agg.sh to do the metadmg 
  METADMG_DIR <- paste0("output/", core, "/metadmg/")
  PLOT_DIRECTORY <- paste0(outdir, "/plots/")
  dir.create(PLOT_DIRECTORY, showWarnings = FALSE)

  # have we already run this? if so we can save a lot of time by using save taxonomy info
  TAXONOMY_TSV <- paste0(outdir, "/saved_taxids.tsv")

  # save damage model here 
  GLOBAL_DMG_MODEL <- paste0(outdir, "/qdata.tsv")

  # final excel 
  EXCEL_OUT <- paste0(outdir, "/report.xlsx")


  # ----------------- READ NUMS  ----------------- #


  # ------------- HELPER TO READ LOTS OF FILES  ------------- #

  read_file <- function(f) {
    df <- fread(f, header=T, sep="\t", fill=T, nThread=20)
    df$sample_id <- f
    return(df)
  }

  # ------------- METADMG  ------------- #

  file_list <- list.files(path = METADMG_DIR, pattern = "*gz", recursive = FALSE, full.name=TRUE)

  holi_data <- do.call(rbind, lapply(file_list, read_file))
  holi_data <- holi_data %>%
    mutate(library_id = sapply(strsplit(basename(sample_id), "_"), `[`, 2))


  # merge metadmg data with metadata 
  holi_data <- inner_join(holi_data, df)


  # filter the metadmg results to just euks 
  # adds plant/animal column 
  filter_metadmg <- function(df, samples){
    holi_data_sp_euk <- df |>
      filter(rank == "species") |> # must be classified at species level 
      filter(grepl("Eukaryota", taxa_path)) |> # must be a Eukaryote 
      mutate(PlantAnimal = case_when( # define plant/animal based on taxpath 
        grepl("Viridiplantae", taxa_path) ~ "plant",
        grepl("Metazoa", taxa_path) ~ "animal",
      )) |>
      rename(tax_name = taxid, n_reads = nreads) # rename these for ease later

  }

  holi_data$label <- holi_data$library_id
  dat_filt <- filter_metadmg(holi_data, holi_data$label |> unique())



  # ----------- C) damage all data ----------- #


  dat_filt <- dat_filt %>%
  filter(grepl(':\"[^\"]+\":\"genus\"', taxa_path)) %>%
    mutate(genus = gsub('.*:(\\"[^"]+\\"):\"genus\".*', '\\1', taxa_path) %>% gsub('\"', '', .))



  agg_stats <- dat_filt %>%
    filter(!is.na(genus) & !is.na(PlantAnimal)) %>%  # Remove rows with missing genus or PlantAnimal
    mutate(
      flt = paste(label, name, sep = "--"),  # Combine label and name
      num_alns = sum(nalign, na.rm = TRUE),  # Total number of alignments
      weight = nalign / num_alns            # Weight based on nalign
    ) %>%
    filter(nalign >=1 & n_reads >= 1) %>% 
    group_by(label, genus, PlantAnimal) %>%
    summarise(
      n = n(),  # Count rows
      median_A_b = median(A_b, na.rm = TRUE),  # Median of A_b
      penalized_weighted_median_A_b = penalized_weighted_median(A_b, n_reads, A_b), # Custom function
      median_c_b = median(c_b, na.rm = TRUE),  # Median of c_b
      median_n_reads = median(n_reads, na.rm = TRUE),  # Median of n_reads
      mean_n_reads = mean(n_reads, na.rm = TRUE),      # Mean of n_reads
      total_n_reads = sum(n_reads, na.rm = TRUE),       # Sum of n_reads
      median_readlen = median(mean_rlen, na.rm = TRUE)       # Sum of n_reads

    ) %>%
    ungroup() 

  agg_stats$library_id <- agg_stats$label 
  agg_stats <- inner_join(agg_stats,metadata)

  ggplot(agg_stats[agg_stats$total_n_reads>100,], aes(x = median_A_b, y = time, color = median_readlen, size=total_n_reads)) +
    geom_point(alpha=0.6) + facet_wrap(~PlantAnimal) + 
    labs(y = "Depth", title = " Damage over time") + 
    scale_y_reverse()
  ggsave(paste0(PLOT_DIRECTORY, "/dmg.alldata.100reads.png"))



  filter_occ_nopass <- function(df, reads, occ){
    filtered_data <- df %>%
    group_by(genus, label) %>%
    filter(total_n_reads >= reads) %>%
    group_by(genus) %>%
    filter(n_distinct(label) >= occ) %>%
    ungroup()
    unique(filtered_data$genus)
  }

  good_genera <- filter_occ_nopass(agg_stats, 50, 1)

  agg_stats <- agg_stats %>%
    filter(genus %in% good_genera)
  agg_stats$n_reads <- agg_stats$total_n_reads

  # ----------- F) dmg in filtered data  ----------- #

  plot_filtered <- function(df, plotsave, plotname){
    # get annotations about the number of reads in blanks and number of reads in genus
    annotations <- df %>%
      group_by(genus) %>%
      summarise(
        n_reads = sum(n_reads),
      ) 


    # find out how many pages we need 
    nplots <- length(unique(df$genus))
    rows <- 4
    cols <- 3
    num_pages <- ceiling(nplots / (rows * cols))


    pdf(file = plotsave)

    # loop through each page 
    for (page in 1:num_pages) {
      if (nplots < 0) {
        facet_params <- list(facet_wrap_paginate(~genus, page = page)) # if its the last page, we dont care how many rows etc 
      } else {
        facet_params <- list(facet_wrap_paginate(~genus, page = page, nrow = rows, ncol = cols))
      }
    nplots <- nplots - (rows * cols)

    # creating the plot
    plot <- ggplot(df) +
      geom_point(aes(y = time, x = median_A_b, color = median_readlen, size = n_reads)) +
      scale_alpha_identity() +
      geom_text( 
          data = annotations, # how many reads annotation
              aes(
                  label = paste(
                      "Reads: ", n_reads
                  ),
                  x = Inf, y = Inf
              ),
              hjust = 1.1, vjust = 1.1,
              inherit.aes = FALSE, size = 3
          ) +
      labs(
        x = "Median Damage (A_b)",
        y = "Depth "
      ) + scale_y_reverse() +

      xlim(0,0.5)+  
      facet_params[[1]]

    print(plot)

    }
    dev.off()
  }

  # seperate animals and plants 
  plant_df <- agg_stats %>% filter(PlantAnimal == "plant" & n_reads >= 50)
  animal_df <-  agg_stats %>% filter(PlantAnimal == "animal" & n_reads >= 50)

  plot_filtered(animal_df, paste0(PLOT_DIRECTORY, "/damage.filt.animals.pdf"), "Damage in filtered animals")
  plot_filtered(plant_df, paste0(PLOT_DIRECTORY, "/damage.filt.plants.pdf"), "Damage in filtered plants")


  # ----------- G) ani, gini, breadth in filtered data  ----------- #


  # cant do any of this with no bf 



  # ----------- H) % strat plots  ----------- #

  do_strat_percentage <- function(dat){
    aggregated_data <- dat %>%
        group_by(time, genus) %>%
        summarize(n_reads = sum(n_reads, na.rm = TRUE), .groups = 'drop') %>%
        group_by(time) %>%
        mutate(total_reads = sum(n_reads),
               pct_reads = (n_reads / total_reads) * 100) %>%
        ungroup() %>%
        select(time, genus, pct_reads)

    data_wide <- aggregated_data %>%
        pivot_wider(names_from = genus, values_from = pct_reads, values_fill = 0) %>%
        arrange(time)

    data_matrix <- as.matrix(data_wide %>% select(-time))
    rownames(data_matrix) <- data_wide$time

    time_values <- as.numeric(rownames(data_matrix))

    strat.plot(data_matrix, y.rev=TRUE, plot.line=TRUE, plot.poly=FALSE, plot.bar=TRUE,
               lwd.bar=4, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01,
               x.pc.lab=TRUE, x.pc.omit0=TRUE, srt.xlabel=45, las=2,
               exag=TRUE, exag.mult=5, ylabel = "time", yvar = time_values)
  }

  get_top <- function(df, n, planimal){
    tmp <- df %>%
    filter(PlantAnimal == planimal) %>%
    group_by(genus) %>%
    summarize(total_reads = sum(total_n_reads, na.rm = TRUE)) %>%
    arrange(desc(total_reads)) %>%
    slice_head(n = n)
    return(df %>% filter(genus %in% tmp$genus))
  }
  # Save the plot to a PDF
  pdf(file = paste0(PLOT_DIRECTORY, "/strat_plants_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(plant_df, 20, "plant"))
  dev.off()

  pdf(file = paste0(PLOT_DIRECTORY, "/strat_plants_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(plant_df, 50, "plant"))
  dev.off()


  pdf(file = paste0(PLOT_DIRECTORY, "/strat_animals_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(animal_df, 20, "animal"))
  dev.off()


  pdf(file = paste0(PLOT_DIRECTORY, "/strat_animals_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(animal_df, 50, "animal"))
  dev.off()


  # ----------------- EXCEL OUT ----------------- #

  # sum or get proportion of reads, so genera are columns and times are rows 
  transform_tibble <- function(tibble, calculation = c("sum", "proportion")) {
    calculation <- match.arg(calculation)  

    tibble <- tibble %>%
      select(time, genus, n_reads) %>%  
      group_by(time, genus) %>%        
      summarize(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop")  

    if (calculation == "proportion") {
      tibble <- tibble %>%
        group_by(time) %>%
        mutate(total_reads = sum(n_reads)) %>%  
        mutate(proportion = n_reads / total_reads) %>%  
        select(-total_reads, -n_reads) %>%  
        pivot_wider(names_from = genus, values_from = proportion)
    } else {
      tibble <- tibble %>%
        pivot_wider(names_from = genus, values_from = n_reads)
    }

    tibble <- tibble %>%
      arrange((time)) %>%             
      select(time, tidyselect::peek_vars()) %>%  
      select(time, order(names(.)[-1]) + 1)  

    return(tibble)
  }

  # transform data 
  out_a <- transform_tibble(animal_df, "sum")
  out_p <- transform_tibble(plant_df, "sum")
  p_out_a <- transform_tibble(animal_df, "proportion")
  p_out_p <- transform_tibble(plant_df, "proportion")

  # create a new workbook
  wb <- createWorkbook()

  # add worksheets
  addWorksheet(wb, "Animals")
  writeData(wb, "Animals", out_a)

  addWorksheet(wb, "Plants")
  writeData(wb, "Plants", out_p)

  addWorksheet(wb, "ProportionAnimals")
  writeData(wb, "ProportionAnimals", p_out_a)

  addWorksheet(wb, "ProportionPlants")
  writeData(wb, "ProportionPlants", p_out_p)

  #addWorksheet(wb, "GeneraInBlanks")
  #writeData(wb, "GeneraInBlanks", contam_genera)

  # save the workbook to a file
  saveWorkbook(wb, EXCEL_OUT, overwrite = TRUE)

}