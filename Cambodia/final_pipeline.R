
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



args <- commandArgs(trailingOnly = TRUE)

# Check if the correct number of arguments are provided
if (length(args) < 1) {
  stop("Usage: Rscript script.R <core>")
}

# Assign arguments to variables
core <- args[1]



metadata <- read.csv("cambodiadb.csv")
metadata <- metadata %>% select(field_sample_parent_id, field_sample_id, archive_sample_id, library_id, "Master.Depth..cm.", "ArchiveSampleDepthCalTape")
metadata <- metadata %>% distinct()
metadata$time <- metadata$ArchiveSampleDepthCalTape # change this 

finished <- read.csv("../files/241129.files.csv", header=FALSE)
names(finished) <- c("basedir", "basename")
finished$library_id <- sapply(finished$basename, function(x) strsplit(x,"_")[[1]][2])


big <- inner_join(metadata, finished)

  print(core)

  df <- big %>% filter(field_sample_parent_id == core)
  outdir <- paste0("output/", core)

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



  write_tsv(df, paste0(outdir, "/data.csv"))
  # when this is made, us it as input to get_agg.sh to do the metadmg 
  METADMG_DIR <- paste0(outdir, "/metadmg/")

  PLOT_DIRECTORY <- paste0(outdir, "/plots/")
  dir.create(PLOT_DIRECTORY, showWarnings = FALSE)

  # have we already run this? if so we can save a lot of time by using save taxonomy info
  TAXONOMY_TSV <- paste0(outdir, "/saved_taxids.tsv")

  # save damage model here 
  GLOBAL_DMG_MODEL <- paste0(outdir, "/qdata.tsv")

  # final excel 
  EXCEL_OUT <- paste0(outdir, "/report.xlsx")


  # ----------------- READ NUMS  ----------------- #


  get_total_seqs <- function(file_path) {
    cmd <- paste("unzip -p", shQuote(file_path), "\"*.txt\" | grep \"Total Seq\" | cut -f2")
    total_seqs <- as.numeric(system(cmd, intern = TRUE))
    if (is.na(total_seqs)) {
      return(0)
    }
    return(total_seqs)
  }

  get_lane_seqs <- function(fastqc_files) {
    lanes <- str_extract(fastqc_files, "L00[1-4]")
    lane_seqs <- sapply(unique(lanes), function(lane) {
      files_by_lane <- fastqc_files[lanes == lane]
      sum(sapply(files_by_lane, get_total_seqs), na.rm = TRUE)
    }, simplify = FALSE)
    return(lane_seqs)
  }

  # get number of raw reads per lane 
  lane_df <- df %>%
    rowwise() %>%
    mutate(
      fastqc_files = list(
        as.character(list.files(
          path = paste0(basedir, "/stats/reads/fastqc_raw/"),
          pattern = "*R1_fastqc.zip$",
          full.names = TRUE
        )) # list of relevant fastqc files per lib 
      ),
      lane_seqs = list(get_lane_seqs(fastqc_files))
    ) %>%
    ungroup()

  lane_df <- lane_df %>%
    unnest_longer(lane_seqs) %>% # unnest so we have a new col, lane 
    mutate(lane = names(lane_seqs),  
           total_seqs = unlist(lane_seqs)) %>%  
    select(-lane_seqs, -fastqc_files, -lane_seqs_id) 

  # same for collapsed reads 
  collapsed_lane_df <- df %>%
    rowwise() %>%
    mutate(
      trim_fastqc_files = list(
        as.character(list.files(
          path = paste0(basedir, "/stats/reads/fastqc_trim/"),
          pattern = "*collapsed_fastqc.zip$",
          full.names = TRUE
        ))
      ),
      lane_seqs = list(get_lane_seqs(trim_fastqc_files))
    ) %>%
    ungroup()

  collapsed_lane_df <- collapsed_lane_df %>%
    unnest_longer(lane_seqs) %>%  
    mutate(lane = names(lane_seqs),  
           collapsed_seqs = unlist(lane_seqs)) %>%  
    select(-lane_seqs, -trim_fastqc_files, -lane_seqs_id) 

  # join raw and collapsed 
  lane_df <- inner_join(lane_df, collapsed_lane_df)

  lane_df <- lane_df %>%
    arrange(-time) %>%  
    mutate(label = factor(library_id, levels = unique(library_id))) 

  ggplot(lane_df, aes(y=label, x=total_seqs, color=lane)) +
    geom_point() + 
    ggtitle("Total Sequenced Reads per Lane")
  ggsave(paste0(PLOT_DIRECTORY, "/reads_sequenced_per_lane.png"))

  ggplot(lane_df, aes(y=label, x=collapsed_seqs, color=lane)) +
    geom_point() + 
    ggtitle("Collapsed Reads per Lane")
  ggsave(paste0(PLOT_DIRECTORY, "/reads_collapsed_per_lane.png"))

  # now sum across lanes 
  summed_df <- lane_df %>%
    group_by(label, library_id) %>%
    summarize(
      across(c(basedir,basename, time), first),
      total_seqs = sum(total_seqs),
      collapsed_seqs = sum(collapsed_seqs),
      .groups = "drop"
    ) %>% 
    arrange(-time) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id)))  


  # we can also get the number of assigned readsd from metadmg 
  get_ass_seqs <- function(file_path) {
    cmd <- paste("zcat", shQuote(file_path), "|sed 1d|head -n1|cut -f4")
    total_seqs <- as.numeric(system(cmd, intern = TRUE))
    if (is.na(total_seqs)) {
      return(0)
    }
    return(total_seqs)
  }

  ass_df <- df %>%
    rowwise() %>%
    mutate(
      agg_file = 
        as.character(list.files(
          path = paste0(basedir, "/stats/metadmg/aggregate/"),
          pattern = "*_collapsed.stat.gz$",
          full.names = TRUE
        )
      ),
      ass_seqs = get_ass_seqs(agg_file)
    ) %>%
    ungroup()

  # join with other data
  summed_df <- inner_join(ass_df, summed_df)

  # get % 
  summed_df <- summed_df %>%
    mutate(collapsed_pct = (collapsed_seqs / total_seqs) * 100,
    assigned_pct = (ass_seqs / collapsed_seqs) * 100)

  long_reads <- summed_df %>% 
    pivot_longer(
      cols = c(total_seqs, collapsed_seqs, ass_seqs),
      names_to = "read_type",
      values_to = "reads"
    ) %>%
    arrange(-time) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id))) 


  ggplot(long_reads, aes(y=library_id, x=reads, fill=read_type)) +
    geom_bar(stat="identity", position="dodge") +
    ggtitle("Sequenced, trimmed and assigned reads")
  ggsave(paste0(PLOT_DIRECTORY, "reads_seq_collapsed_ass.png"))

  long_pct <- summed_df %>%
    pivot_longer(
      cols = c(collapsed_pct, assigned_pct),
      names_to = "read_type_pct",
      values_to = "reads_pct"
    ) %>%
    arrange(-time) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id))) 

  ggplot(long_pct, aes(y=library_id, x=reads_pct, color=read_type_pct)) +
    geom_point() +
    ggtitle("Trimmed and assigned reads as percent")
  ggsave(paste0(PLOT_DIRECTORY, "reads_collapsed_ass_perc.png"))

  # ----------------- READ IN DATA  ----------------- #


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

# will this work? we dont want big numbers because of --howmany issue 
print(names(holi_data))

holi_data <- holi_data %>%
  select(
    which(
      colnames(.) %>% 
        map_lgl(~ {
          # Extract numeric part from column name
          num_part <- as.numeric(str_extract(.x, "\\d+"))
          # Return TRUE for columns with NA or numbers <= 15
          is.na(num_part) || num_part <= 15
        })
    )
  )

print(names(holi_data))
  # filter the metadmg results to just euks 
  # adds plant/animal column 
  # adds "good" or "bad" fit
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

  holi_data$label <- holi_data$library_id
  dat_filt <- filter_metadmg(holi_data, holi_data$label |> unique())



  # ----------------- BASIC PLOTS  ----------------- #
  # rl 
  median_fill <- dat_filt %>%
    filter(n_reads >= 10) %>%
    group_by(time) %>%
    summarize(median_n_reads = median(n_reads, na.rm = TRUE), .groups = 'drop')

  # Join the median values back to the filtered dataset
  dat_filt_with_medians <- dat_filt %>%
    filter(n_reads >= 10) %>%
    left_join(median_fill, by = "time")

  # Plot with the updated median fill
  ggplot(dat_filt_with_medians, aes(y = as.factor(time), x = mean_rlen, fill = median_n_reads,color=median_n_reads)) +
    geom_boxplot() +
    scale_y_discrete(limits = rev(levels(as.factor(dat_filt$time)))) + 
    facet_wrap(~PlantAnimal) +
    ggtitle("Mean read length over depth") +
    scale_fill_gradient(name = "Median n_reads") # Optional: to add a color legend

  # Save the plot
  ggsave(paste0(PLOT_DIRECTORY, "/read_lengths.png"))
  # ----------- A) ranks classified  ----------- #

  # get the top 10 most commonly occuring ranks reads are classified at 
  ranks_to_plot <- holi_data %>%
    group_by(rank) %>%
    summarize(totalreads = sum(nreads)) %>%
    arrange(desc(totalreads)) %>% slice_head(n=10) %>% select(rank)

  # this is just based on metadmg data, not the aggregated stats 
  df <- holi_data %>%
    filter(grepl("Eukaryota", taxa_path)) %>% # get just eukaryotes 
    filter(rank %in% ranks_to_plot$rank) %>% # get ranks we already decided to plot
    group_by(library_id, time, rank) %>%
    summarize(numreads = sum(nreads)) %>% # number of reads per rank 
    mutate(proportion = numreads / sum(numreads)) %>% # proportion 
    arrange(-time) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id))) 


  # plot ranks classified (raw numbers)
  plot_total <- ggplot(df, aes(y = as.factor(time), x = numreads, fill = rank)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(levels(as.factor(df$time)))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Number of reads per rank")

  # plot proportion
  plot_proportion <- ggplot(df, aes(y = as.factor(time), x = proportion, fill = rank)) +
    geom_bar(stat = "identity") +
    scale_y_discrete(limits = rev(levels(as.factor(df$time)))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    ggtitle("Proportion of reads per rank")

  # combine 
  combined_plot <- plot_total / plot_proportion +
      plot_layout(guides = "collect") & theme(legend.position = "bottom")
  combined_plot
  ggsave(paste0(PLOT_DIRECTORY, "/euk.ranks.png"))

  # superkingdoms 

  extract_superkingdom <- function(taxa_path) {
    match <- str_extract(taxa_path, "(?<=:)[^;]+(?=:\"superkingdom\")")
    match <- str_replace_all(match, '[\"/]', '')
    return(ifelse(is.na(match), NA, match))
  }

  holi_data$superkingdom <- sapply(holi_data$taxa_path, extract_superkingdom)

  df_prop <- holi_data %>%
    group_by(library_id, time, superkingdom) %>%
    summarise(total_reads = sum(nreads, na.rm = TRUE)) %>%  
    mutate(prop = total_reads / sum(total_reads)) %>% # and also the proportion    
    arrange(-time) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id))) 

  ggplot(df_prop, aes(y = as.factor(time), x = total_reads, fill = superkingdom)) + 
      geom_bar(stat = "identity") + 
      scale_y_discrete(limits = rev(levels(as.factor(df$time)))) + 
      labs(x = "Number of reads", fill = "Superkingdom", title = " Number of reads classified to each superkingdom")
  ggsave(paste0(PLOT_DIRECTORY, "/superkingdoms.png"))


  # ----------------- CONSRTUCTING THE DAMAGE MODEL  ----------------- #
  # now we make the damage model based on damage vs depth 

  # ----------- B) good and bad example fits ----------- #
  plot_example_fits <- function(dat, goodbad, plotname="", nreads = 1, howmany){
    tax <- dat |>
      ungroup() |>
      filter(fit == goodbad & n_reads > nreads) |> # get entries with either good or bad fit and over read limit (depending on func arg) 
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


  plot_example_fits(dat_filt, "good", paste0(PLOT_DIRECTORY, "/good_examples.100reads.pdf"), 100, 100)
  plot_example_fits(dat_filt, "bad", paste0(PLOT_DIRECTORY, "/bad_examples.100reads.pdf"), 100, 100)



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
    group_by(label, genus, fit, PlantAnimal) %>%
    summarise(
      n = n(),  # Count rows
      median_A_b = median(A_b, na.rm = TRUE),  # Median of A_b
      penalized_weighted_median_A_b = penalized_weighted_median(A_b, n_reads, A_b), # Custom function
      median_c_b = median(c_b, na.rm = TRUE),  # Median of c_b
      median_n_reads = median(n_reads, na.rm = TRUE),  # Median of n_reads
      mean_n_reads = mean(n_reads, na.rm = TRUE),      # Mean of n_reads
      total_n_reads = sum(n_reads, na.rm = TRUE)       # Sum of n_reads
    ) %>%
    ungroup() 

  agg_stats$library_id <- agg_stats$label 
  agg_stats <- inner_join(agg_stats,metadata)

  ggplot(agg_stats[agg_stats$total_n_reads>100,], aes(x = median_A_b, y = time, color = fit, size=total_n_reads)) +
    geom_point(alpha=0.6) + facet_wrap(~PlantAnimal) + 
    labs(y = "Depth", title = " Damage over time") + 
    scale_y_reverse()
  ggsave(paste0(PLOT_DIRECTORY, "/dmg.alldata.100reads.png"))


  # ----------- D) damage in "good" plants 500  ----------- #

  # small function to get the n most abundant genera in either plant/animal 
  get_top <- function(df, n, planimal){
    tmp <- df %>%
    filter(PlantAnimal == planimal) %>%
    group_by(genus) %>%
    summarize(total_reads = sum(total_n_reads, na.rm = TRUE)) %>%
    arrange(desc(total_reads)) %>%
    slice_head(n = n)
    return(df %>% filter(genus %in% tmp$genus))
  }

  # plot 20 most abundant plants, with good and bad fit, only those with >200 reads 
  df <- get_top(agg_stats, 20, "plant") %>% filter(total_n_reads >= 200)
  ggplot(df, aes(x = median_A_b, y = time, color = fit, size=total_n_reads)) +
    geom_point(alpha=0.8) +
    facet_wrap(~genus) + 
    labs(y = "Depth", title = " Damage over time in confident plants") + 
    scale_y_reverse()

  ggsave(paste0(PLOT_DIRECTORY, "/dmg.plants.200reads.png"))


  # ----------- E) damage model  ----------- #


  get_conditional_quantile <- function(df) {
    # get 3 first layers
    latest_times <- df %>%
      distinct(time) %>%
      arrange(time) %>%
      slice(c(1:3))


    # find damage quantiles 
    qdata <- df %>%
      mutate(is_first_times = time %in% latest_times$time) %>%
      group_by(time, is_first_times) %>%
        summarise(
          QFILT = if (first(is_first_times)) {
            quantile(median_A_b, probs = 0.05, na.rm = TRUE) # keep all taxa 
          } else {
            quantile(median_A_b[fit == 'good'], probs = 0.05, na.rm = TRUE) # just good fit for the middle times 
          }, .groups = 'drop'
        ) 

    return(qdata)
  }

  # get plants with > 500 reads 
  plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & total_n_reads >= 200)

  # get damage line 
  qdata <- get_conditional_quantile(plants500)
  loess_fit <- loess(QFILT ~ time, data = qdata, span = 0.3) # loess fit to smooth it 
  pred <- predict(loess_fit, newdata = qdata, se = TRUE)
  qdata <- qdata %>%
    mutate(
      fit_q = pred$fit,
      lwr = pred$fit - 1.96 * pred$se.fit,
      upr = pred$fit + 1.96 * pred$se.fit,
    ) # get 5% percentile, use this as the limit 
  write.table(qdata, GLOBAL_DMG_MODEL)

  # assign pass or fail if passed line, if not one of the first or last times it must also be good fit 
  agg_stats <- inner_join(agg_stats, qdata) %>% mutate(status = ifelse(median_A_b < lwr, "fail", "pass"),
    status = ifelse(!is_first_times & fit == "bad", "fail", status))

  # for plotting, assign a low  alpha val to the failed taxa 
  plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & total_n_reads >= 200) %>% 
    mutate(alpha = ifelse(status=="pass",1,0.2))

  # plot damage model and the points that made it 
  ggplot(plants500) +
    geom_point(aes(y = time, x = median_A_b, color=fit, alpha=alpha)) +
    geom_ribbon(data = qdata, aes(y = time, xmin = lwr, xmax = upr), alpha = 0.2) +
    geom_path(data = qdata, aes(y = time, x = QFILT)) +
    geom_path(data = qdata, color="blue", aes(y = time, x = fit_q)) +
    scale_alpha_identity() + scale_y_reverse() +
      labs(y = "time",
         x = "A_b",
         title = "Damage model")
  ggsave(paste0(PLOT_DIRECTORY, "/dmgmodel.png"))



  # ----------------- FILTERING USING THE DAMAGE MODEL ----------------- #


  # find passing genera with at least reads at occ different times 
  filter_occ <- function(df, reads, occ){
    filtered_data <- df %>%
    group_by(genus, label) %>%
    filter(total_n_reads >= reads & status == "pass") %>%
    group_by(genus) %>%
    filter(n_distinct(label) >= occ) %>%
    ungroup()
    unique(filtered_data$genus)
  }



  good_genera <- filter_occ(agg_stats, 50, 1)

  agg_stats <- agg_stats %>%
    filter(genus %in% good_genera)
  agg_stats$n_reads <- agg_stats$total_n_reads

  # ----------- F) dmg in filtered data  ----------- #

  plot_filtered <- function(qdata, df, plotsave, plotname){
    # get annotations about the number of reads in blanks and number of reads in genus
    annotations <- df %>%
      group_by(genus) %>%
      summarise(
        n_reads = sum(n_reads),
        pass_proportion = mean(status == "pass")
      ) %>%
      #left_join(contam_genera) %>% # cant do negatives 
      #mutate(reads_blanks = scales::comma(reads_blanks), pass_proportion = scales::percent(pass_proportion), n_reads = scales::comma(n_reads))
      mutate(pass_proportion = scales::percent(pass_proportion), n_reads = scales::comma(n_reads))


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
      geom_ribbon(data = qdata, aes(y = time, xmin = lwr, xmax = upr), alpha = 0.4) + # damage model area 
      geom_path(data = qdata, aes(y = time, x = fit_q)) + # damage model mean 
      geom_point(aes(y = time, x = median_A_b, color = status, size = n_reads, alpha = ifelse(status == "pass", 1, 0.4)), show.legend = c(color = TRUE, alpha = FALSE)) +
      scale_alpha_identity() +
      geom_text( 
          data = annotations, # how many reads annotation
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
        y = "time "
      ) + scale_y_reverse() +
  #    geom_text(
  #       data = annotations, # reads in blanks annotation 
  #           aes(
  #               label = paste(
  #
  #                   ifelse(!is.na(reads_blanks), paste("\nReads in Blanks: ", reads_blanks), "")
  #               ),
  #               x = Inf, y = min(df$time)
  #           ),
  #           hjust = 1.1, vjust = 0,
  #           inherit.aes = FALSE, size = 3, color="red" # red for angry 
  #       ) +

      xlim(0,0.5)+  
      facet_params[[1]]

    print(plot)

    }
    dev.off()
  }

  # seperate animals and plants 
  plant_df <- agg_stats[agg_stats$PlantAnimal == "plant",]
  animal_df <- agg_stats[agg_stats$PlantAnimal == "animal",]

  plot_filtered(qdata, animal_df, paste0(PLOT_DIRECTORY, "/damage.filt.animals.pdf"), "Damage in filtered animals")
  plot_filtered(qdata, plant_df, paste0(PLOT_DIRECTORY, "/damage.filt.plants.pdf"), "Damage in filtered plants")

  # get just passing genera 
  pass_animal_df <- animal_df[animal_df$status == "pass",]
  pass_plant_df <- plant_df[plant_df$status == "pass",]


  # ----------- G) ani, gini, breadth in filtered data  ----------- #


  # cant do any of this with no bf 


  # plot a given variable 
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
        geom_point(aes(y = time, x = !!sym(var), color = median_A_b, size = n_reads)) +
        labs(
          y = "time (CE)"
        ) +
        facet_params[[1]]

      print(plot)

      }
      dev.off()
  }

  #plot_filtered_var(qdata, pass_animal_df, "mean_read_ani_median", paste0(PLOT_DIRECTORY, "/G.ani.filt.animals.pdf"), "ANI in filtered animals")
  #plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_breadth_exp_ratio", paste0(PLOT_DIRECTORY, "/G.breadth.filt.animals.pdf"), "Breadth in filtered animals")
  #plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_gini", paste0(PLOT_DIRECTORY, "/G.gini.filt.animals.pdf"), "Gini in filtered animals")

  #plot_filtered_var(qdata, pass_plant_df, "mean_read_ani_median", paste0(PLOT_DIRECTORY, "/G.ani.filt.plants.pdf"), "ANI in filtered plants")
  #plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_breadth_exp_ratio", paste0(PLOT_DIRECTORY, "/G.breadth.filt.plants.pdf"), "Breadth in filtered plants")
  #plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_gini", paste0(PLOT_DIRECTORY, "/G.gini.filt.plants.pdf"), "Gini in filtered plants")


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


  # Save the plot to a PDF
  pdf(file = paste0(PLOT_DIRECTORY, "/strat_plants_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_plant_df, 20, "plant"))
  dev.off()

  pdf(file = paste0(PLOT_DIRECTORY, "/strat_plants_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_plant_df, 50, "plant"))
  dev.off()


  pdf(file = paste0(PLOT_DIRECTORY, "/strat_animals_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_animal_df, 20, "animal"))
  dev.off()


  pdf(file = paste0(PLOT_DIRECTORY, "/strat_animals_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_animal_df, 50, "animal"))
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
  out_a <- transform_tibble(pass_animal_df, "sum")
  out_p <- transform_tibble(pass_plant_df, "sum")
  p_out_a <- transform_tibble(pass_animal_df, "proportion")
  p_out_p <- transform_tibble(pass_plant_df, "proportion")

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
