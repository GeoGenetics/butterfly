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
source("getReads.R")

# ----------------- METADATA  ----------------- #

outdir <- "output/"
PLOT_DIRECTORY <- paste0(outdir, "/plots/")
dir.create(PLOT_DIRECTORY, showWarnings = FALSE)

metadata <- fread("250113.cambodia.megatable.tsv")

metadata <- metadata %>%
  select(field_sample_parent_id, archive_sample_id, library_id, archive_sample_master_depth, archive_sample_depth_cal_tape) %>%
  distinct() %>%
  rename(core = field_sample_parent_id, depth = archive_sample_master_depth) %>%
  mutate(
    # Assign site based on core values
    site = case_when(
      core == "CAM2313M" ~ "Reflecting",
      core == "CAM2311M" ~ "Reflecting",
      core == "CAM2201M" ~ "WestBaray",
      core == "CAM2208M" ~ "AngkorThom",
      TRUE ~ NA_character_
    ),
    depth = if_else(is.na(depth) | depth == "", archive_sample_depth_cal_tape, depth) # this has to be triple checked and not used blindly
  )

finished <- read.csv("250113.files.csv", header=FALSE)
names(finished) <- c("basedir", "basename")
finished$library_id <- sapply(finished$basename, function(x) strsplit(x,"_")[[1]][2])


df <- inner_join(metadata, finished) %>% distinct()

ages <- fread("age_depth_model.tsv") %>%
  filter(!is.na(depth))

ages_interpolated <- ages %>%
  arrange(site, depth) %>% 
  group_by(site) %>% 
  mutate(
    depth_next = lead(depth),               
    interpolated_depth = (depth + depth_next) / 2, 
    interpolated_min = (min + lead(min)) / 2,      
    interpolated_max = (max + lead(max)) / 2,      
    interpolated_median = (median + lead(median)) / 2, 
    interpolated_mean = (mean + lead(mean)) / 2    
  ) %>%
  filter(!is.na(depth_next)) %>%                  # ok in this case
  select(site, interpolated_depth, interpolated_min, interpolated_max, interpolated_median, interpolated_mean) %>%
  rename(
    depth = interpolated_depth,
    min = interpolated_min,
    max = interpolated_max,
    median = interpolated_median,
    mean = interpolated_mean
  ) %>%
  bind_rows(ages, .) %>% 
  select(-min, -max, -mean) %>% 
  rename(age_median = median) %>%                           
  arrange(site, age_median)   

df <- inner_join(ages_interpolated, df)  


# ----------------- PLOTTING STUFF  ----------------- #

min_age <- round(min(df$age_median), -2)
max_age <- round(max(df$age_median), -2)
centuries <- seq(min_age, max_age, +100)

core_cols <- c(wes_palette("IsleofDogs1"))

events <- fread("historical_events.tsv")
events$line_label <- as.character(1:nrow(events))  # Create numeric labels 1-7

# ----------------- READ NUMS  ----------------- #

reads_df <- get_reads_df(df)

long_reads <- reads_df %>% 
  pivot_longer(
    cols = c(total_seqs, collapsed_seqs, assigned_seqs),
    names_to = "read_type",
    values_to = "reads"
  ) %>%
  arrange(age_median) %>%  
  mutate(library_id = factor(library_id, levels = unique(library_id))) 

plot <- ggplot(long_reads, aes(x = reads, y = as.factor(age_median), fill = read_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  ggtitle("Sequenced, trimmed and assigned reads") + 
  scale_fill_manual(values=wes_palette("Darjeeling1")) +
  facet_wrap(~site)  +  
  labs(
      x = "Number of Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_seq_collapsed_assigned.png"),plot=plot)


# Raw read numbers 
seq <- ggplot(reads_df, aes(y = age_median, x = total_seqs, color = site)) +
  geom_point(size=5,alpha=0.8) + 
  ggtitle("Sequenced reads") + 
  scale_color_manual(values=core_cols) +
  labs(
      x = "Number of Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_sequenced.png"),plot=seq)


coll <- ggplot(reads_df, aes(y = age_median, x = collapsed_seqs, color = site, size = total_seqs)) +
  geom_point(alpha=0.8) + 
  ggtitle("Collapsed reads") + 
  scale_color_manual(values=core_cols) +
  labs(
      x = "Number of Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_collapsed.png"),plot=coll)


ass <- ggplot(reads_df, aes(y = age_median, x = assigned_seqs, color = site, size = total_seqs)) +
  geom_point(alpha=0.8) + 
  ggtitle("Assigned reads") + 
  scale_color_manual(values=core_cols) +
  labs(
      x = "Number of Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_assigned.png"),plot=ass)

plot <- seq + coll + ass +  plot_layout(guides = "collect")
ggsave(paste0(PLOT_DIRECTORY, "reads_overview.png"),plot=plot)


# Percent read numbers 


coll_pct <- ggplot(reads_df, aes(y = age_median, x = collapsed_pct, color = site, size = total_seqs)) +
  geom_point(alpha=0.8) + 
  ggtitle("Collapsed reads") + 
  scale_color_manual(values=core_cols) +
  labs(
      x = "% Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_collapsed_pct.png"),plot=coll_pct)


ass_pct <- ggplot(reads_df, aes(y = age_median, x = assigned_pct, color = site, size = total_seqs)) +
  geom_point(alpha=0.8) + 
  ggtitle("Assigned reads") + 
  scale_color_manual(values=core_cols) +
  labs(
      x = "% Reads",
      y = "Date",
      )
ggsave(paste0(PLOT_DIRECTORY, "reads_assigned_pct.png"),plot=ass_pct)

plot <- coll_pct + ass_pct +  plot_layout(guides = "collect")
ggsave(paste0(PLOT_DIRECTORY, "reads_overview_pct.png"),plot=plot)


# ----------------- READ IN DATA  ----------------- #


# ------------- HELPER TO READ LOTS OF FILES  ------------- #

read_file <- function(f) {
  df <- fread(f, header=T, sep="\t", fill=T, nThread=20)
  df$sample_id <- f
  return(df)
}

# ------------- METADMG  ------------- #

file_list <- list.files(
  path = "metadmg/raw/",
  pattern = "\\.gz$",  
  recursive = TRUE,
  full.names = TRUE
)

holi_data <- do.call(rbind, lapply(file_list, read_file))
holi_data <- holi_data %>%
  mutate(library_id = sapply(strsplit(basename(sample_id), "_"), `[`, 2),)


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

holi_data <- inner_join(holi_data, reads_df)
dat_filt <- filter_metadmg(holi_data, holi_data$library_id |> unique())


# ------------- READ LENGTH  ------------- #


plot <- ggplot(dat_filt %>% filter(n_reads >= 50),aes(y = as.factor(age_median), x = mean_rlen, fill = assigned_seqs, color = assigned_seqs))+
  geom_boxplot() +
  labs(
      x = "Read length",
      y = "Date",
  ) + 
  facet_wrap(~site)
ggsave(paste0(PLOT_DIRECTORY, "/read_length.pdf"), plot = plot)


# ------------- SUPERKINGDOM  ------------- #

extract_superkingdom <- function(taxa_path) {
  match <- str_extract(taxa_path, "(?<=:)[^;]+(?=:\"superkingdom\")")
  match <- str_replace_all(match, '[\"/]', '')
  return(ifelse(is.na(match), NA, match))
}

species <- holi_data %>% filter(rank == "species")
species$superkingdom <- sapply(species$taxa_path, extract_superkingdom)

df_prop <- species %>% filter(!is.na(superkingdom)) %>%
  group_by(library_id, age_median, superkingdom, site) %>%
  summarise(total_reads = sum(nreads, na.rm = TRUE)) %>%  
  group_by(age_median, site) %>%
  mutate(prop = total_reads / sum(total_reads, na.rm = TRUE)) %>%
  ungroup()

# bars 
plot <- ggplot(df_prop, aes(y = as.factor(age_median), x = total_reads, fill = superkingdom)) +
    geom_bar(stat="identity") +
    labs(
      x = "Number of Reads",
      y = "Date",
      fill = "Superkingdom"
    ) +
    facet_wrap(~site) + 
    theme(axis.text.y = element_text(size = 5))
ggsave(paste0(PLOT_DIRECTORY, "/superkingdom.png"), plot = plot)

# area plots 
plot <- ggplot(df_prop, aes(x = age_median, y = prop, fill = superkingdom)) +
  geom_area() +
  labs(
    x = "Date",
    y = "Proportion",
    fill = "Superkingdom"
  ) +
  scale_x_continuous(breaks = centuries) +
  facet_wrap(~site) +
  geom_vline(
    data = events,
    aes(xintercept = age_median),
    linetype = "dotted",
    color = "black"
  ) +
geom_text(
  data = events,
  aes(x = age_median, label = event, y = 0),
  angle = 90, 
  hjust = -0.1,
  vjust = 0.5,
  color = "black",
  size = 6,
  inherit.aes = FALSE
)
ggsave(paste0(PLOT_DIRECTORY, "/superkingdom_proportions.png"), plot = plot)


# ------------- PLANT ANIMAL  ------------- #

df_prop <- dat_filt %>% filter(!is.na(PlantAnimal)) %>%
  group_by(library_id, age_median, PlantAnimal, site) %>%
  summarise(total_reads = sum(n_reads, na.rm = TRUE)) %>%  
  group_by(age_median, site) %>%
  mutate(prop = total_reads / sum(total_reads, na.rm = TRUE)) %>%
  ungroup()


plot <- ggplot(df_prop, aes(y = as.factor(age_median), x = total_reads, fill = PlantAnimal)) +
    geom_bar(stat="identity") +
    labs(
      x = "Number of Reads",
      y = "Date",
      fill = "Plant/Animal"
    ) +
    facet_wrap(~site) + 
    theme(axis.text.y = element_text(size = 5))
ggsave(paste0(PLOT_DIRECTORY, "/plantanimal.png"), plot = plot)

plot <- ggplot(df_prop, aes(x = age_median, y = prop, fill = PlantAnimal)) +
  geom_area() +
  labs(
    x = "Date",
    y = "Proportion",
    fill = "Plant/Animal"
  ) +
  scale_x_continuous(breaks = centuries) +
  facet_wrap(~site) +
  geom_vline(
    data = events,
    aes(xintercept = age_median),
    linetype = "dotted",
    color = "black"
  ) +
geom_text(
  data = events,
  aes(x = age_median, label = event, y = 0),
  angle = 90, 
  hjust = -0.1,
  vjust = 0.5,
  color = "black",
  size = 6,
  inherit.aes = FALSE
)
ggsave(paste0(PLOT_DIRECTORY, "/plantanimal_proportions.png"), plot = plot)


# -------------  AGGREGATE TO GENUS LEVEL  ------------- #

dat_filt <- dat_filt %>%
filter(grepl(':\"[^\"]+\":\"genus\"', taxa_path)) %>%
  mutate(genus = gsub('.*:(\\"[^"]+\\"):\"genus\".*', '\\1', taxa_path) %>% gsub('\"', '', .))


agg_stats <- dat_filt %>%
  filter(!is.na(genus) & !is.na(PlantAnimal)) %>%  
  mutate(
    flt = paste(label, name, sep = "--"), 
    num_alns = sum(nalign, na.rm = TRUE), 
    weight = nalign / num_alns            
  ) %>%
  filter(nalign >=1 & n_reads >= 1) %>% 
  group_by(site, library_id, genus, PlantAnimal, age_median) %>%
  summarise(
    n = n(),  # Count rows
    median_A_b = median(A_b, na.rm = TRUE),  
    penalized_weighted_median_A_b = penalized_weighted_median(A_b, n_reads, A_b), 
    median_c_b = median(c_b, na.rm = TRUE), 
    median_n_reads = median(n_reads, na.rm = TRUE),  
    mean_n_reads = mean(n_reads, na.rm = TRUE),      
    total_n_reads = sum(n_reads, na.rm = TRUE),       
    median_readlen = median(mean_rlen, na.rm = TRUE)  

  ) %>%
  ungroup() 


# -------------  DAMAGE  ------------- #

get_top <- function(df, n){
  tmp <- df %>%
  group_by(site,genus) %>%
  summarize(total_reads = sum(total_n_reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = n)
  return(df %>% filter(genus %in% tmp$genus))
}

get_conditional_quantile <- function(df) {
  latest_times <- df %>% filter(age_median>2000)
   
  # find damage quantiles 
  qdata <- df %>%
    group_by(age_median) %>%
      summarise(
        QFILT = 
          quantile(median_A_b, probs = 0.05, na.rm = TRUE))
         

  return(qdata)
}

# get plants with > 500 reads 
plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & total_n_reads >= 500)

# get damage line 
qdata <- get_conditional_quantile(plants500)
loess_fit <- loess(QFILT ~ age_median, data = qdata, span = 0.3) # loess fit to smooth it 
pred <- predict(loess_fit, newdata = qdata, se = TRUE)
qdata <- qdata %>%
  mutate(
    fit_q = pred$fit,
    lwr = pred$fit - 1.96 * pred$se.fit,
    upr = pred$fit + 1.96 * pred$se.fit,
  ) # get 5% percentile, use this as the limit 

plot <- ggplot(plants500) +
  geom_point(aes(y = age_median, x = median_A_b, color = site, size = total_n_reads)) +
  scale_color_manual(values = core_cols) +
  geom_ribbon(data = qdata, aes(y = age_median, xmin = lwr, xmax = upr), alpha = 0.2) +
  geom_path(data = qdata, aes(y = age_median, x = QFILT)) +
  geom_path(data = qdata, color = "blue", aes(y = age_median, x = fit_q)) +
  annotate("rect", xmin = -Inf, xmax = 0.1, ymin = -Inf, ymax = 2000,
           fill = "red", alpha = 0.2) +
  labs(y = "time",
       x = "A_b",
       title = "Damage model")
ggsave(paste0(PLOT_DIRECTORY, "/damage_model.png"), plot = plot)


agg_stats <- agg_stats %>% 
  mutate(status = ifelse(age_median > 2000, "pass", 
                         ifelse(age_median < 2000 & median_A_b < 0.1, "fail", "pass")))

agg_stats <- agg_stats %>%
  group_by(site, genus) %>%
  mutate(
    oldest_pass_date = min(age_median[status == "pass" & total_n_reads >= 50], na.rm = TRUE),
    lenient_status = if_else(age_median >= oldest_pass_date, "pass", status)
  ) %>%
  ungroup()


# -------------  TOP GENERA  ------------- #

df <- get_top(agg_stats %>% filter(PlantAnimal=="plant"), 20) %>% filter(total_n_reads >= 50)
plot <- ggplot(agg_stats %>% filter(genus=="Oryza"), aes(x= median_A_b, y = age_median, fill = status, color = lenient_status, size=total_n_reads)) +
  geom_point(alpha=0.8,shape=21) +
  facet_wrap(~site) 
ggsave(paste0(PLOT_DIRECTORY, "/rice_effectofdmgmodel.png"), plot = plot)



pass_plant <- get_top(agg_stats %>% filter(PlantAnimal=="plant" & lenient_status=="pass"), 20) %>% filter(total_n_reads >= 20)
plot <- ggplot(pass_plant, aes(x = median_A_b, y = age_median, color = site, size = total_n_reads)) +
  facet_wrap(~genus) +
  scale_color_manual(values = core_cols) +
  geom_hline(
    data = events,
    aes(yintercept = age_median),
    linetype = "dotted",
    color = "black"
  ) +
  geom_text(
    data = events,
    aes(y = age_median, label = line_label, x = 1),
    angle = 0,
    hjust = -0.1,
    vjust = 0.5,
    color = "black",
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = centuries) +
  geom_point(alpha = 0.8) 
ggsave(paste0(PLOT_DIRECTORY, "/20_plants_pass.png"))


pass_animal <- get_top(agg_stats %>% filter(PlantAnimal=="animal" & lenient_status=="pass"), 20) %>% filter(total_n_reads >= 20)
plot <- ggplot(pass_animal, aes(x = median_A_b, y = age_median, color = site, size = total_n_reads)) +
  facet_wrap(~genus) +
  scale_color_manual(values = core_cols) +
  geom_hline(
    data = events,
    aes(yintercept = age_median),
    linetype = "dotted",
    color = "black"
  ) +
  geom_text(
    data = events,
    aes(y = age_median, label = line_label, x = 1),
    angle = 0,
    hjust = -0.1,
    vjust = 0.5,
    color = "black",
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_y_continuous(breaks = centuries) +
  geom_point(alpha = 0.8) 
ggsave(paste0(PLOT_DIRECTORY, "/20_animals_pass.png"))


# -------------  ALL GENERA  ------------- #

plot_filtered <- function(df, plotsave, plotname){
  # get annotations about the number of reads in blanks and number of reads in genus
annotations <- df %>%
  group_by(genus, site) %>%
  summarise(
    n_reads = sum(total_n_reads),  # Total reads per site and genus
    .groups = "drop"
  ) %>%
  group_by(genus) %>%
  summarise(
    annotation_text = paste0(site, ": ", n_reads, collapse = "\n"),  # Combine site and reads for each genus
    .groups = "drop"
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
    geom_point(aes(y = age_median, x = median_A_b, color = site, size = total_n_reads)) +
    geom_text(
  data = annotations,
  aes(
    label = annotation_text,  # Use the pre-computed annotation text
    x = Inf, y = Inf          # Position at the top-right of the plot/panel
  ),
  hjust = 1.1, vjust = 1.1,    # Adjust alignment
  inherit.aes = FALSE,         # Do not inherit global plot aesthetics
  size = 3                     # Text size
)+
    labs(
      x = "Median Damage (A_b)",
      y = "Date") +
      scale_color_manual(values = core_cols) +
  geom_hline(
    data = events,
    aes(yintercept = age_median),
    linetype = "dotted",
    color = "black"
  ) +
  geom_text(
    data = events,
    aes(y = age_median, label = line_label, x = 1),
    angle = 0,
    hjust = -0.1,
    vjust = 0.5,
    color = "black",
    size = 3,
    inherit.aes = FALSE
  ) +
    scale_y_continuous(breaks = centuries) +

    facet_params[[1]]

  print(plot)

  }
  dev.off()
}

plot_filtered(pass_animal, paste0(PLOT_DIRECTORY, "/damage.filt.animals.pdf"), "Damage in filtered animals")
plot_filtered(pass_plant, paste0(PLOT_DIRECTORY, "/damage.filt.plants.pdf"), "Damage in filtered plants")


  do_strat_percentage <- function(dat){
    aggregated_data <- dat %>%
        group_by(age_median, genus) %>%
        summarize(n_reads = sum(total_n_reads, na.rm = TRUE), .groups = 'drop') %>%
        group_by(age_median) %>%
        mutate(total_reads = sum(n_reads),
               pct_reads = (n_reads / total_reads) * 100) %>%
        ungroup() %>%
        select(age_median, genus, pct_reads)

    data_wide <- aggregated_data %>%
        pivot_wider(names_from = genus, values_from = pct_reads, values_fill = 0) %>%
        arrange(age_median)

    data_matrix <- as.matrix(data_wide %>% select(-age_median))
    rownames(data_matrix) <- data_wide$age_median

    time_values <- as.numeric(rownames(data_matrix))

    strat.plot(data_matrix, y.rev=FALSE, plot.line=TRUE, plot.poly=FALSE, plot.bar=TRUE,
               lwd.bar=4, sep.bar=TRUE, scale.percent=TRUE, xSpace=0.01,
               x.pc.lab=TRUE, x.pc.omit0=TRUE, srt.xlabel=45, las=2,
               exag=TRUE, exag.mult=5, ylabel = "age", yvar = time_values)
  }




  transform_tibble <- function(tibble, calculation = c("sum", "proportion")) {
    calculation <- match.arg(calculation)  

    tibble <- tibble %>%
      select(age_median, genus, total_n_reads) %>%  
      group_by(age_median, genus) %>%        
      summarize(n_reads = sum(total_n_reads, na.rm = TRUE), .groups = "drop")  

    if (calculation == "proportion") {
      tibble <- tibble %>%
        group_by(age_median) %>%
        mutate(total_reads = sum(n_reads)) %>%  
        mutate(proportion = n_reads / total_reads) %>%  
        select(-total_reads, -n_reads) %>%  
        pivot_wider(names_from = genus, values_from = proportion)
    } else {
      tibble <- tibble %>%
        pivot_wider(names_from = genus, values_from = n_reads)
    }

    tibble <- tibble %>%
      arrange(desc(age_median)) %>%             
      select(age_median, tidyselect::peek_vars()) %>%  
      select(age_median, order(names(.)[-1]) + 1)  

    return(tibble)
  }


for (current_site in unique(pass_plant$site)) {
  # Generate PDF reports
  pdf(file = paste0(PLOT_DIRECTORY, "/", current_site, "_strat_plants_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_plant %>% filter(site == current_site), 20))
  dev.off()

  pdf(file = paste0(PLOT_DIRECTORY, "/", current_site, "_strat_plants_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_plant %>% filter(site == current_site), 50))
  dev.off()

  pdf(file = paste0(PLOT_DIRECTORY, "/", current_site, "_strat_animals_top20_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_animal %>% filter(site == current_site), 20))
  dev.off()

  pdf(file = paste0(PLOT_DIRECTORY, "/", current_site, "_strat_animals_top50_percentage.pdf"), width = 60, height = 15)
  do_strat_percentage(get_top(pass_animal %>% filter(site == current_site), 50))
  dev.off()

  # Transform data
  out_a <- transform_tibble(pass_animal %>% filter(site == current_site), "sum")
  out_p <- transform_tibble(pass_plant %>% filter(site == current_site), "sum")
  p_out_a <- transform_tibble(pass_animal %>% filter(site == current_site), "proportion")
  p_out_p <- transform_tibble(pass_plant %>% filter(site == current_site), "proportion")

  # Create a new workbook
  wb <- createWorkbook()

  # Add worksheets
  addWorksheet(wb, "Animals")
  writeData(wb, "Animals", out_a)

  addWorksheet(wb, "Plants")
  writeData(wb, "Plants", out_p)

  addWorksheet(wb, "ProportionAnimals")
  writeData(wb, "ProportionAnimals", p_out_a)

  addWorksheet(wb, "ProportionPlants")
  writeData(wb, "ProportionPlants", p_out_p)

  # Save the workbook to a file
  EXCEL_OUT <- paste0(PLOT_DIRECTORY, "/", current_site, "_report.xlsx")
  saveWorkbook(wb, EXCEL_OUT, overwrite = TRUE)
}
