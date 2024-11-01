
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
source("butterfly/dmg.R")
source("butterfly/get_calculate_plot_grid.R")
source("butterfly/perk.R")
source("butterfly/damage_est_function.R")
source("butterfly/perk_wrapper.R")
source("butterfly/perk_wrapper_function.R")
source("butterfly/get_dmg_decay_fit.R")
source("butterfly/median.R")
source("butterfly/filter.R")

# ----------------- ARGS  ----------------- #

# probably don't need to change 
NAMES <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/names.dmp"
NODES <- "/projects/caeg/data/db/aeDNA-refs/resources/20230825/ncbi/taxonomy/nodes.dmp"
ACC2TAXID <- "/projects/lundbeck/scratch/for_antonio/ncbi_taxonomy_01Oct2022/combined_accession2taxid_20221112.gz"

# change per project
METADATA <- "new_metadata.tsv"
PLOT_DIRECTORY <- "finalplotsV2/"
METADMG_DIR <- "raw_results/stats/metadmg/aggregate/"
BAMFILTER_DIR <- "raw_results/results/bamfilter/"

# negative data 
NEG_METADMG_DIR <- "raw_results/stats/metadmg/aggregate/lib/"
NEG_BAMFILTER_DIR <- "raw_results/results/bamfilter/lib"

# have we already run this? if so we can save a lot of time by using save taxonomy info
TAXONOMY_TSV <- "saved_taxids.tsv" 

# save agg stats here so we can skip this if needed
AGG_TSV <- "agg_stats.tsv"
NEGAGG_TSV <- "neg-agg_stats.tsv"

# save damage model here 
GLOBAL_DMG_MODEL <- "qdata.tsv"

# final excel 
EXCEL_OUT <- "TJORNIN.XLSX"

# ----------------- READ IN DATA  ----------------- #

# ------------- METADATA  ------------- #

metadata <- read.table(METADATA,he=T,as.is=T)
metadata$label <- gsub("_", "", metadata$CGG_ID)
names(metadata)[2] <- "date" # we know this column should be called date, this must be changed per project 
metadata <- metadata[, c("label","date")]



# ------------- HELPER TO READ LOTS OF FILES  ------------- #

read_file <- function(f) {
	df <- fread(f, header=T, sep="\t", fill=T, nThread=20)
	df$sample_id <- f
	return(df)
}



# ------------- METADMG  ------------- #

file_list <- list.files(path = METADMG_DIR, pattern = "*gz", recursive = FALSE, full.name=TRUE)
holi_data <- do.call(rbind, lapply(file_list, read_file))
holi_data$label <- sub(".*(CGG[0-9]+)\\..*", "\\1", holi_data$sample_id) # this will have to change when we have non-CGG ids 


# merge metadmg data with metadata 

holi_data <- inner_join(holi_data, metadata, by = "label")


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


dat_filt <- filter_metadmg(holi_data, metadata$label |> unique())



# ------------- BAMFILTER  ------------- #

read.names.sql(NAMES, sqlFile = "nameNode.sqlite", overwrite=TRUE)
read.nodes.sql(NODES, sqlFile = "nameNode.sqlite", overwrite=TRUE)

file_list <- list.files(path = BAMFILTER_DIR, pattern = "*stats.tsv.gz", recursive = FALSE, full.name=TRUE)
fb_data <- do.call(rbind, lapply(file_list, read_file))
fb_data$label <- sub(".*(CGG[0-9]+)\\..*", "\\1", fb_data$sample_id)

# take the long or fast route depending on if we have the taxid info saved already 
if (file.exists(TAXONOMY_TSV)){
	tax_ids <- read_table(TAXONOMY_TSV)
} else { 
	# unique accession IDs 
	accs <- fb_data |>
		select(accession.version = reference) |>
		distinct() |>
		pull(accession.version)

	# read acc2taxid file 
	acc2taxid <- fread(ACC2TAXID, 
		tmpdir = "temp/",
		nThread = 48,
		showProgress = TRUE)

	# filter just those used
	tax_ids <- acc2taxid %>%
		filter(accession.version %in% accs)

	# save for next time 
	write_tsv(tax_ids, TAXONOMY_TSV)
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



# ----------------- MERGE AND AGGREGATE TO GENUS LEVEL  ----------------- #


agg_stats <- get_stats(tax_ids, tax_data, fb_data, dat_filt, metadata=metadata, mode="library")

if (!file.exists(AGG_TSV)){
	write.csv(agg_stats, AGG_TSV)
}



# ----------------- READ NEGATIVE DATA  ----------------- #

# read negative metadmg data 
file_list <- list.files(path = NEG_METADMG_DIR, pattern = "*gz", recursive = FALSE, full.name=TRUE)
neg_holi_data <- do.call(rbind, lapply(file_list, read_file))
neg_holi_data$label <- sapply(neg_holi_data$sample_id, function(x) strsplit(x,"_")[[1]][3])
neg_dat_filt <- filter_metadmg(neg_holi_data, neg_holi_data$label |> unique())

# read negative bam
file_list <- list.files(path = NEG_BAMFILTER_DIR, pattern = "*stats.tsv.gz", recursive = FALSE, full.name=TRUE)
neg_fb_data <- do.call(rbind, lapply(file_list, read_file))
neg_fb_data$label <- sapply(neg_fb_data$sample_id, function(x) strsplit(x,"_")[[1]][3])

# aggregate to genus level 
neg_agg_stats <- get_stats(tax_ids, tax_data, neg_fb_data, neg_dat_filt, mode="negatives") %>%
  mutate(genus = factor(genus, levels = sort(unique(genus),decreasing=TRUE))) %>% 
  group_by(label) %>%
  ungroup()

# write negative tsv 
if (!file.exists(NEGAGG_TSV)){
	write.csv(neg_agg_stats, NEGAGG_TSV)
}

# plot a nice heatmap of the number of reads in each genus in the blanks 
ggplot(neg_agg_stats, aes(x = label, y = genus, fill = n_reads)) +
  geom_tile() +
  geom_text(aes(label = n_reads), color = "black", size = 3) +
  scale_fill_gradientn(colors = wes_palette("Zissou1"),values=scales::rescale(c(min(neg_agg_stats$n_reads),median(neg_agg_stats$n_reads),max(neg_agg_stats$n_reads)))) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ggtitle("Number of reads for each genus in negative controls")
ggsave(paste0(PLOT_DIRECTORY, "/negatives_heatmap.png"))

# save all the genera in the blanks for later 
contam_genera <- neg_agg_stats %>%
	group_by(genus) %>%
	summarise(
	  reads_blanks = sum(n_reads) 
	) %>% arrange(desc(reads_blanks))



# ----------------- BASIC PLOTS  ----------------- #

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
  group_by(date, rank) %>%
  summarize(numreads = sum(nreads)) %>% # number of reads per rank 
  mutate(proportion = numreads / sum(numreads)) # proportion 

# plot ranks classified (raw numbers)
plot_total <- ggplot(df, aes(x = as.factor(date), y = numreads, fill = rank)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + labs(x = NULL)

# plot proportion
plot_proportion <- ggplot(df, aes(x = as.factor(date), y = proportion, fill = rank)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
  scale_y_continuous(labels = scales::percent)

# combine 
combined_plot <- plot_total / plot_proportion +
    plot_layout(guides = "collect") & theme(legend.position = "bottom")

ggsave(paste0(PLOT_DIRECTORY, "/A.euk.ranks.png"))



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


plot_example_fits(dat_filt, "good", paste0(PLOT_DIRECTORY, "/B.good_examples.100reads.pdf"), 100, 100)
plot_example_fits(dat_filt, "bad", paste0(PLOT_DIRECTORY, "/B.bad_examples.100reads.pdf"), 100, 100)



# ----------- C) damage all data ----------- #

ggplot(agg_stats[agg_stats$n_reads>100,], aes(x = median_A_b, y = date, color = fit, size=n_reads)) +
  geom_point(alpha=0.6) + facet_wrap(~PlantAnimal)

ggsave(paste0(PLOT_DIRECTORY, "/C.dmg.alldata.100reads.png"))


# ----------- D) damage in "good" plants 500  ----------- #

# small function to get the n most abundant genera in either plant/animal 
get_top <- function(df, n, planimal){
  tmp <- df %>%
  filter(PlantAnimal == planimal) %>%
  group_by(genus) %>%
  summarize(total_reads = sum(n_reads, na.rm = TRUE)) %>%
  arrange(desc(total_reads)) %>%
  slice_head(n = n)
  return(df %>% filter(genus %in% tmp$genus))
}

# plot 20 most abundant plants, with good and bad fit, only those with >500 reads 
df <- get_top(agg_stats, 20, "plant") %>% filter(n_reads >= 500)
ggplot(df, aes(x = median_A_b, y = date, color = fit, size=n_reads)) +
  geom_point(alpha=0.8) +
  facet_wrap(~genus)
ggsave(paste0(PLOT_DIRECTORY, "/D.dmg.plants.500reads.png"))


# ----------- E) damage model  ----------- #


get_conditional_quantile <- function(df) {
	# get 10 first and last dates, these should be filtered differently 
	latest_dates <- df %>%
		distinct(date) %>%
		arrange(date) %>%
		slice(c(1:10, (n()-9):n()))


	# find damage quantiles 
	qdata <- df %>%
		mutate(is_first_dates = date %in% latest_dates$date) %>%
		group_by(date, is_first_dates) %>%
		filter(n() >= 10) %>%
	    summarise(
	      QFILT = if (first(is_first_dates)) {
	        quantile(median_A_b, probs = 0.05, na.rm = TRUE) # keep all taxa 
	      } else {
	        quantile(median_A_b[fit == 'good'], probs = 0.05, na.rm = TRUE) # just good fit for the middle dates 
	      }, .groups = 'drop'
	    ) 

  return(qdata)
}

# get plants with > 500 reads 
plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & n_reads >= 500)

# get damage line 
qdata <- get_conditional_quantile(plants500)
loess_fit <- loess(QFILT ~ date, data = qdata, span = 0.3) # loess fit to smooth it 
pred <- predict(loess_fit, newdata = qdata, se = TRUE)
qdata <- qdata %>%
  mutate(
    fit_q = pred$fit,
    lwr = pred$fit - 1.96 * pred$se.fit,
    upr = pred$fit + 1.96 * pred$se.fit,
  ) # get 5% percentile, use this as the limit 
write.table(qdata, GLOBAL_DMG_MODEL)

# assign pass or fail if passed line, if not one of the first or last dates it must also be good fit 
agg_stats <- inner_join(agg_stats, qdata) %>% mutate(status = ifelse(median_A_b < lwr, "fail", "pass"),
  status = ifelse(!is_first_dates & fit == "bad", "fail", status))

# for plotting, assign a low  alpha val to the failed taxa 
plants500 <- agg_stats %>% filter(PlantAnimal == "plant" & n_reads >= 500) %>% 
	mutate(alpha = ifelse(status=="pass",1,0.2))

# plot damage model and the points that made it 
ggplot(plants500) +
  geom_point(aes(y = date, x = median_A_b, color=fit, alpha=alpha)) +
  geom_ribbon(data = qdata, aes(y = date, xmin = lwr, xmax = upr), alpha = 0.2) +
  geom_path(data = qdata, aes(y = date, x = QFILT)) +
  geom_path(data = qdata, color="blue", aes(y = date, x = fit_q)) +
  scale_alpha_identity() +
    labs(y = "Date",
       x = "A_b",
       title = "Damage model")
ggsave(paste0(PLOT_DIRECTORY, "/E.dmgmodel.png"))



# ----------------- FILTERING USING THE DAMAGE MODEL ----------------- #


# find passing genera with at least reads at occ different dates 
filter_occ <- function(df, reads, occ){
  filtered_data <- df %>%
  group_by(genus, label) %>%
  filter(n_reads >= reads & status == "pass") %>%
  group_by(genus) %>%
  filter(n_distinct(label) >= occ) %>%
  ungroup()
  unique(filtered_data$genus)
}


good_genera <- filter_occ(agg_stats, 100, 2) # at least 2 occurences of 100 reads, pass

agg_stats <- agg_stats %>%
  filter(genus %in% good_genera)


# ----------- F) dmg in filtered data  ----------- #

plot_filtered <- function(qdata, df, plotsave, plotname){
	# get annotations about the number of reads in blanks and number of reads in genus
	annotations <- df %>%
		group_by(genus) %>%
		summarise(
		  n_reads = sum(n_reads),
		  pass_proportion = mean(status == "pass")
		) %>%
		left_join(contam_genera) %>%
		mutate(reads_blanks = scales::comma(reads_blanks), pass_proportion = scales::percent(pass_proportion), n_reads = scales::comma(n_reads))


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
	  geom_ribbon(data = qdata, aes(y = date, xmin = lwr, xmax = upr), alpha = 0.4) + # damage model area 
	  geom_path(data = qdata, aes(y = date, x = fit_q)) + # damage model mean 
	  geom_point(aes(y = date, x = median_A_b, color = status, size = n_reads, alpha = ifelse(status == "pass", 1, 0.4)), show.legend = c(color = TRUE, alpha = FALSE)) +
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
	    y = "Date (CE)"
	  ) +
	   geom_text(
	      data = annotations, # reads in blanks annotation 
	          aes(
	              label = paste(

	                  ifelse(!is.na(reads_blanks), paste("\nReads in Blanks: ", reads_blanks), "")
	              ),
	              x = Inf, y = min(df$date)
	          ),
	          hjust = 1.1, vjust = 0,
	          inherit.aes = FALSE, size = 3, color="red" # red for angry 
	      ) +

	  xlim(0,0.5)+  
	  facet_params[[1]]

	print(plot)

	}
	dev.off()
}

# seperate animals and plants 
plant_df <- agg_stats[agg_stats$PlantAnimal == "plant",]
animal_df <- agg_stats[agg_stats$PlantAnimal == "animal",]

plot_filtered(qdata, animal_df, paste0(PLOT_DIRECTORY, "/F.damage.filt.animals.pdf"), "Damage in filtered animals")
plot_filtered(qdata, plant_df, paste0(PLOT_DIRECTORY, "/F.damage.filt.plants.pdf"), "Damage in filtered plants")

# get just passing genera 
pass_animal_df <- animal_df[animal_df$status == "pass",]
pass_plant_df <- plant_df[plant_df$status == "pass",]


# ----------- G) ani, gini, breadth in filtered data  ----------- #

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
      geom_point(aes(y = date, x = !!sym(var), color = median_A_b, size = n_reads)) +
      labs(
        y = "Date (CE)"
      ) +
      facet_params[[1]]

    print(plot)

    }
    dev.off()
}

plot_filtered_var(qdata, pass_animal_df, "mean_read_ani_median", paste0(PLOT_DIRECTORY, "/G.ani.filt.animals.pdf"), "ANI in filtered animals")
plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_breadth_exp_ratio", paste0(PLOT_DIRECTORY, "/G.breadth.filt.animals.pdf"), "Breadth in filtered animals")
plot_filtered_var(qdata, pass_animal_df, "penalized_weighted_median_gini", paste0(PLOT_DIRECTORY, "/G.gini.filt.animals.pdf"), "Gini in filtered animals")

plot_filtered_var(qdata, pass_plant_df, "mean_read_ani_median", paste0(PLOT_DIRECTORY, "/G.ani.filt.plants.pdf"), "ANI in filtered plants")
plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_breadth_exp_ratio", paste0(PLOT_DIRECTORY, "/G.breadth.filt.plants.pdf"), "Breadth in filtered plants")
plot_filtered_var(qdata, pass_plant_df, "penalized_weighted_median_gini", paste0(PLOT_DIRECTORY, "/G.gini.filt.plants.pdf"), "Gini in filtered plants")


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
pdf(file = paste0(PLOT_DIRECTORY, "/H.strat_plants_top20_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_plant_df, 20, "plant"))
dev.off()

pdf(file = paste0(PLOT_DIRECTORY, "/H.strat_plants_top50_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_plant_df, 50, "plant"))
dev.off()


pdf(file = paste0(PLOT_DIRECTORY, "/H.strat_animals_top20_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_animal_df, 20, "animal"))
dev.off()


pdf(file = paste0(PLOT_DIRECTORY, "/H.strat_animals_top50_percentage.pdf"), width = 60, height = 15)
do_strat_percentage(get_top(pass_animal_df, 50, "animal"))
dev.off()


# ----------------- EXCEL OUT ----------------- #

# sum or get proportion of reads, so genera are columns and dates are rows 
transform_tibble <- function(tibble, calculation = c("sum", "proportion")) {
  calculation <- match.arg(calculation)  

  tibble <- tibble %>%
    select(date, genus, n_reads) %>%  
    group_by(date, genus) %>%        
    summarize(n_reads = sum(n_reads, na.rm = TRUE), .groups = "drop")  

  if (calculation == "proportion") {
    tibble <- tibble %>%
      group_by(date) %>%
      mutate(total_reads = sum(n_reads)) %>%  
      mutate(proportion = n_reads / total_reads) %>%  
      select(-total_reads, -n_reads) %>%  
      pivot_wider(names_from = genus, values_from = proportion)
  } else {
    tibble <- tibble %>%
      pivot_wider(names_from = genus, values_from = n_reads)
  }

  tibble <- tibble %>%
    arrange(desc(date)) %>%             
    select(date, tidyselect::peek_vars()) %>%  
    select(date, order(names(.)[-1]) + 1)  

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

addWorksheet(wb, "GeneraInBlanks")
writeData(wb, "GeneraInBlanks", contam_genera)

# save the workbook to a file
saveWorkbook(wb, paste0(PLOT_DIRECTORY, "/", EXCEL_OUT), overwrite = TRUE)
