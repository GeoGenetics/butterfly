library(data.table)
library(dplyr)

# ------------ args  ------------ #
args <- commandArgs(trailingOnly = TRUE)
country <- toupper(args[1])
latest_filelist <- args[2]
out_csv <- args[3]

# ------------ comlicated stuff to get the DB directory location ------------ #
get_script_dir <- function() {
  cmd_args <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmd_args)
  return(normalizePath(dirname(sub(needle, "", cmd_args[match]))))
}

script_dir <- get_script_dir()
butterfly_root <- normalizePath(file.path(script_dir, ".."))
db_path <- file.path(butterfly_root, "db")
db_path <- normalizePath(db_path)


# ------------ function for cleaning up colnames ------------ #
read_clean <- function(f){
	df <- fread(f,header=T) 
	names(df) <- make.names(names(df), unique=TRUE)
	return(df)
}


# ------------ do the stuff ------------ #

# get field samples, and just get the country we want
field_samples <- read_clean(paste0(db_path, "/field_sample.tsv"))
field_samples <- field_samples %>%
  filter(toupper(Country.Ocean) == country) %>%
  select(Unique.Sample.ID) 
names(field_samples)[1] <- "BulkSampleID"

# get archive samples, where the bulkid matches the ones we want 
archive_samples <- read_clean(paste0(db_path, "/edna_archive_sample.tsv"))

samples <- archive_samples %>% inner_join(field_samples) %>% 
	select(BulkSampleID,ArchiveSampleID, DepthSampledCalTape) 


# get library names of interesting archive samples 
wetlab_report <- read_clean(paste0(db_path, "/edna_wetlab_report.tsv"))

names(wetlab_report)[11] <- "ArchiveSampleID"
samples <- samples %>% inner_join(wetlab_report) %>% 
	filter(FastQ.File.ID != "") %>% 
	select(BulkSampleID,ArchiveSampleID, DepthSampledCalTape, Library.ID) %>% 
	mutate(Filebase=paste0(ArchiveSampleID,"_",Library.ID)) 

# load automated filelists 
file_df <- read.csv(latest_filelist, header=F)
names(file_df) <- c("Basedir", "Filebase")
  
samples <- inner_join(samples, file_df) %>% 
	filter(!grepl("_HNLW5DSX5_", Basedir)) # this is needed due to small bug in production pipeline, see email to filipe 23-07

# out 
write.csv(samples, out_csv, row.names=FALSE, quote = FALSE)
print(paste0("Written csv to: ", out_csv))