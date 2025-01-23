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

get_ass_seqs <- function(file_path) {
  cmd <- paste("zcat", shQuote(file_path), "|sed 1d|head -n1|cut -f4")
  total_seqs <- as.numeric(system(cmd, intern = TRUE))
  if (is.na(total_seqs)) {
    return(0)
  }
  return(total_seqs)
}

get_reads_df <- function(df) {

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
    arrange(-depth) %>%  
    mutate(label = factor(library_id, levels = unique(library_id))) 

  # now sum across lanes 
  summed_df <- lane_df %>%
    group_by(label, library_id) %>%
    summarize(
      across(c(basedir,basename, depth), first),
      total_seqs = sum(total_seqs),
      collapsed_seqs = sum(collapsed_seqs),
      .groups = "drop"
    ) %>% 
    arrange(-depth) %>%  
    mutate(library_id = factor(library_id, levels = unique(library_id)))  


  # Assigned sequences

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
      assigned_seqs = get_ass_seqs(agg_file)
    ) %>%
    ungroup()

  # join with other data
  summed_df <- inner_join(ass_df, summed_df)

  # get % 
  summed_df <- summed_df %>%
    mutate(collapsed_pct = (collapsed_seqs / total_seqs) * 100,
    assigned_pct = (assigned_seqs / collapsed_seqs) * 100)

  summed_df
}