
# ----------------- THE BIG FILTER & STATS MACHINE ----------------- #

# ANI must be 95+

# Adds a column "rm"
# Keep if penalized_weighted_median_breadth_exp_ratio > 0.8 and either penalized_weighted_median_gini < 0.6 or penalized_weighted_median_entropy > 0.75
# Remove taxa without a genus

# Join bamfilter data and metadmg data together with taxonomy info
# Get this at genus level

get_stats <- function(tax_ids, tax_data, fb_data, dat_filt, metadata = NULL, mode="library"){
  agg_stats <- tax_ids |>
  inner_join(tax_data) |>
  inner_join(fb_data |> rename(accession.version = reference)) |>
  mutate(flt = paste(label, species, sep = "--")) |>
  inner_join(dat_filt |>
               select(label, name, A_b, c_b, fit, PlantAnimal) |>
               mutate(flt = paste(label, name, sep = "--")) |>
               select(-label, -name)) |>
  filter(read_ani_median >= 95, n_reads >= 1) |>
  as_tibble() |>
  select(-flt) |>
    group_by(label, genus, superkingdom, phylum, class, order, family, fit, PlantAnimal) |>
  mutate(
    num_alns = sum(n_alns),
    weight = n_alns / num_alns,
    scaled_breadth_exp_ratio = breadth_exp_ratio * weight
  ) |>
  ungroup() |>
    group_by(label, genus, superkingdom, phylum, class, order, family, fit, PlantAnimal) |>
  summarise(
    n = n(),
    median_A_b = median(A_b),
    penalized_weighted_median_A_b =penalized_weighted_median(A_b, n_reads, A_b),
    median_c_b = median(c_b),
    mean_read_ani_median = mean(read_ani_median),
    mean_read_ani_std = mean(read_ani_std),
    median_reference_length = median(reference_length),
    sum_reference_length = sum(reference_length),
    mean_breadth_exp_ratio = mean(breadth_exp_ratio),
    median_breadth_exp_ratio = median(breadth_exp_ratio),
    penalized_weighted_median_entropy = penalized_weighted_median(norm_entropy, n_reads, norm_entropy),
    penalized_weighted_median_gini = penalized_weighted_median(norm_gini, n_reads, norm_gini),
    penalized_weighted_median_breadth_exp_ratio = penalized_weighted_median(breadth_exp_ratio, n_reads, breadth_exp_ratio),
    median_entropy = median(norm_entropy),
    mean_entropy = mean(norm_entropy),
    median_gini = median(norm_gini),
    mean_gini = mean(norm_gini),
    median_n_reads = median(n_reads),
    mean_n_reads = mean(n_reads),
    n_reads = sum(n_reads)
  ) |>
  ungroup() |>
  mutate(rm = ifelse(penalized_weighted_median_breadth_exp_ratio > 0.8 & (penalized_weighted_median_gini < 0.6 | penalized_weighted_median_entropy > 0.75), "keep", "remove")) |>
  filter(!is.na(genus) & !is.na(PlantAnimal)) |> filter(rm == "keep")
  if (mode == "replicates") {
    agg_stats$cgg <- sapply(agg_stats$label, function(x) strsplit(x,"_")[[1]][1])
    agg_stats <- inner_join(agg_stats, metadata, by="cgg")
  }
  if (mode == "library"){ 
    agg_stats <- inner_join(agg_stats, metadata) 
  }

  return(agg_stats)
}