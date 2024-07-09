library(tidyverse)

# set working directory
setwd("/Users/npl206/Google Drive/My Drive/MediterreneanProject/sapropels/")

# metadata
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/get-metadata.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/parse_metadata.R")

# Execute R scripts to enable functions for the Concordance Correlation Coefficient calculations and more from the perk repo
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/damage_est_function.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk_wrapper.R")
source("~/Library/CloudStorage/GoogleDrive-mikkel.w.pedersen@gmail.com/My Drive/MediterreneanProject/sapropels/perk_wrapper_function.R")

# Let's get the damage at local level
dmg <- read_tsv("tp-damage.tsv.gz") %>%
  mutate(label = paste(str_extract(label, "(MED)-(\\d+)-(\\d+|NC\\d+|NL\\d+)"), sep = "-")) %>%
  filter(label %in% cdata$label) # %>%
#### inner_join(tax_data %>% select(label, taxid, name, rank))

# For species
# Let's identify all those taxa with a bad fit
# We set a damage
dmg_sp <- dmg %>%
  filter(rank == "species") %>%
  filter(nreads >= 100) %>%
  rename(tax_name = taxid, n_reads = nreads)

samples <- cdata$label %>% unique()
dat <- dmg_fwd_CCC(dmg_sp, samples, ci = "asymptotic", nperm=100, nproc=10)


# dat %>%
#     ungroup() %>%
#         inner_join(dmg_sp) %>%
#     inner_join(cdata %>% filter(geographic_location == "South of Cyprus")) %>%
#     filter(A_b < 0.1) %>%
#     ggplot(aes(rho_c, A_b, color = label)) +
#     geom_point() 

dat1 <- dat %>%
  ungroup() %>%
  inner_join(dmg_sp) %>%
  mutate(fit = ifelse(rho_c >= 0.85 & C_b > 0.9 & round(rho_c_perm_pval, 3) < 0.1 & !is.na(rho_c), "good", "bad")) %>%
  mutate(fit = ifelse(q_CI_h >= 1 | c_CI_l <= 0, "bad", fit))

dat2 <- dat1 %>%
  filter(fit == 'good' & nalign/n_reads > 1 & grepl("Eukaryota", taxa_path))

dat_bad <- dat1 %>% filter(fit == 'bad')