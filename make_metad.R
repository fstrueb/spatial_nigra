library(tidyverse)
pwd <- '/home/fstruebi/projects/spatial_nigra_v3/'
make_metad <- function() {
    excludesamples <- read_lines(paste0(pwd, 'resources/excludesamples.txt'))
    metad <- read_tsv(paste0(pwd, 'resources/essential_metad_combined_cohorts.txt'), col_types = cols()) %>% 
        mutate(case = case_when(
            group %in% c('AD', 'ADLBD') ~ paste0('NBB', str_replace(case, '-', '_')),
            TRUE ~ case
        )) %>% 
        mutate(braak_lewy_inferred = case_when(
            is.na(braak_lewy) ~ '6', # one PD case, assume it's Braak 6
            braak_lewy == 'Amygdala focal' ~ '3', # one iLBD case, assume it's B3 if Amygdala shows LBP
            TRUE ~ braak_lewy
        )) %>% 
        mutate(braak_lewy_inferred = as.numeric(braak_lewy_inferred)) %>% 
        mutate(braak_tau = ifelse(is.na(braak_tau), 1, braak_tau))
    saa_auc <- read_tsv(paste0(pwd, 'resources/SAA_AUC_table.txt'), col_types = cols())
    saa_binned <- saa_auc %>% 
        mutate(saa_bin = case_when(median_auc < 5e5 ~ 'low',
                                   median_auc > 5e5 & median_auc <= 1e6 ~ 'medium',
                                   median_auc > 1e6 ~ 'high')) %>% 
        mutate(saa_bin_factor = case_when(
            saa_bin == 'low' ~ 1,
            saa_bin == 'medium' ~2,
            saa_bin == 'high' ~3
        ))
    metadqc <- read_tsv(paste0(pwd, 'resources/metad_qcmeasures.tsv'), col_types = cols()) %>% select(1, 3:26)
    metadbatch <- read_tsv(paste0(pwd, 'resources/spatial_metad_with_batch_info.txt'), col_types = cols()) %>% 
        mutate(kit_batch = ifelse(batch == 2 & grepl('T', slide_id), 1, 2)) %>% 
        mutate(kit_batch = as.factor(ifelse(batch == 1, 1, kit_batch))) %>% 
        select(sample, kit_batch)
    metad <- metad %>% 
        left_join(saa_binned, by = c('case' = 'sample')) %>% 
        left_join(metadbatch, by = c('case' = 'sample')) %>% 
        left_join(metadqc, by = c('case' = 'sample')) %>% 
        filter(!case %in% excludesamples)
}
