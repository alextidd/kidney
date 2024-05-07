# dirs
wd <- "/lustre/scratch126/casm/team154pc/at31/kidney/"
data_dir <- "/nfs/cancer_ref01/nst_links/live/2731"
setwd(wd)

# libraries
library(magrittr)

# columns from decode       col_n   (rename?)
# existing_pd_id_for_donor  23       -> sample_id
# gender                    24       -> donor_gender
# donor_age_at_diagnosis    25       -> donor_age
# tissue_phenotype          26       -> sample_phenotype
# tumour or normal?         27       -> sample_tumour_status
# tissue_histology          28       -> sample_histology
# PD_ID                     30       -> biopsy_id
# creating three levels of the metadata:
# 1) donor-level, 2) sample-level, 3) biopsy-level
# all metadata columns will be prefixed with their level

# read decode files
sample_sheet <-
    list.files("data/decode", pattern = "xlsx", full.names = TRUE) %>%
    purrr::set_names(., basename(.)) %>%
    purrr::map(function(file) {
        file %>%
            readxl::read_xlsx(skip = 1) %>%
            dplyr::filter(
                dplyr::row_number() > 12,
                !is.na(existing_pd_id_for_donor)) %>%
            dplyr::transmute(
                sample_id = existing_pd_id_for_donor,
                donor_id = stringr::str_sub(sample_id, end = -2),
                donor_gender = gender,
                donor_age = donor_age_at_diagnosis,
                sample_phenotype = tissue_phenotype,
                sample_tumour_status = `tumour or normal?`,
                sample_histology = tissue_histology,
                biopsy_id = PD_ID,
                biopsy_bam_file = paste0(data_dir, "/", biopsy_id, "/", biopsy_id, ".sample.merged.bam")
            ) 
    }) %>%
    # collapse duplicate entries across decode files
    dplyr::bind_rows(.id = 'decode_file') %>%
    dplyr::group_by(dplyr::across(c(-decode_file))) %>%
    dplyr::summarise(decode_file = paste(decode_file, collapse = ',')) %>%
    dplyr::ungroup() %>%
    # filter out biopsies with missing bam files
    # cp19@sanger.ac.uk talking about missing BAMs:
    # "The goal of LCMing is to get the minimum amount of tissue needed to reach the DNA 
    # input threshold.  My cutting philosophy is to be right at that threshold, which 
    # means some samples will end up failing library prep.  These proformas are all the 
    # samples I submitted for sequencing, not all of the samples that passed sequencing."
    dplyr::filter(file.exists(biopsy_bam_file)) %>%
    # remove donors that have no normal biopsies
    dplyr::group_by(donor_id) %>%
    dplyr::filter('Normal' %in% sample_tumour_status) %>%
    dplyr::ungroup()

# check that biopsy_id and sample_id and donor_id match up
check <-
    sample_sheet %>%
    dplyr::distinct(sample_id, donor_id, biopsy_id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        donor_id_in_sample_id = stringr::str_sub(sample_id, end = -2) == donor_id,
        donor_id_in_biopsy_id = stringr::str_sub(biopsy_id, end = -9) == donor_id,
        sample_id_in_biopsy_id = stringr::str_sub(biopsy_id, end = -8) == sample_id) %>% 
    dplyr::filter(
        donor_id_in_sample_id != T | is.na(donor_id_in_sample_id) |
        donor_id_in_biopsy_id != T | is.na(donor_id_in_biopsy_id) |
        sample_id_in_biopsy_id != T | is.na(sample_id_in_biopsy_id))
if(nrow(check)) {
    message(paste('Check failed!', nrow(check), 'IDs are inconsistent!'))
}

# write sample sheet for bams that exist
sample_sheet %>%
    readr::write_csv("data/sample_sheet.csv")

# create sample sheet for nanoseq
sample_sheet_nanoseq <-
    sample_sheet %>%
    dplyr::transmute(
        id = biopsy_id,
        d_bam = biopsy_bam_file,
        n_bam = file.path(wd, 'out/merged_normal_bams', paste0(donor_id, '_merged.bam'))) %>% 
    # check files exist
    dplyr::filter(file.exists(d_bam) & file.exists(n_bam)) %>%
    readr::write_csv("data/sample_sheet_nanoseq.csv") 

# read patient metadata from Chloe and Tom
patient_metadata1 <-
    readxl::read_xlsx("data/LCM_log_byStructure_kidney.xlsx") %>%
    dplyr::filter(Sample != "Total") %>%
    dplyr::transmute(
        sample_id = Sample, 
        sample_tumour_status = `Tumor?`, 
        donor_id = sample_id %>% stringr::str_sub(end = -1),
        donor_predisposition = `Predisposition?`,
        donor_diagnosis = Diagnosis)
patient_metadata2 <-
    readr::read_csv("data/patientManifest.txt")
patient_metadata3 <-
    list.files("data/tjm", pattern = "TBC", full.names = TRUE) %>%
    purrr::set_names(., basename(.)) %>%
    purrr::map(function(file) {
        readxl::read_xlsx(file, sheet = "Sheet1") %>%
        dplyr::filter(dplyr::row_number() > 15) %>%
        dplyr::transmute(
            sample_id = `Donor ID`,
            donor_id = stringr::str_sub(sample_id, end = -2),
            donor_gender = gender,
            donor_age = donor_age_at_diagnosis,
            sample_phenotype = phenotype,
            sample_tumour_status = `tumour?`
        ) 
    }) %>%
    dplyr::bind_rows(.id = 'decode_file') 

