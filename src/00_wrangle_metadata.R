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
    dplyr::bind_rows() %>%
    dplyr::distinct() %>%
    # remove problematic samples for now
    dplyr::filter(!(sample_id %in% c("PD44966b", "PD43948s"))) %>%
    # check if bam file exists
    dplyr::mutate(
        biopsy_bam_file = ifelse(file.exists(biopsy_bam_file), biopsy_bam_file, NA)
    )

# write sample sheet for bams that exist
sample_sheet %>%
    dplyr::filter(!is.na(biopsy_bam_file)) %>%
    readr::write_csv("data/sample_sheet.csv")

# read patient metadata
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