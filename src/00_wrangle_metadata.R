# dirs
wd <- "/lustre/scratch126/casm/team154pc/at31/kidney/"
data_dir <- "/nfs/cancer_ref01/nst_links/live/2731"
setwd(wd)
dir.create("out/merged_normal_bams/", showWarnings = F)
dir.create("out/nanoseq/", showWarnings = F)

# libraries
library(magrittr)

# create clean-up function
cleanup_metadata <- function(df) {
    if ("donor_sex" %in% colnames(df)) {
        df <-
            df %>%
            dplyr::mutate(
                donor_sex = dplyr::case_when(
                    donor_sex %in% c("male", "Male") ~ "M",
                    donor_sex %in% c("female", "Female") ~ "F",
                    donor_sex == "Unknown" ~ NA_character_,
                    TRUE ~ donor_sex)) 
    }
    df %>%
        # re-check col types
        readr::type_convert() %>%
        dplyr::distinct()
}

# columns from decode       col_n   (rename?)
# existing_pd_id_for_donor  23       -> sample_id
# gender                    24       -> donor_sex
# donor_age_at_diagnosis    25       -> donor_age
# tissue_phenotype          26       -> sample_phenotype
# tumour or normal?         27       -> sample_tumour_status
# tissue_histology          28       -> sample_histology
# PD_ID                     30       -> biopsy_id
# creating three levels of the metadata:
# 1) donor-level, 2) sample-level, 3) biopsy-level
# all metadata columns will be prefixed with their level

# read decode files
decodes <-
    list.files("data/decode/cp19/", pattern = "xlsx", full.names = TRUE) %>%
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
                donor_sex = gender,
                donor_age = donor_age_at_diagnosis,
                sample_phenotype = tissue_phenotype,
                sample_tumour_status = `tumour or normal?`,
                sample_histology = tissue_histology,
                biopsy_id = PD_ID,
                biopsy_bam_file = paste0(data_dir, "/", biopsy_id, "/", biopsy_id, ".sample.merged.bam")
            )
    })

# convert to sample sheet
sample_sheet <-
    decodes %>%
    # collapse duplicate entries across decode files
    dplyr::bind_rows(.id = "decode_file") %>%
    dplyr::group_by(dplyr::across(c(-decode_file))) %>%
    dplyr::summarise(decode_file = paste(decode_file, collapse = ",")) %>%
    dplyr::ungroup() %>%
    # filter out biopsies with missing bam files
    # cp19@sanger.ac.uk talking about missing BAMs:
    # "The goal of LCMing is to get the minimum amount of tissue needed to reach the DNA
    # input threshold.  My cutting philosophy is to be right at that threshold, which
    # means some samples will end up failing library prep. These proformas are all the
    # samples I submitted for sequencing, not all of the samples that passed sequencing."
    dplyr::filter(
        file.exists(biopsy_bam_file),
        sample_phenotype != "blank") %>%
    # remove donors that have no normal biopsies
    dplyr::group_by(donor_id) %>%
    dplyr::filter("Normal" %in% sample_tumour_status) %>%
    dplyr::ungroup() %>%
    # standardise sample phenotype annotations
    dplyr::mutate(
        sample_phenotype = dplyr::case_when(
            sample_phenotype %in% c("glomerulus", "glomeruls", "Glomeruli") ~ "glomerulus",
            sample_phenotype %in% c("lymphoid", "lymphoidAggregate", "lymphTumor") ~ "lymphoid_aggregate",
            sample_phenotype %in% c("arterial", "venous", "vessel") ~ "endothelial",
            sample_phenotype == "tumor" ~ "tumour",
            sample_phenotype == "hybridStromaDCT" ~ "hybrid_stroma_DCT",
            TRUE ~ sample_phenotype)) %>%
    # standardise donor metadata
    cleanup_metadata()

# check that biopsy_id and sample_id and donor_id match up
check <-
    sample_sheet %>%
    dplyr::distinct(sample_id, donor_id, biopsy_id) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
        donor_id_in_sample_id = stringr::str_sub(sample_id, end = -2) == donor_id,
        donor_id_in_biopsy_id = stringr::str_sub(biopsy_id, end = -9) == donor_id,
        sample_id_in_biopsy_id = stringr::str_sub(biopsy_id, end = -8) == sample_id
    ) %>%
    dplyr::filter(
        donor_id_in_sample_id != T | is.na(donor_id_in_sample_id) |
            donor_id_in_biopsy_id != T | is.na(donor_id_in_biopsy_id) |
            sample_id_in_biopsy_id != T | is.na(sample_id_in_biopsy_id)
    )
if (nrow(check)) {
    message(paste("Check failed!", nrow(check), "IDs are inconsistent!"))
}

# write sample sheet for bams that exist
sample_sheet %>%
    readr::write_csv("out/merged_normal_bams/sample_sheet.csv")

# create sample sheet for nanoseq
sample_sheet_nanoseq <-
    sample_sheet %>%
    dplyr::transmute(
        id = biopsy_id,
        d_bam = biopsy_bam_file,
        n_bam = file.path(wd, "out/merged_normal_bams", paste0(donor_id, "_merged.bam"))
    ) %>%
    # check files exist
    dplyr::filter(file.exists(d_bam) & file.exists(n_bam))
sample_sheet_nanoseq %>%
    readr::write_csv("out/nanoseq/sample_sheet.csv")

# create sample sheet for nanoseq testing
sample_sheet_nanoseq %>%
    head(2) %>%
    readr::write_csv("test/out/sample_sheet.csv")

# read metadata from LCM log
lcm_log <-
    readxl::read_xlsx("data/LCM_log_byStructure_kidney.xlsx") %>%
    dplyr::filter(Sample != "Total") %>%
    dplyr::transmute(
        sample_id = Sample,
        sample_tumour_status = `Tumor?` > 0,
        donor_id = sample_id %>% stringr::str_sub(end = -2),
        donor_predisposition = dplyr::case_when(
            `Predisposition?` == "NA" ~ NA_character_,
            TRUE ~ `Predisposition?`),
        donor_diagnosis = Diagnosis) %>%
    # clean up metadata
    cleanup_metadata()

# read metadata from patient manifest
patient_manifest <- 
    readr::read_tsv("data/patientManifest.txt") %>%
    dplyr::transmute(
        donor_id = Sanger,
        donor_predisposition = Syndrome,
        donor_diagnosis = Dx,
        donor_age = Age,
        donor_sex = Sex
    ) %>%
    # clean up metadata
    cleanup_metadata() 

# read sample info from tjm
samples_1_tjm <-
    "data/decode/tjm/Samples_Tom_24-02-2021_sent.xlsx" %>%
    readxl::read_xlsx() %>%
    dplyr::transmute(
        tube_barcode_identifier = `G number`,
        sample_tumour_status = Region,
        donor_predisposition = Syndrome,
        donor_age = Age,
        donor_sex = `...8`,
        sample_id_tjm = ID
    )

# read decode supplement from tjm
samples_2_tjm <-
    "data/decode/tjm/TBC_TomM_Royal_Free_Renal_Sanger_Human_Proforma_Rec20210309_WithCGP_IDS_20210315.xlsx" %>%
    readxl::read_xlsx(sheet = "Sheet2") %>%
    dplyr::transmute(
        sample_type = `Derivative Type/SampleType`,
        tube_barcode_identifier = `ISBT aliquot 15Digits`,
        donor_diagnosis = `Tumour type`,
        donor_predisposition = Syndrome,
        sample_id_tjm = ID
    )

# read decode files from tjm
sample_sheet_tjm <-
    list.files("data/decode/tjm", pattern = "TBC", full.names = TRUE) %>%
    purrr::set_names(., basename(.)) %>%
    purrr::map(function(file) {
        readxl::read_xlsx(file, sheet = "Sheet1") %>%
            dplyr::filter(dplyr::row_number() > 15) %>%
            dplyr::transmute(
                sample_id = `Donor ID`,
                donor_id = stringr::str_sub(sample_id, end = -2),
                donor_sex = gender,
                donor_age = donor_age_at_diagnosis,
                sample_phenotype = phenotype,
                sample_tumour_status = `tumour?`,
                tube_barcode_identifier,
                sample_id_tjm = Supplier_Donor_ID
            )
    }) %>%
    dplyr::bind_rows(.id = "decode_file")

# harmonise tjm metadata
samples_1_tjm %>%
    dplyr::inner_join(samples_2_tjm, by = "sample_id_tjm") 

# create list of all sources of donor metadata
collapse_donor_metadata <- function(df, suffix) {
    df %>%
        dplyr::select(dplyr::starts_with("donor")) %>%
        dplyr::group_by(donor_id) %>%
        dplyr::summarise(
            dplyr::across(everything(), 
            ~ stringr::str_c(unique(.x), collapse = ","))) %>%
        dplyr::distinct() %>%
        dplyr::rename_with(~ paste0(., ".", suffix), -donor_id)
}
dm <-
    list(
        lcm_log = lcm_log,
        patient_manifest = patient_manifest,
        sample_sheet = sample_sheet) %>%
    purrr::map2(., names(.), ~ collapse_donor_metadata(.x, .y))

# join
dm$sample_sheet %>%
    dplyr::left_join(dm$patient_manifest) %>%
    dplyr::select(order(colnames(.)))

# create list of all sources of sample metadata
collapse_sample_metadata <- function(df, suffix) {
    df %>%
        dplyr::select(dplyr::starts_with("sample"), dplyr::starts_with("donor")) %>%
        dplyr::group_by(donor_id, sample_id) %>%
        dplyr::summarise(
            dplyr::across(everything(), 
            ~ stringr::str_c(unique(.x), collapse = ","))) %>%
        dplyr::distinct() %>%
        dplyr::rename_with(~ paste0(., ".", suffix), -c(sample_id, donor_id))
}
sm <-
    list(
        lcm_log = lcm_log,
        sample_sheet = sample_sheet %>% dplyr::mutate(donor_decode_file = decode_file),
        sample_sheet_tjm = sample_sheet_tjm) %>%
    purrr::map2(., names(.), ~ collapse_sample_metadata(.x, .y))
