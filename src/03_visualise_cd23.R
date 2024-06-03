#!/usr/bin/env Rscript

# this script will process the nanoseqV2 pipeline output and generate:
# -> table of corrected mutations burdens - identify samples for exclusion due
#    to contamination, sorting and sequencing problems
# -> table of corrected tinucleotide mutations frequencies
# -> table of all "PASS" indels
# -> table of observed mutation burdens

# load libraries
library(magrittr)

# load sample sheet
sample_sheet <-
  readr::read_csv("out/nanoseq/sample_sheet.csv")

# get output files
post_files <-
  list(".mut_burden.tsv$", ".trint_subs_obs_corrected.tsv$") %>%
  purrr::set_names(., .) %>%
  purrr::map(function(ext) {
    sample_sheet$id %>%
      purrr::set_names(., .) %>%
      purrr::map(function(biopsy_id) {
        list.files(paste0("out/nanoseq/outNextflow/", biopsy_id, "/post"),
                   full.names = TRUE, pattern = ext)
      }) %>%
      # remove empty elements
      Filter(length, .)
  })

# read in files
post_files %>%
  purrr::map(function(file_ext) {
    file_ext %>%
      purrr::map_df(function(file) {
        file %>% read.table() %>% tibble::as_tibble(rownames = "row")
      }, .id = "biopsy_id")
  })

# # get indels # TODO: fix environment to allow vcfR install
# post_files[[".indel.vcf.gz$"]] %>%
#   purrr::map_df(function(file) {
#     indels <- vcfR::read.vcfR(file)@fix %>% as.data.frame()
#   }, .id = "biopsy_id")

