#!/usr/bin/env Rscript

# load libraries
library(magrittr)
library(ggplot2)

# load sample sheet
sample_sheet <-
  readr::read_csv("out/nanoseq/sample_sheet.csv")

# set order of nephron pathway
# "The correct pathway within the nephron starts with the glomerulus followed
# by the Bowman's capsule, proximal convoluted tubule, the descending loop of
# Henle, the ascending loop of Henle, distal convoluted tubule, and collecting
# duct. The glomerulus is found within the Bowman's capsule.
phenotype_order <-
  c("glomerulus" = "#66A61E",
    "PCT" = "#E6AB02",
    "DCT" = "#D95F02",
    "lymphoid_aggregate" = "#7570B3",
    "endothelial" = "#E7298A",
    "tumour" = "black")

# load biopsy metadata and reorder phenotypes
biopsy_metadata <-
  readr::read_tsv("out/metadata/biopsy_metadata.tsv") %>%
  dplyr::mutate(biopsy_phenotype = factor(biopsy_phenotype,
                                          levels = names(phenotype_order)))

# function to get files from prefix
get_files <- function(biopsy_ids, exts) {
  exts %>%
    purrr::set_names(., .) %>%
    purrr::map(function(ext) {
      biopsy_ids %>%
        purrr::set_names(., .) %>%
        purrr::map(function(biopsy_id) {
          list.files(paste0("out/nanoseq/outNextflow/", biopsy_id, "/post"),
                    full.names = TRUE, pattern = ext) 
        }) %>%
        # remove empty elements
        Filter(length, .)
    })
}

# get output tsvs
post_tsvs <-
  get_files(sample_sheet$id,
            list(".mut_burden.tsv$", ".trint_subs_obs_corrected.tsv$")) %>%
  purrr::map(function(file_ext) {
    file_ext %>%
      purrr::map_df(function(file) {
        file %>% read.table() %>% tibble::as_tibble(rownames = "row") %>%
          dplyr::filter(row != "observed")
      }, .id = "biopsy_id")
  })

# get indels
post_vcfs <-
  get_files(sample_sheet$id, list(".indel.vcf.gz$", ".muts.vcf.gz$")) %>%
  purrr::map(function(file_ext) {
    file_ext %>%
      purrr::map_df(function(file) {
        vcfR::read.vcfR(file)@fix %>%
          tibble::as_tibble() %>%
          dplyr::filter(FILTER == "PASS")
      }, .id = "biopsy_id") %>%
      # count variants per biopsy
      dplyr::group_by(biopsy_id) %>%
      dplyr::mutate(n = dplyr::n()) %>%
      dplyr::ungroup()
  })

# get all burdens
burdens <-
  post_tsvs[[".mut_burden.tsv$"]] %>%
  # standardise names
  dplyr::transmute(biopsy_id, snv_count = muts, snv_burden = burden, 
                   total, snv_total = total, indel_total = total) %>%
  dplyr::left_join(post_vcfs[[".indel.vcf.gz$"]] %>%
                     dplyr::select(biopsy_id, indel_count = n) %>%
                     dplyr::distinct()) %>%
  # calculate indel burden
  dplyr::mutate(indel_burden = indel_count / total) %>%
  # convert to long format
  tidyr::pivot_longer(
    cols = c(dplyr::ends_with("_count"), 
             dplyr::ends_with("_burden"),
             dplyr::ends_with("_total")),
    names_pattern = "([^.]+)\\_([^.]+)$",
    names_to = c("variant_type", "value_type")) %>%
  # add metadata
  dplyr::left_join(biopsy_metadata)

# plot burdens, counts, and totals vs biopsy id
burdens %>%
  ggplot(aes(x = biopsy_id, y = value, fill = biopsy_phenotype)) +
  geom_col() +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  scale_fill_manual(values = phenotype_order) +
  ggh4x::facet_nested(
    value_type ~ variant_type + biopsy_phenotype,
    scales = "free", space = "free_x") +
  theme_bw()

# plot burdens, counts, and totals vs age
burdens %>%
  dplyr::filter(!is.na(donor_age)) %>%
  ggplot(aes(x = donor_age, y = value, colour = biopsy_phenotype)) +
  geom_point() +
  scale_x_discrete(guide = guide_axis(angle = -90)) +
  scale_colour_manual(values = phenotype_order) +
  ggh4x::facet_nested(
    value_type ~ variant_type,
    scales = "free", space = "free_x") +
  theme_bw()

# set trinucleotide signature colours
colours_vec <- c("C>A" = "deepskyblue",
                 "C>G" = "black",
                 "C>T" = "firebrick2",
                 "T>A" = "gray",
                 "T>C" = "darkolivegreen3",
                 "T>G" = "rosybrown2")

# get trinucleotide substitutions
trinucs <-
  post_tsvs[[".trint_subs_obs_corrected.tsv$"]] %>%
  dplyr::rename(substitution = row) %>%
  # add metadata
  dplyr::left_join(biopsy_metadata) %>%
  dplyr::mutate(facet = paste0(stringr::str_sub(substitution, 2, 2),
                               ">",
                               stringr::str_sub(substitution, 5, 5)),
                trinuc = stringr::str_sub(substitution, 1, 3))

# function: plot trinuc distributions
plot_trinucs <- function(trinucs) {
  trinucs %>%
    ggplot(aes(x = trinuc, y = trint_subst_obs, fill = facet)) +
    geom_col() +
    scale_fill_manual(values = colours_vec) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, family = "courier"),
          legend.position = "none",
          panel.grid = element_blank()) +
    labs(x = "Substitution", y = "Observed mutation counts")
}

# plot by patient
trinucs %>%
  plot_trinucs() +
  facet_grid(donor_id ~ facet, scales = "free_x")

# plot by phenotype
trinucs %>%
  plot_trinucs() +
  facet_grid(biopsy_phenotype ~ facet, scales = "free_x")

# plot by patient x phenotype
trinucs %>%
  {split(., .$donor_id)} %>%
  purrr::map(function(df) {
    df %>%
      plot_trinucs() +
      facet_grid(biopsy_phenotype ~ facet, scales = "free_x") +
      ggtitle(unique(df$donor_id))
  })