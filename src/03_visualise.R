#!/usr/bin/env Rscript

# directories
nanoseq_dir <- "out/nanoseq/outNextflow/"

# libraries
library(magrittr)
library(ggplot2)

# load metadata
biopsy_metadata <-
  readr::read_tsv("out/metadata/biopsy_metadata.tsv")

# get burdens
burdens <-
  biopsy_metadata$biopsy_id %>%
  purrr::map(function(biopsy_id) {
    file <- paste0(nanoseq_dir, "/", biopsy_id, "/post/", 
                   biopsy_id, ".mut_burden.tsv")
    if (file.exists(file)) {
      file %>%
        read.table() %>%
        tibble::as_tibble(rownames = "type") %>%
        dplyr::mutate(biopsy_id = biopsy_id)
    }
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(biopsy_metadata)

# plot burdens vs age
pdf("test.pdf", height = 5) ; burdens %>%
  dplyr::filter(!is.na(donor_age), type == "corrected") %>%
  ggplot(aes(x = donor_age, y = burden,
             group = donor_id)) +
  geom_point(aes(colour = biopsy_phenotype), size = 5) ; dev.off()

# get trinuc
trinucs <-
  biopsy_metadata$biopsy_id %>%
  purrr::map(function(biopsy_id) {
    file <- paste0(nanoseq_dir, "/", biopsy_id, "/post/", 
                   biopsy_id, ".trint_subs_obs_corrected.tsv")
    if (file.exists(file)) {
      file %>%
        read.table() %>%
        tibble::as_tibble(rownames = "substitution") %>%
        dplyr::mutate(biopsy_id = biopsy_id)
    }
  }) %>%
  dplyr::bind_rows() %>%
  dplyr::left_join(biopsy_metadata)

pdf("test.pdf", height = 20, width = 20) ; trinucs %>%
  dplyr::mutate(facet = paste0(stringr::str_sub(substitution, 2, 2),
                               ">",
                               stringr::str_sub(substitution, 5, 5))) %>%
  ggplot(aes(x = substitution, y = trint_subst_obs, fill = facet)) +
  geom_col() +
  facet_grid(biopsy_id ~ facet) +
  theme(panel.spacing = unit(0, "lines"),
        axis.text.x = element_text(size = 2)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values = colours) ; dev.off()

colours <- rep(c("deepskyblue", "black", "firebrick2", "gray", "darkolivegreen3", "rosybrown2"), each = 16)
sub_vec <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
ctx_vec <- paste(rep(c("A", "C", "G", "T"), each = 4), rep(c("A", "C", "G", "T"), times = 4), sep = "-")
full_vec <- paste(rep(sub_vec, each = 16), rep(ctx_vec, times = 6), sep = ",")
xstr <- paste(substr(full_vec, 5, 5), substr(full_vec, 1, 1), substr(full_vec, 7, 7), sep = "")
ordered_names <- paste(xstr, ">", rep(c("A", "G", "T", "A", "C", "G"), each = 16), sep = "")
# tmp_ <- table(unique_variants[which(unique_variants$ismasked == 0), "pyrsub"])
# tri_obs = as.vector(tmp_)
# names(tri_obs) = names(tmp_)
# tri_obs[setdiff(ordered_names, names(tri_obs))] = 0
# tri_obs = tri_obs[ordered_names]

colours_vec <- c("C>A" = "deepskyblue",
                 "C>G" = "black",
                 "C>T" = "firebrick2",
                 "T>A" = "gray",
                 "T>C" = "darkolivegreen3",
                 "T>G" = "rosybrown2")

p <- list()
pdf(width = 9, height = 3.5, file = "test.pdf", onefile = T)

trinucs %>%
  {split(., .$donor_id)} %>%
  purrr::map(function(df) {

    # summarise by phenotype
    df %>%
      dplyr::group_by(donor_id, biopsy_phenotype, substitution) %>%
      dplyr::summarise(trint_subst_obs = sum(trint_subst_obs)) %>%
      {split(., .$biopsy_phenotype)} %>%
      purrr::map(function(df2) {
        id <- paste0(unique(df2$donor_id), "_", unique(df2$biopsy_phenotype))
        tri_obs <- df2$trint_subst_obs
        names(tri_obs) <- df2$substitution

        # observed counts
        y <- tri_obs
        maxy <- max(y)
        h <- barplot(y, las = 2, col = colours, border = NA, ylim = c(0, maxy * 1.5),
                    space = 1, cex.names = 0.6, #names.arg = xstr,
                    ylab = "Observed mutation counts",
                    main = id)
        for (j in 1:length(sub_vec)) {
          xpos <- h[c((j - 1) * 16 + 1, j * 16)]
          rect(xpos[1] - 0.5, maxy * 1.2,
              xpos[2] + 0.5, maxy * 1.3, 
              border = NA,
              col = colours[j * 16])
          text(x = mean(xpos), y = maxy * 1.3, pos = 3, label = sub_vec[j])
        }

        # with ggplot
        df2_p <-
          df2 %>%
          dplyr::mutate(
            colours = colours,
            trinuc = stringr::str_sub(substitution, 1, 3),
            facet = factor(
              paste0(stringr::str_sub(substitution, 2, 2), ">",
              stringr::str_sub(substitution, 5, 5)), levels = sub_vec))
        p[[id]] <<-
          ggplot(df2_p, aes(x = trinuc, y = trint_subst_obs, fill = facet)) +
          geom_col() +
          scale_fill_manual(values = colours_vec) +
          theme_minimal() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1),
                legend.position = "none") +
          labs(title = id, x = "Substitution", y = "Observed mutation counts") +
          facet_grid(~ facet, scales = "free_x")
      })
  })

dev.off()

pdf("test.pdf", width = 15) ; purrr::walk(p, print) ; dev.off()
