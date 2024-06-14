#!/usr/bin/env Rscript
#----------Corrected for WGS, SBS---------

# read args
args <- commandArgs(trailingOnly = TRUE)
mut_burden_corrected <- read.table(args[1], check.names = FALSE, header = FALSE)
sigprofiler_matrix <- read.table(args[2], check.names = FALSE, header = TRUE)
output_file <- args[3]

selected <- sigprofiler_matrix[, -1]
for (i in 1:ncol(selected)) {
  patient <- substr(colnames(selected)[i], 1, 15)
  selected[, i] <-
    selected[, i] * 6e9 *
    mut_burden_corrected$V5[grep(patient, mut_burden_corrected$V1)] /
    mut_burden_corrected$V3[grep(patient, mut_burden_corrected$V1)]
}
selected <- round(selected)
sigprofiler_matrix_corrected <- cbind(sigprofiler_matrix[, 1], selected)
colnames(sigprofiler_matrix_corrected)[1] <- "MutationType"

write.table(sigprofiler_matrix_corrected, output_file,
            quote = FALSE, sep = "\t", row.names = FALSE)
