#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; ~/bin/jsub lsf -q week -n merge_bams -m 2g -l log "bash src/run.sh" | bsub

# # merge bams
# nextflow run nf-kidney \
#     --sample_sheet data/sample_sheet.csv \
#     --out_dir out \
#     -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config 

# # set up sample sheet for nanoseq
# head -3 data/sample_sheet_nanoseq.csv > data/sample_sheet_nanoseq_test.csv

# run nanoseq pipeline
# increase maximum normal VAF (var_v) because query variants will be in the matched normal (due to pseudobulking)
# increase coverage (var_z) because matched normal is huge (was 12 before)
nextflow run ./NanoSeq_develop/Nextflow/NanoSeq_main.nf  \
  --jobs 200 -qs 20000 -profile lsf_singularity \
  -w work/ \
  --remap false \
  --ref /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh37d5/genome.fa \
  --dsa_d 2 \
  --cov_Q 15 --var_b 0 --var_n 3  \
  --sample_sheet data/sample_sheet_nanoseq_test.csv \
  --var_a 50 --var_d 2 \
  --var_r 144 --var_x 8 \
  --indel_rb 2 \
  --var_q 40 \
  --var_v 0.05 \
  --var_z 25 \
  --outDir out/nanoseq/ \
  -resume 