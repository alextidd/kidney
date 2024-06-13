#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; ~/bin/jsub lsf -q week -n merge_bams -m 2g -l log "bash src/run.sh" | bsub

# merge bams
nextflow run nf-merge-normal-bams \
    --sample_sheet out/merge_normal_bams/sample_sheet.csv \
    --out_dir out/merged_normal_bams/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config 
