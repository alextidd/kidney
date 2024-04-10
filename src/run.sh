#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; ~/bin/jsub lsf -q week -n merge_bams  -m 2g -l log "bash nf-kidney/run.sh" | bsub

# merge bams
nextflow run nf-kidney \
    --sample_sheet data/sample_sheet.csv \
    --out_dir out \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config 