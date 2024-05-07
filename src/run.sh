#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; ~/bin/jsub lsf -q week -n merge_bams -m 2g -l log "bash src/run.sh" | bsub

# merge bams
nextflow run nf-kidney \
    --sample_sheet data/sample_sheet.csv \
    --out_dir out \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config 

# run nanoseq pipeline
nextflow run ./NanoSeq/Nextflow/NanoSeq_main.nf \
    -qs 300 -profile lsf_singularity -resume \
    --jobs 100 \
    --ref /lustre/scratch126/casm/team273jn/share/pileups/reference_data/hg38/genome.fa \
    --sample_sheet data/sample_sheet_nanoseq.csv \
    --cov_Q 15 --var_b 0 --var_n 2 --var_z 15

# modules
module purge
module load singularity

head -3 data/sample_sheet_nanoseq.csv > data/sample_sheet_nanoseq_test.csv

# run nanoseq pipeline
nextflow run ./NanoSeq/Nextflow/NanoSeq_main.nf  \
  --jobs 200 -qs 20000 -profile lsf_singularity \
  -w work/ \
  --remap false \
  --ref /lustre/scratch124/casm/team78pipelines/canpipe/live/ref/Homo_sapiens/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
  --dsa_d 2 \
  --cov_Q 15 --var_b 0 --var_n 3 --var_z 12  \
  --sample_sheet data/sample_sheet_nanoseq_test.csv \
  --var_v 0.01 --var_a 50 --var_d 2 \
  --var_r 144 --var_x 8 \
  --indel_rb 2 \
  --var_q 40 \
  --outDir out/nanoseq/ \
  -resume 