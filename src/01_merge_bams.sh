#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/kidney; . ~/.bashrc; mamba activate trees; ~/bin/jsub lsf -q week -n merge_bams --ncpu 48 -m 10g -l log "bash src/01_merge_bams.sh" | bsub

# modules
module load samtools

# dir
wd=/lustre/scratch126/casm/team154pc/at31/kidney/
cd $wd
data_dir=/nfs/cancer_ref01/nst_links/live/2731
out_dir=out/merged_bams/
mkdir -p $out_dir

(
    cd $out_dir

    # we merge all bam files from the same donor 
    while read -r donor_id; do
        echo "merging bams for $donor_id"
        file_pattern=$data_dir/${donor_id}*/*.sample.merged.bam
        
        if ls $file_pattern 1> /dev/null 2>&1; then
            echo "$(ls $file_pattern | wc -l) bams found for $donor_id"
            # merge and index bams
            samtools merge -f \
                --threads 48 \
                $donor_id.bam \
                $file_pattern
            samtools index $donor_id.bam 
        else 
            echo "[WARNING: 0 bams found for $donor_id]"
        fi
        
        echo
    done < <(sed 1d $wd/data/lcm_patient_manifest.tsv | cut -f1 | sort -u)
)
