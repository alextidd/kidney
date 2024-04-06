#!/bin/bash

# dir
wd=/lustre/scratch126/casm/team154pc/at31/kidney/
cd $wd
data_dir=/nfs/cancer_ref01/nst_links/live/2731

# get file paths for all samples from all donors
(
    # initialise file header
    echo -e 'donor_id,sample_id,bam_id,bam_path' ;

    # we merge all bam files from the same sample 
    while read -r donor_id sample_id; do
        file_pattern=$data_dir/${sample_id}*/*.sample.merged.bam 
        
        if ls $file_pattern 1> /dev/null 2>&1; then
        
            for file in $file_pattern; do
                bam_id=$(basename $file) ; bam_id=${bam_id/.*/}
                echo -e "$donor_id,$sample_id,$bam_id,$file" ;
            done

        fi  

    done < <(sed 1d $wd/data/lcm_patient_manifest.tsv | cut -f1,2 | sort -u)
) | cat > data/sample_sheet.csv

# columns from decode       col_n   (rename?)
# tube_barcode_identifier   23       -> sample_id
# gender                    24       -> donor_gender
# donor_age_at_diagnosis    25       -> donor_age
# tissue_phenotype          26       -> sample_phenotype
# tumour or normal?         27       -> sample_tumour_status
# tissue_histology          28       -> sample_histology
# PD_ID                     30       -> biopsy_id
# creating three levels of the metadata:
# 1) donor-level, 2) sample-level, 3) biopsy-level
# all metadata columns will be prefixed with their level
(
    echo -e 'donor_id,sample_id,donor_gender,donor_age,sample_phenotype,sample_tumour_status,sample_histology,biopsy_id,biopsy_bam_file,decode_file' ;
    for file in $wd/data/decode/*.tsv; do

        sed 1d $file |
        awk -F'\t' -v OFS=',' \
            -v biopsy_bam_file=$biopsy_bam_file \
            -v decode_file=$(basename $file) \
            '{print substr($23, 1, length($23) - 1), $23, $24, $25, $26, $27, $28, $30, biopsy_bam_file, decode_file}' ;

    done
) | cat > data/decode/sample_sheet.csv



# check bams exist
while read -r biopsy_id ; do 
    file_pattern=$data_dir/${biopsy_id}*/*.sample.merged.bam
    if ! ls $file_pattern 1> /dev/null 2>&1; then 
        echo "0 files found for $biopsy_id" 
    else
        echo "$(ls $file_pattern | wc -l) files found for $biopsy_id" 
    fi  
done < <(sed 1d data/decode/sample_sheet.csv | cut -d, -f8 | sort -u)