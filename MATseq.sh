# !/bin/bash

# read -n1 -p 'Do you want to update the data from Genomescan server?
# press Y for yes! ' checkdata
# case $checkdata in 
# 	y|Y) echo''
# 	python library/scripts/updates_data_from_genomescan.py
# 	echo 'Genomescan update completed';;
# esac
# echo ''
# # # # # 
 
snakemake -s library/scripts/1_download_files.snk  # to generate

snakemake -s library/scripts/2_unzip_files.snk -p --rerun-incomplete --cores 16

snakemake -s library/scripts/3_remove_PCR_duplicates.snk # to generate

snakemake -s library/scripts/4_trim_polyg_tails.snk -p --latency-wait 1000 --cores 16

snakemake -s library/scripts/5_quality_control.snk -p --cores 24

snakemake -s library/scripts/6_map_gene_counts.snk -p --cores 24 --resources load=2

snakemake -s library/scripts/7_classification.shk # to generate


