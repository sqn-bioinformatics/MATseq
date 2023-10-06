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
# snakemake -s library/scripts/gunzip.snk -p --rerun-incomplete --cores 4 

python library/scripts/removes_PCR_duplicates.py

# snakemake -s library/scripts/fastqc.snk -p --cores 4

snakemake -s library/scripts/fastqc_separate_runs.snk -p --cores 4

# snakemake -s library/scripts/mapping.snk -p --cores 4 --resources load=2