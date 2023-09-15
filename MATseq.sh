# #!/bin/bash

# read -n1 -p 'Do you want to update the data from Genomescan server?
# press Y for yes! ' checkdata
# case $checkdata in 
# 	y|Y) echo''
# 	python library/scripts/updates_data_from_genomescan.py
# 	echo 'Genomescan update completed';;
# esac
# echo ''
# # 
# snakemake -s library/scripts/gunzip.snk -p --cores 4


# python library/scripts/removes_PCR_duplicates.py

python library/scripts/removes_PCR_duplicates_bruno.py


# snakemake -s library/scripts/compare_dudepers.snk -p --cores 4


# snakemake -s library/scripts/fastqc.snk -p --cores 4

# snakemake -s library/scripts/maping.snk -p --cores 4 --resources load=2