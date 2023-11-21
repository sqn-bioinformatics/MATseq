import os
import re

from glob import glob
from parses_experiment_info import get_experiment_info


FILES=glob(config['SampleDir']+'/*.fastq.gz')
SAMPLES=sorted(list(set([ re.split('_R1|_R2', os.path.basename(x))[0] for x in FILES])))

exp_info, _ = get_experiment_info()

SQN_project_number = exp_info['SQN_project_number']


"""
Target rules 

"""
rule all:
    input:
    	expand(["temp/sam/{sample}.sam", 'temp/bed/{sample}.bed',
        'temp/read_counts/{sample}_genecount.csv','logs/bowtie2/{sample}.log'], sample=SAMPLES),




from glob import glob
import os
import re


# path of full files
FILES=glob(config['SampleDir']+'/*.fastq.gz')
# samples
SAMPLES=sorted(list(set([ re.split('_R1|_R2', os.path.basename(x))[0] for x in FILES])))

#### target rules ####
rule all:
    input:
        expand(['fastp/report/{sample}.html',
                'fastqc/{sample}_fastqc.html',
                'fastqc/{sample}_trimmed_paired_fastqc.html',
                'star/{sample}_Aligned.sortedByCoord.out.bam'],
                sample=SAMPLES),



"""
Load rules 

"""
include: '2_unzip_files.smk'
include: '3_control_quality.smk'
include: '4_remove_PCR_duplicates.smk'
include: '5_trim_polyg_tails.smk'
include: '6_map_gene_counts.smk'





# # 		'experiment/results/read_counts_analysis/'+SQN_project_number+'_raw_reads.xlsx', 
# 		'experiment/results/read_counts_analysis/'+SQN_project_number+'_raw_reads.csv', 
# 		'experiment/results/stats/'+SQN_project_number+'_seqstats.xlsx',
# 		'experiment/results/MLPC_prediction/'+SQN_project_number+'_MLPC_prediction.xlsx',
# 		'experiment/results/MLPC_prediction/'+SQN_project_number+'_MLPC_prediction.csv'
