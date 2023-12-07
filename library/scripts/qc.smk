#### Sanquin Bioinformatics ####
# description: RNA-seq pipeline
# author: Margo Schuller
# date: 2021-03-10

rule fastqc:
# FastQC - A high throughput sequence QC analysis tool
# https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
    input:
        read=config['SampleDir']+'/{sample}_{paired}.fastq.gz'
    output:
        html='fastqc/{sample}_{paired}_fastqc.html',
        zip='fastqc/{sample}_{paired}_fastqc.zip'
    conda:
        'environment.yml'
    threads: 1
    resources:
        mem_mb=2000
    shell:
        'mkdir -p fastqc/tmp && \
         fastqc \
         --threads {threads} \
         --outdir fastqc \
         --dir fastqc/tmp \
         {input.read}'

use rule fastqc as fastqc_fastp with:
    input:
        read='fastp/reads/{sample}_{paired}_trimmed_paired.fastq.gz'
    output:
        html='fastqc/{sample}_{paired}_trimmed_paired_fastqc.html',
        zip='fastqc/{sample}_{paired}_trimmed_paired_fastqc.zip'


rule qualimap:
# qualimap - RNA-seq QC reports quality control metrics and bias estimations
# which are specific for whole transcriptome sequencing
# http://qualimap.conesalab.org/
    input:
        sam='star/{sample}_Aligned.sortedByCoord.out.bam',
        GTFfile=glob(config['GenomeDir']+'/*.gtf')
    output:
        counts='qualimap/{sample}.counts'
    conda:
        'environment.yml'
    threads: 1
    resources:
        mem_mb=8000
    shell:
        'unset DISPLAY && \
	 qualimap rnaseq \
         -bam {input.sam} \
         -gtf {input.GTFfile} \
         -oc  {output.counts} \
	 -pe \
         -outdir qualimap/{wildcards.sample}/ \
	 --java-mem-size={resources.mem_mb}M'