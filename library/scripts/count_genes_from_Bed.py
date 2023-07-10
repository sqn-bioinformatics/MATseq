# import os
import pandas as pd

# loads dict
dictDF = pd.read_csv('library/support_files/GRCh38_PCT_ENST_geneNAME_DICT.csv')
dictDF.set_index('ENST', inplace=True)
transcript_to_gene = dictDF.to_dict()
transcript_to_gene = transcript_to_gene['gene_name']

# loads unique gene name list
uniquegenelist = open('library/support_files/GRCh38_uniquegenelist.txt').read().split('\n')
genecount = pd.DataFrame(uniquegenelist, columns=['gene'])
genecount.set_index('gene', inplace=True)

# loads bed file
bedDF=pd.read_csv(snakemake.input[0],sep='\t',header=None)
#bedDF = pd.read_csv('D:\WSL_dir\MAT AI pipeline\bed\SRR5698973.bed', sep='\t', header=None)
geneDF = bedDF[0].map(transcript_to_gene)
geneDF = geneDF.value_counts()
geneDF.rename('raw_reads', inplace=True)
genecount = genecount.join(geneDF)

genecount.fillna('0').to_csv(snakemake.output[0])
#genecount.fillna('0').to_csv('D:\WSL_dir\MAT AI pipeline\results\RPK\SRR5698973.bed')
# print(snakemake.output[0]+' done!')
