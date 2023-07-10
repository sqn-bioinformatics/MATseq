import pandas as pd
import os


log_filename=snakemake.log[0]

##loads dataframe with deduping stats
## the path is given by the snakemake pipeline

deduping_stats=pd.read_csv(snakemake.input.dedupstats, index_col='sample')

## creates list with sample references

bowtie_log_path=snakemake.input.seqlogs #path to individual bowtie logs

#extracts sample name from the log path.
sample_list=[]
for sample_path in bowtie_log_path:
      sample_list.append(sample_path.replace('logs/bowtie2/','').replace('.log',''))

    
    
## parses deduped stats
stats=[]

for sample in sample_list:
    temp=deduping_stats.loc[sample].set_index('run')
    unique_reads=temp.loc['merged']['unique']
    temp.drop('merged',axis=0,inplace=True)
    initial_reads=temp['initial'].sum()
    stats.append(dict(zip(['intial_reads','unique_reads'],[initial_reads,unique_reads])))

stats=pd.DataFrame(stats)
stats.index=sample_list

### loads mapping information and extracts mapped reads
mapped_reads=[]

for i in range(len(sample_list)):

    sample=sample_list[i]
    log=open(bowtie_log_path[i]).readlines()
    
    #removes lines until the mapping report starts on log[0]
    try:
        while not log[0].startswith(str(stats.loc[sample]['unique_reads'])):
            log.remove(log[0])
    except:
        print(sample+' log is not what it should be!!')

        
    mapped_reads.append(int(log[3].replace(' ','').split('(')[0])+
                       int(log[4].replace(' ','').split('(')[0]))
    
    

stats['mapped']=mapped_reads

stats.index.name='GS_sample_code'

### loads sample name information
from load_exp_info import load_exp_info

exp_info,sample_dict=load_exp_info(log_filename)


###changes index names
oldnames=stats.index.to_list()
newnames=[]
for name in oldnames:
    newnames.append(sample_dict[name])

stats['sample name']=newnames


stats.to_excel(snakemake.output[0])
