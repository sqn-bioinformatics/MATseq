import pandas as pd
import os
from load_exp_info import load_exp_info
import pickle

from scipy.stats import zscore as zscore
from sklearn.neural_network import MLPClassifier


def raw_to_zscore(DataFrame):
    outputDF=pd.DataFrame(zscore(DataFrame, axis=0)
                          ,index=DataFrame.index, columns=DataFrame.columns)
    return(outputDF)






log_filename=snakemake.log[0]

exp_info,sample_dict=load_exp_info(log_filename)


#loads MLPC model
MLPC_model=pickle.load(open('library/MLPC_model/MATseq_MLPC.model', 'rb'))

#loads list of under 10 rpm reads
under10rpm=open('library/MLPC_model/MLPC_genes_under_10.txt').read().split('\n')

#loads raw reads
reads_raw=pd.read_csv(snakemake.input[0], index_col='gene')


#creates test dataframe
x_test=raw_to_zscore(reads_raw.drop(under10rpm,axis=0)).transpose()

#runs prediction
MLPC_results=pd.DataFrame(MLPC_model.predict_proba(x_test),columns=list(MLPC_model.classes_),index=list(x_test.index))

## saves prediction as excel and csv
MLPC_results.to_excel(snakemake.output[0])
MLPC_results.to_csv(snakemake.output[1])
