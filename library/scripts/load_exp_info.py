import pandas as pd
#from library.scripts.log_entry import log_entry
from log_entry import log_entry
def load_exp_info(log_filename):
    conditions_pd=pd.read_excel('experiment/exp_settings.xlsx',index_col=('field')).fillna('')

    exp_info=conditions_pd.iloc[:conditions_pd.index.get_loc('GS_sample_id')
                      ].to_dict()['values']
    
    if exp_info['GS_project_number']=='':
        log_entry('error absence of project number\n parse_exp_info.py TERMINATED',True,log_filename)
        exit()

    sample_dict=conditions_pd.iloc[conditions_pd.index.get_loc('GS_sample_id')+1:
                      ].to_dict()['values']
    
    log_entry(('Project '+str(exp_info['GS_project_number'])+' information loaded'),True,log_filename)
    
    return(exp_info,sample_dict)