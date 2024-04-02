from glob import glob
import os
import re

files = glob(config["SampleDir"]+'/*.fastq.gz')

file_names = sorted(list(set([re.split('_R1|_R2', os.path.basename(file))[0] for file in files])))
 
