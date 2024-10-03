import os, sys
import glob
import numpy as np
import re
import pandas as pd

dataset = ['Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein','Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K']

run_df = []
for each in dataset:
    log_file = glob.glob(each + '/*.out')[0]
    method = log_file.split('_')[1]
    #print(log_file)
    with open(log_file,'r') as f:
        for line in f:
            if re.search('Run time',line):
                print(line)
                run_time = line.split(':')[1].strip()
                print(run_time)
            if re.search('Max Memory',line):
                print(line)
                max_mem = line.split(':')[1].strip()
                print(max_mem)

    tmp_df = [each,run_time,max_mem]
    run_df = run_df + [tmp_df]


run_df = pd.DataFrame(run_df,columns=['dataset_id','RunTime','MaxMemory'])
run_df

run_df.to_csv(method +'_RunSummary.txt',sep='\t')