import sys, os
import glob

import pandas as pd
import numpy as np

ordered_dataset_ids = ['Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein',
                       'Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K']

true_n_cluster = {
    'Yan':8,
    'Goolam':5,
    'Buettner':3,
    'Pollen':11,
    'Chung':4,
    'Usoskin':4,
    'Kolod':3,
    'Klein':4,
    'Zeisel':7,
    'PBMC14K':7,
    'HMC76K':20,
    'Covid97K':28,
    'Covid650K':51
}

best_pred_label = []
for each in ordered_dataset_ids:
    #print(each)
    if glob.glob(each + '/*Silouette_k*'):
        sil_file = glob.glob(each + '/*Silouette_k*')[0]
        #print(sil_file)
        pred_k = sil_file.split('_k')[1]
        pred_k = int(pred_k.split('.')[0])
    else:
        continue


    for metric_file in glob.glob(each + '/*metrics.csv'):
        if metric_file == None:
            break    
        eval_df = pd.read_table(metric_file)    

    metric_row = list(eval_df.iloc[0,:].values)[1:]   
    metric_row.insert(0,pred_k)
    metric_row.insert(0,each)
    
    best_pred_label.append(metric_row)

#print(best_pred_label)

best_pred_label_df = pd.DataFrame(best_pred_label)
best_pred_label_df.rename(columns = {0:'dataset_id',1:'pred_k',2:'ARI',
                                    3:'AMI',4:'NMI',5:'avg_sil',
                                    6:'ASW',7:'AvgBIO'},inplace=True)

best_pred_label_df  = best_pred_label_df.set_index('dataset_id').reindex(ordered_dataset_ids).reset_index()
best_pred_label_df.to_csv('scDeepCluster_clustering_evaluation.csv',sep='\t',header=True)