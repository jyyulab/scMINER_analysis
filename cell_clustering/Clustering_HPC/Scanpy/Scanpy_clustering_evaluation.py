import os,sys
import glob
import shutil

import anndata as ad
import pandas as pd
import numpy as np

from sklearn.metrics import (silhouette_score, silhouette_samples,   
                             adjusted_rand_score, adjusted_mutual_info_score,   
                             normalized_mutual_info_score)


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


# evaluation
best_pred_label = []
for each in ordered_dataset_ids:
    true_k = true_n_cluster[each]
    #print(true_k) 
    if glob.glob(each + '/*metrics.csv'):
        metric_file = glob.glob(each + '/*metrics.csv')[0]
        eval_df = pd.read_table(metric_file)
    else:
        continue
    
    best_ari = 0
    best_pred_file = ''
    best_leiden = ''
    best_k = 100 # tag for checking correct number of cluster

    for label_file in glob.glob(each + '/*leiden*'):
        parts = label_file.split("_")  
      
        leiden = parts[2]+'_'+parts[3] 
        k_value = parts[-1].split('.')[0]
        pred_k = int(k_value[1:])

        tmp_ari = eval_df.loc[eval_df['res'] == leiden,'ARI'].values[0]

        if best_k != true_k: # number of cluster not matching true label
            if pred_k == true_k:
                best_k = pred_k
                best_ari = tmp_ari
                best_pred_file = label_file
                best_leiden = leiden
            else:
                if abs(pred_k-true_k) < abs(best_k-true_k):
                    best_ari = tmp_ari
                    best_k = pred_k
                    best_pred_file = label_file 
                    best_leiden = leiden
                    
        else:
            if pred_k == true_k and tmp_ari > best_ari:
                best_ari = tmp_ari
                best_pred_file = label_file
                best_leiden = leiden
                
    metric_row = list(eval_df.loc[eval_df['res'] == best_leiden].values[0])[1:]

    
    new_row = [each,best_pred_file] + metric_row
    #print(new_row)
    best_pred_label.append(new_row)


# write to df
best_pred_label_df = pd.DataFrame(best_pred_label)
best_pred_label_df.rename(columns = {0:'dataset_id',1:'best_pred_file',2:'pred_k',
                                    3:'ARI',4:'AMI',5:'NMI'},inplace=True)
best_pred_label_df  = best_pred_label_df .set_index('dataset_id').reindex(ordered_dataset_ids).reset_index()
best_pred_label_df.to_csv('Scanpy_clustering_evaluation.csv',sep='\t',header=True)


# Copy the selected pred label to output folder
output_dir = './best_pred_label'
os.makedirs(output_dir,exist_ok=True)
for each in ordered_dataset_ids:
    if each == 'Covid650K':
       continue

    pred_file = best_pred_label_df.loc[best_pred_label_df['dataset_id']==each,'best_pred_file'].values[0]
    pred_pdf = pred_file.replace('.txt','.pdf')
    pred_pdf = pred_pdf.replace('clustering','Umap')
    #print(pred_file)
    #print(pred_pdf)
    #break
    output_pred_file = output_dir + '/' + each + '_' +  os.path.basename(pred_file)
    output_pdf_file = output_dir + '/' + each + '_' +  os.path.basename(pred_pdf) 
    
    shutil.copy(pred_file,output_pred_file)
    shutil.copy(pred_pdf,output_pdf_file)
    #break
