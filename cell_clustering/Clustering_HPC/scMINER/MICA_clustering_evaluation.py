import os,sys
import glob

import anndata as ad
import pandas as pd
import numpy as np

from sklearn.metrics import (silhouette_score, silhouette_samples,   
                             adjusted_rand_score, adjusted_mutual_info_score,   
                             normalized_mutual_info_score)


# to be run within MICA/scMINER_13datasets_batch_run_MICA/

h5ad_root = '../../Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
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

MDS_datasets = ['Yan','Goolam','Buettner','Pollen','Chung','Usoskin','Kolod','Klein']
GE_datasets = ['Zeisel','PBMC14K','HMC76K','Covid97K','Covid650K']

# evaluation
eval_metrics = []
for dataset in ordered_dataset_ids:
    #if dataset == 'Covid650K':
    #    continue
        
    #print(dataset)
    if glob.glob(h5ad_root + '/' + dataset + '*.h5ad'):
        h5ad_file = glob.glob(h5ad_root + '/' + dataset + '*.h5ad')[0]
    else:
        continue

    adata = ad.read_h5ad(h5ad_file)
    #print(adata.obs.columns)
    true_label_df = pd.DataFrame(adata.obs[['cell_id','true_label']])
    #print(true_label_df)

    #break    
    # evaluation for MDS 
    if dataset in MDS_datasets:
        print(dataset)
        
        pred_label_file = glob.glob(dataset + '/' + 'clustering_umap'+ '*.txt')[0]
        pred_df = pd.read_table(pred_label_file)
        #print(pred_df)

        merge_df = pred_df.merge(true_label_df,left_on = 'ID',right_on='cell_id', how='inner') 
        #print(merge_df)

        ARI = np.round(adjusted_rand_score(merge_df['true_label'],merge_df['label']),3)
        AMI = np.round(adjusted_mutual_info_score(merge_df['true_label'],merge_df['label']),3)  
        NMI = np.round(normalized_mutual_info_score(merge_df['true_label'],merge_df['label']),3)

        # sil based metrics not available now
        pred_k = len(set(merge_df['label']))
        new_row = [dataset,pred_k,ARI,pred_label_file,AMI,NMI]
        eval_metrics.append(new_row) 

    # evaluation for GE
    if dataset in GE_datasets:   
        print(dataset)
        
        best_ari = 0
        best_pred_file = ''
        best_leiden = ''
        best_k = 100 # number of cluster for the final selected pred label
        
        true_k =  true_n_cluster[dataset]
        pred_AMI = 0
        pred_NMI = 0
            
        # check silhouette files to extract pred_k (number of cluster in pred labels)
        sil_files = glob.glob(dataset + '/' + 'silhouette'+ '*.pdf')
        #print(sil_files)
        for each in sil_files:
            sil_file_name = os.path.basename(each)
            sil_file_name = os.path.splitext(sil_file_name)[0]
            pred_k = int(sil_file_name.split('_')[2])
            res = sil_file_name.split('_')[3]
            #print(res)
            # check corresponding clustering labels for ARI 
            pred_label_file = glob.glob(dataset + '/*' + res + '.txt')[0]
            #print(pred_label_file)
            pred_df = pd.read_table(pred_label_file)
            #print(pred_df)

            merge_df = pred_df.merge(true_label_df,left_on = 'ID',right_on='cell_id', how='inner')    
            #print(merge_df)
            #break

            ARI = np.round(adjusted_rand_score(merge_df['true_label'],merge_df['label']),3)
            AMI = np.round(adjusted_mutual_info_score(merge_df['true_label'],merge_df['label']),3)  
            NMI = np.round(normalized_mutual_info_score(merge_df['true_label'],merge_df['label']),3)

            # select pred label 
            if best_k != true_k: # n_cluster not matching true n_cluster
                if pred_k == true_k:
                    best_k = pred_k
                    best_ari = ARI
                    best_pred_file = pred_label_file
                    best_leiden = 'leiden_' + res 

                    pred_AMI = AMI
                    pred_NMI = NMI
                else:
                    if abs(pred_k-true_k) < abs(best_k-true_k):
                        best_ari = ARI
                        best_k = pred_k
                        best_pred_file = pred_label_file 
                        best_leiden = 'leiden_' + res 
                        pred_AMI = AMI
                        pred_NMI = NMI
                        
            else: # n_cluster matches true n_cluster already
                if pred_k == true_k and ARI > best_ari:
                    best_ari = ARI
                    best_pred_file = pred_label_file
                    best_leiden = 'leiden_' + res 
                    pred_AMI = AMI
                    pred_NMI = NMI

        new_row = [dataset,best_k,best_ari,best_pred_file,pred_AMI,pred_NMI]
        eval_metrics.append(new_row) 
   # break


# write evaluation metrics to csv file
eval_metrics = pd.DataFrame(eval_metrics)
eval_metrics.rename(columns = {0:'dataset_id',1:'n_cluster',2:'ARI',
                            3:'pred_label_file',4:'AMI',5:'NMI'},inplace=True)
eval_metrics.to_csv("scMINER_clustering_evaluation.csv",sep='\t')