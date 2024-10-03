import os,sys
import logging
import time

import argparse

import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams


from sklearn.metrics import (silhouette_score, silhouette_samples,   
                             adjusted_rand_score, adjusted_mutual_info_score,   
                             normalized_mutual_info_score)



#######################################################################################################################
#   input: .h5ad
#
#   output: clustering_label.txt, tab-separted table file with columns=['cell_id','X','Y','Scanpy_label','Sil_index'])
#
#           Umap.pdf, 
#           Silouette_plot.pdf, 
#           evaluation_metrics.csv,
#           running_summary.csv
#
#           Scanpy_clustered.h5ad, slots adata.uns['eval_metrics'],adata.uns['avg_sil'] 
########################################################################################################################

def main():

    # argument parser
    parser = argparse.ArgumentParser(description='Scanpy cell clustering.')
    parser.add_argument('-i','--input_file',type=str,required=True,help='.h5ad format input file path')
    parser.add_argument('-o','--save_dir',type=str,required=True,help='valid output folder path.')
    parser.add_argument('-t','--input_type',type=str,default='raw',help='type of input count matrix, can be raw,tpm,log1p')
    parser.add_argument('-r','--max_res',type=float,default=2,help='max resolution for leiden algorithm, default = 4')
    parser.add_argument('-s','--step_size',default=0.01,help='step size of sweeping resolution.')

    global args
    args = parser.parse_args()

    # IO
    adata = ad.read_h5ad(args.input_file)
    if not os.path.exists(args.save_dir):
        os.makedirs(args.save_dir)
        
    args.save_dir = args.save_dir.rstrip('/')

    # setup logging
    logname = args.save_dir + "/Scanpy_clustering.log"
    logging.basicConfig(filename=logname,
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

    logging.info("Running Scanpy cell clustering.")

    # normalization and log transformation
    #adata.raw = adata
    if args.input_type == 'raw':
        sc.pp.normalize_total(adata)
        sc.pp.log1p(adata)    
    elif args.input_type == 'tpm':
        sc.pp.log1p(adata)
    elif args.input_type == 'log1p':
        print("Already log1p transformed.")
    else:
        print("Error,data type not accepted.")
        sys.exit()

    # select highly variable genes
    logging.info('Read gene count matrix and do PCA DimRed...')
    start = time.time()
    
    sc.pp.highly_variable_genes(adata,n_top_genes=2000,subset=True)
    
    # scale each gene to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata,max_value=10)

    # PCA for dimension reduction
    sc.tl.pca(adata, svd_solver='arpack')

    pca_matrix = adata.obsm['X_pca']
    pca_file = args.save_dir + '/Scanpy_pca_matrix.txt'
    np.savetxt(pca_file,pca_matrix,delimiter='\t')


    # neighborhood graph and Umap embedding
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    end = time.time()
    PCA_DimRed_runtime = end - start
    logging.info('PCA and DimRed done. Runtime: {} seconds'.format(PCA_DimRed_runtime))

    # community detection with Leiden algorithm
    logging.info('Begin Leiden community detection ...')
    start = time.time()

    resolution = np.round(np.arange(0.01,args.max_res,args.step_size),2).tolist()

    global clustering_key # globally declear clustering_key
    clustering_key = [] 

    for res in resolution:
        key = 'leiden_' + str(res)
        clustering_key.append(key)
        sc.tl.leiden(adata,key_added = key, resolution = res)

    end = time.time()
    Leiden_runtime = end - start    
    logging.info('Leiden community detection done. Runtime: {} seconds'.format(Leiden_runtime))

    # clustering evaluation
    logging.info('Begin clustering evaluation ...')
    start = time.time()

    eval_metrics = clustering_eval(adata)
    eval_metrics_file = args.save_dir + "/Scanpy_clustering_eval_metrics.csv"
    eval_metrics.to_csv(eval_metrics_file,sep='\t',index=False)

    end = time.time()
    clustering_eval_runtime = end - start
    logging.info('Clustering eval done. Runtime: {} seconds'.format(clustering_eval_runtime))   

    # write clustered and evaluated adata to file
    result_file = args.save_dir + "/Scanpy_clustered.h5ad"
    adata.write(result_file)

    # UMAP
    logging.info('Begin plotting Umap ...')
    start = time.time()

    for key in clustering_key:
        fig,ax = plt.subplots(figsize=(8,8))
        sc.pl.umap(adata, color=key, 
                   title='UMAP Scanpy ' + key,
                   show=False,
                   ax=ax)
        k = len(set(adata.obs[key]))
        filename =  args.save_dir + '/Scanpy_Umap_' + key + '_k'+str(k)+ '.pdf'
        plt.savefig(filename,dpi=300,bbox_inches="tight")

    end = time.time()
    Umap_runtime = end - start
    logging.info('Umap plotting done. Runtime: {} seconds'.format(Umap_runtime))


    # save clustering label to .txt file
    logging.info('Writing clustering labels ...')
    start = time.time() 

    umap_coords = adata.obsm['X_umap']
    X_coord = umap_coords[:,0]
    Y_coord = umap_coords[:,1]

    for key in clustering_key:
        pred_label = adata.obs[key].values
        df = pd.DataFrame(columns=['cell_id','UMAP_X','UMAP_Y','Scanpy_label'])
        df['cell_id'] = adata.obs['cell_id'].values
        df['UMAP_X'] = X_coord
        df['UMAP_Y'] = Y_coord
        df['Scanpy_label'] = adata.obs[key].values

        k = len(set(adata.obs[key]))# number of cluster
        df_file = args.save_dir + "/Scanpy_clustering_" + key + '_k' + str(k) + '.txt'

        df.to_csv(df_file,sep='\t')

    end = time.time()
    Writing_label_runtime = end - start
    logging.info('Writing clustering label done. Runtime: {} seconds'.format(Writing_label_runtime))

    logging.info('Scanpy cell clustering is finished.')
        



def clustering_eval(adata):

    true_label = adata.obs['true_label']

    eval_metrics = []
    for key in clustering_key: # clustering_key globally declared before
        pred_label = adata.obs[key].values

        pred_k = len(set(pred_label))
        ARI = np.round(adjusted_rand_score(true_label,pred_label),3)
        AMI = np.round(adjusted_mutual_info_score(true_label,pred_label),3)
        NMI = np.round(normalized_mutual_info_score(true_label,pred_label),3)

        new_row = [key,pred_k,ARI,AMI,NMI]
        eval_metrics.append(new_row)
        
    eval_metrics = pd.DataFrame(eval_metrics,columns=['res','pred_k','ARI','AMI','NMI'])    
    
    return eval_metrics


if __name__ == "__main__":
    main()


