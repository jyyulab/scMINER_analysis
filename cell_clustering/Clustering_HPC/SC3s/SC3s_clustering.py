import os,sys
import logging
import time

import argparse

import scanpy as sc
import sc3s

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
#   output: clustering_label.txt, tab-separted table file with columns=['cell_id','X','Y','sc3s_label','Sil_index'])
#
#           Umap.pdf, 
#           Silouette_plot.pdf, 
#           evaluation_metrics.csv,
#           running_summary.csv
#
#           sc3s_clustered.h5ad, slots adata.uns['eval_metrics'],adata.uns['avg_sil'] 
########################################################################################################################

def main():

    # argument parser
    parser = argparse.ArgumentParser(description='sc3s cell clustering.')
    parser.add_argument('-i','--input_file',type=str,required=True,help='.h5ad format input file path')
    parser.add_argument('-o','--save_dir',type=str,required=True,help='valid output folder path.')
    parser.add_argument('-t','--input_type',type=str,default='raw',help='type of input count matrix, can be raw,tpm,log1p')
    parser.add_argument('-k','--n_cluster',type=int,default=0,help='number of clusters,single value')


    global args
    args = parser.parse_args()

    # IO
    adata = ad.read_h5ad(args.input_file)
    if not os.path.exists(args.save_dir):
        os.makedirs(args.save_dir)
        
    args.save_dir = args.save_dir.rstrip('/')

    # setup logging
    logname = args.save_dir + "/SC3s_clustering.log"
    logging.basicConfig(filename=logname,
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

    logging.info("Running sc3s cell clustering.")

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
    pca_file = args.save_dir + '/SC3s_pca_matrix.txt'
    np.savetxt(pca_file,pca_matrix,delimiter='\t')

    # neighborhood graph and Umap embedding
    #sc.pp.neighbors(adata, n_neighbors=20,n_pcs=40)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    end = time.time()
    PCA_DimRed_runtime = end - start
    logging.info('PCA and DimRed done. Runtime: {} seconds'.format(PCA_DimRed_runtime))

    # kmeans clustering
    logging.info('Begin kmeans clustering ...')
    start = time.time()


    global clustering_key # globally declear clustering_key

    clustering_key = ['sc3s_' + str(args.n_cluster)]

    # sc3s consensus k-means clustering
    sc3s.tl.consensus(adata,n_clusters = args.n_cluster)

    end = time.time()
    kmeans_runtime = end - start    
    logging.info('Consensus kmeans clustering done. Runtime: {} seconds'.format(kmeans_runtime))

    # Silhouette evluation
    logging.info('Begin Silhouette evaluation ...')
    start = time.time()

    global avg_sil, sil_index
    avg_sil,sil_index = silhouette_eval(adata)  

    end = time.time()
    Silhouette_runtime = end - start
    logging.info('Silhouette eval done. Runtime: {} seconds'.format(Silhouette_runtime))

    # clustering evaluation
    logging.info('Begin clustering evaluation ...')
    start = time.time()

    eval_metrics = clustering_eval(adata)
    eval_metrics_file = args.save_dir + "/SC3s_clustering_eval_metrics.csv"
    eval_metrics.to_csv(eval_metrics_file,sep='\t',index=False)

    end = time.time()
    clustering_eval_runtime = end - start
    logging.info('Clustering eval done. Runtime: {} seconds'.format(clustering_eval_runtime))   

    # write clustered and evaluated adata to file
    #result_file = args.save_dir + "/SC3s_clustered.h5ad"
    #adata.write(result_file)

    # UMAP
    logging.info('Begin plotting Umap ...')
    start = time.time()

    for key in clustering_key:
        fig,ax = plt.subplots(figsize=(8,8))
        sc.pl.umap(adata, color=key, 
                   title='UMAP ' + key,
                   show=False,
                   ax=ax)
        k = len(set(adata.obs[key]))
        filename =  args.save_dir + '/SC3s_Umap_' + key + '_' + str(avg_sil[key]) + '.pdf'
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
        df = pd.DataFrame(columns=['cell_id','UMAP_X','UMAP_Y','SC3s_label','Sil_index'])
        df['cell_id'] = adata.obs['cell_id'].values
        df['UMAP_X'] = X_coord
        df['UMAP_Y'] = Y_coord
        df['SC3s_label'] = adata.obs[key].values
        df['Sil_index'] = sil_index[key]

        k = len(set(adata.obs[key]))# number of cluster
        df_file = args.save_dir + "/SC3s_clustering_" + key + "_" + str(avg_sil[key]) + '.txt'

        df.to_csv(df_file,sep='\t')

    end = time.time()
    Writing_label_runtime = end - start
    logging.info('Writing clustering label done. Runtime: {} seconds'.format(Writing_label_runtime))

    logging.info('SC3s cell clustering is finished.')
        

def silhouette_eval(adata):

    # mean Silhouette coefficient of predicted clustering
    avg_sil = {} # average silhouette width
    sil_index = {} # silhouette index for all cells

    X = adata.obsm['X_pca']
    for key in clustering_key:
        pred_label = adata.obs[key].astype('category')

        k = len(set(adata.obs[key]))# number of cluster     

        if (k >= 2) & (k<= len(pred_label)-1): # where Silhouette Coefficient is defined

            avg_sil[key] = np.round(silhouette_score(X,pred_label),3)
            sil_index[key]= silhouette_samples(X,pred_label)

            fig = viz_silhouette(X, pred_label)

            sil_filename = args.save_dir + '/Silouette_' + key  + '_' + str(avg_sil[key])+'.pdf'
            fig.savefig(sil_filename,dpi=300,bbox_inches='tight')

        else:
            avg_sil[key] = None 
            sil_index[key] = [None] * len(pred_label) 

    # Silhouette Plot with true label projected to ReducedDim
    true_label = adata.obs['true_label'] # make sure this column exists
    avg_sil['true_label'] = np.round(silhouette_score(X,true_label),3)
    avg_sil_true = avg_sil['true_label']
    fig = viz_silhouette(X, true_label)
    sil_filename = args.save_dir + '/Silouette_true_label_projection_'+ str(avg_sil_true) + '.pdf'
    fig.savefig(sil_filename,dpi=300,bbox_inches='tight')

    return avg_sil,sil_index

def clustering_eval(adata):

    true_label = adata.obs['true_label']

    eval_metrics = []
    for key in clustering_key: # clustering_key globally declared before
        pred_label = adata.obs[key].values

        ARI = np.round(adjusted_rand_score(true_label,pred_label),3)
        AMI = np.round(adjusted_mutual_info_score(true_label,pred_label),3)
        NMI = np.round(normalized_mutual_info_score(true_label,pred_label),3)

        if avg_sil[key] != None:
            ASW = np.round((avg_sil[key] + 1)/2,3) # shift Silhouette width to [0,1]
            AvgBIO = np.round((ARI+NMI+ASW)/3,3)
        else:
            ASW = None
            AvgBIO = None

        new_row = [key,ARI,AMI,NMI,avg_sil[key],ASW,AvgBIO]
        eval_metrics.append(new_row)
        
    eval_metrics = pd.DataFrame(eval_metrics,columns=['k','ARI','AMI','NMI','avg_sil','ASW','AvgBIO'])    
    
    return eval_metrics


def viz_silhouette(X, cluster_labels):
    
    sample_silhouette_values = silhouette_samples(X, cluster_labels)
    n_clusters = len(set(cluster_labels))

    fig, ax = plt.subplots(figsize=(8,8))

    min_bound = np.round(min(sample_silhouette_values), 1)
    ax.set_xlim([min_bound, 1])
    # The (n_clusters+1)*10 is for inserting blank space between silhouette
    # plots of individual clusters, to demarcate them clearly.
    ax.set_ylim([0, len(X) + (n_clusters + 1) * 10])


    y_lower = 10

    for i in range(n_clusters):
        # Aggregate the silhouette scores for samples belonging to
        # cluster i, and sort them
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels.cat.codes == i]
        ith_cluster_silhouette_values.sort()
        size_cluster_i = ith_cluster_silhouette_values.shape[0]

        # plot setting
        y_upper = y_lower + size_cluster_i

        color = cm.nipy_spectral(float(i)/ n_clusters)

        ax.fill_betweenx(
            np.arange(y_lower, y_upper),
            0,
            ith_cluster_silhouette_values,
            facecolor=color,
            edgecolor=color,
            alpha=0.7,
        )

        # Label the silhouette plots with their cluster numbers at the middle
        ax.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i))
        # Compute the new y_lower for next plot
        y_lower = y_upper + 10  # 10 for the 0 samples

    ax.set_title("The silhouette plot for the various clusters.")
    ax.set_xlabel("The silhouette coefficient values")
    ax.set_ylabel("Cluster label")

    # The vertical line for average silhouette score of all the values
    silhouette_avg = silhouette_score(X, cluster_labels)
    ax.axvline(x=silhouette_avg, color="red", linestyle="--")

    ax.set_yticks([])  # Clear the yaxis labels / ticks
    ax.set_xticks([min_bound, 0, 0.5, 1])

    return fig
    

if __name__ == "__main__":
    main()


