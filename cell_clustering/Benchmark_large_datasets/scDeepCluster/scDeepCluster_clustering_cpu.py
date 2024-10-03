import os,sys
default_n_threads = 40
os.environ['OPENBLAS_NUM_THREADS'] = f"{default_n_threads}"
os.environ['MKL_NUM_THREADS'] = f"{default_n_threads}"
os.environ['OMP_NUM_THREADS'] = f"{default_n_threads}"

import numpy as np
from sklearn.cluster import MiniBatchKMeans
from sklearn.cluster import KMeans

import logging
import time

scDeepCluster_Pytorch_dir = "/home/xyao/Packages/scDeepCluster_pytorch"
sys.path.append(scDeepCluster_Pytorch_dir) 


import argparse


import anndata as ad
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib import rcParams


from sklearn.metrics import (silhouette_score, silhouette_samples,   
                             adjusted_rand_score, adjusted_mutual_info_score,   
                             normalized_mutual_info_score)

import scanpy as sc

import torch
import torch.nn as nn
from torch.autograd import Variable
from torch.nn import Parameter
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader, TensorDataset
import matplotlib.pyplot as plt

from scipy.sparse import csr_matrix 

# scDeepCluster Pytorch scripts
from scDeepCluster import scDeepCluster
from single_cell_tools import *
import numpy as np
import collections
from sklearn import metrics
import h5py


from preprocess import read_dataset, normalize

#######################################################################################################################
#   input: .h5ad
#
#   output: clustering_label.txt, tab-separted table file with columns=['cell_id','X','Y','scDeepCluster_label','Sil_index'])
#
#           Umap.pdf, 
#           Silouette_plot.pdf, 
#           evaluation_metrics.csv,
#           running_summary.csv
#
#           scDeepCluster_clustered.h5ad, slots adata.uns['eval_metrics'],adata.uns['avg_sil'] 
########################################################################################################################

def main():

    # argument parser
    parser = argparse.ArgumentParser(description='scDeepCluster cell clustering.')
    parser.add_argument('-i','--input_file',type=str,required=True,help='.h5ad format input file path')
    parser.add_argument('-o','--save_dir',type=str,required=True,help='valid output folder path.')
    parser.add_argument('-t','--input_type',type=str,default='raw',help='type of input count matrix, can be raw,tpm,log1p')
    parser.add_argument('-r','--max_res',type=float,default=10,help='max resolution for leiden algorithm, default = 10')
    parser.add_argument('-s','--step_size',default=0.1,help='step size of sweeping resolution.')
   
    # original scDeepCluster Pytorch argument
    parser.add_argument('--n_clusters', default=0, type=int, 
                        help='number of clusters, 0 means estimating by the Louvain algorithm')
    parser.add_argument('--knn', default=20, type=int, 
                        help='number of nearest neighbors, used by the Louvain algorithm')
    parser.add_argument('--resolution', default=.8, type=float, 
                        help='resolution parameter, used by the Louvain algorithm, larger value for more number of clusters')
    parser.add_argument('--select_genes', default=2000, type=int, 
                        help='number of selected genes, 0 means using all genes')
    parser.add_argument('--batch_size', default=256, type=int)
    parser.add_argument('--maxiter', default=2000, type=int)
    parser.add_argument('--pretrain_epochs', default=300, type=int)
    parser.add_argument('--gamma', default=1., type=float,
                        help='coefficient of clustering loss')
    parser.add_argument('--sigma', default=2.5, type=float,
                        help='coefficient of random noise')
    parser.add_argument('--update_interval', default=1, type=int)
    parser.add_argument('--tol', default=0.1, type=float,
                        help='tolerance for delta clustering labels to terminate training stage')
    parser.add_argument('--device', default='cpu')

    parser.add_argument('--ae_weights', default=None,
                    help='file of pretrained weights, None for a new pretraining')


    global args
    args = parser.parse_args()
    print(f'n_clusters = {args.n_clusters}')

    # IO
    adata = ad.read_h5ad(args.input_file)

    if isinstance(adata.X, csr_matrix):  
        adata.X = adata.X.toarray() 
    
    adata.X = adata.X.astype('float64')

    if not os.path.exists(args.save_dir):
        os.makedirs(args.save_dir)
        
    args.save_dir = args.save_dir.rstrip('/')

    # setup logging
    logname = args.save_dir + "/scDeepCluster_clustering.log"
    logging.basicConfig(filename=logname,
                    filemode='a',
                    format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                    datefmt='%H:%M:%S',
                    level=logging.INFO)

    logging.info("Running scDeepCluster cell clustering.")

    adata.raw = adata.copy()

    if args.select_genes > 0: # default changed to 2000 for consistency with Scanpy, Seurat
        importantGenes = geneSelection(adata.X, n=args.select_genes, plot=False)
        adata = adata[:,importantGenes]


    if args.input_type == 'raw':
        print("raw count taken ")
        adata = normalize(adata,
                      size_factors=True,
                      normalize_input=True, # scaled to 0 mean and unit variance
                      logtrans_input=True)
    elif args.input_type == 'tpm':
        print(" tpm taken")
        adata = normalize(adata,
                      size_factors=False,
                      normalize_input=True,
                      logtrans_input=True)
    elif args.input_type == 'log1p':
        print("Already log1p transformed.")
        adata = normalize(adata,
              size_factors=False,
              normalize_input=True,
              logtrans_input=False)
    else:
        print("Error,data type not accepted.")
        sys.exit()


    # model training
    logging.info('Finished normalization and scaling ...')
    ae_weight_file = args.save_dir + '/AE_weights.pth.tar'

    model = scDeepCluster(input_dim=adata.n_vars, z_dim=32, 
                encodeLayer=[256, 64], decodeLayer=[64, 256], sigma=args.sigma, gamma=args.gamma, device=args.device)
    

    if args.ae_weights is None:
        start = time.time()
        model.pretrain_autoencoder(X=adata.X, X_raw=adata.raw.X, size_factor=adata.obs['size_factors'], 
                                batch_size=args.batch_size, epochs=args.pretrain_epochs, ae_weights=ae_weight_file)

        end = time.time()
        model_train_runtime = end - start 
        logging.info('Model training done. Runtime: {} seconds'.format(model_train_runtime))
    else:
        if os.path.isfile(args.ae_weights):
            print("==> loading checkpoint '{}'".format(args.ae_weights))
            checkpoint = torch.load(args.ae_weights)
            model.load_state_dict(checkpoint['ae_state_dict'])
        else:
            print("==> no checkpoint found at '{}'".format(args.ae_weights))
            raise ValueError


    # clustering
    logging.info('Begin clustering ...')
    start = time.time()
    y = adata.obs['true_label'].values
    if args.n_clusters > 0:
        y_pred, _, _, _, _ = model.fit(X=adata.X, X_raw=adata.raw.X, size_factor=adata.obs['size_factors'], n_clusters=args.n_clusters, init_centroid=None, 
                    y_pred_init=None, y=None, batch_size=args.batch_size, num_epochs=args.maxiter, update_interval=args.update_interval, tol=args.tol, save_dir=args.save_dir)
    else:
        ### estimate number of clusters by Louvain algorithm on the autoencoder latent representations
        pretrain_latent = model.encodeBatch(torch.tensor(adata.X, dtype=torch.float64)).cpu().numpy()
        adata_latent = sc.AnnData(pretrain_latent)
        sc.pp.neighbors(adata_latent, n_neighbors=args.knn, use_rep="X")
        sc.tl.louvain(adata_latent, resolution=args.resolution)
        y_pred_init = np.asarray(adata_latent.obs['louvain'],dtype=int)
        features = pd.DataFrame(adata_latent.X,index=np.arange(0,adata_latent.n_obs))
        Group = pd.Series(y_pred_init,index=np.arange(0,adata_latent.n_obs),name="Group")
        Mergefeature = pd.concat([features,Group],axis=1)
        cluster_centers = np.asarray(Mergefeature.groupby("Group").mean())
        n_clusters = cluster_centers.shape[0]
        print('Estimated number of clusters: ', n_clusters)
        y_pred, _, _, _, _ = model.fit(X=adata.X, X_raw=adata.raw.X, size_factor=adata.obs.size_factors, n_clusters=n_clusters, init_centroid=cluster_centers, 
                    y_pred_init=y_pred_init, y=None, batch_size=args.batch_size, num_epochs=args.maxiter, update_interval=args.update_interval, tol=args.tol, save_dir=args.save_dir)

    adata.obs['scDeepCluster_label'] = y_pred
    # declar global clustering key
    global clustering_key
    clustering_key = ['scDeepCluster_label']

    end = time.time()
    clustering_runtime = end - start
    logging.info('Clustering done. Runtime: {} seconds'.format(clustering_runtime))


    # write latent layer to file
    final_latent = model.encodeBatch(torch.tensor(adata.X, dtype=torch.float64).cpu()).numpy()
    np.savetxt(args.save_dir + "/final_latent_file.txt", final_latent, delimiter="\t")

    logging.info('Final latent layer extracted and saved.')


    # 20-NN graph and Umap embedding
    adata.obsm["X_scDeepCluster"] = final_latent
    sc.pp.neighbors(adata, use_rep = "X_scDeepCluster",n_neighbors=20) # consistent with scDeepCluster
    sc.tl.umap(adata)

    # # community detection with Leiden algorithm, skip this part because scDeepCluster implement its own strategy of clustering
    # logging.info('Begin Leiden community detection ...')
    # start = time.time()

    # resolution = np.round(np.arange(0.1,args.max_res,args.step_size),2).tolist()

    # global clustering_key # globally declear clustering_key
    # clustering_key = [] 

    # for res in resolution:
    #     key = 'leiden_' + str(res)
    #     clustering_key.append(key)
    #     sc.tl.leiden(adata,key_added = key, resolution = res)

    # end = time.time()
    # Leiden_runtime = end - start    
    # logging.info('Leiden community detection done. Runtime: {} seconds'.format(Leiden_runtime))

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
    eval_metrics_file = args.save_dir + "/scDeepCluster_clustering_eval_metrics.csv"
    eval_metrics.to_csv(eval_metrics_file,sep='\t',index=False)

    end = time.time()
    clustering_eval_runtime = end - start
    logging.info('Clustering eval done. Runtime: {} seconds'.format(clustering_eval_runtime))   

    # write clustered and evaluated adata to file
    result_file = args.save_dir + "/scDeepCluster_clustered.h5ad"
    adata.write(result_file)

    # UMAP
    logging.info('Begin plotting Umap ...')
    start = time.time()


    # t-SNE plot which is used in scDeepCluster tutorial page
    sc.tl.tsne(adata,use_rep='X_scDeepCluster')

    for key in clustering_key:
        fig,ax = plt.subplots(figsize=(8,8))
        adata.obs[key] = adata.obs[key].astype(str)
        #adata.obs[key] = pd.Series(adata.obs[key])
        #adata.obs[key] = adata.obs[key].astype('category')
        sc.pl.umap(adata, color=key, 
                   title='UMAP ' + key,
                   show=False,
                   ax=ax)
        k = len(set(adata.obs[key]))
        filename =  args.save_dir + '/scDeepCluster_Umap_' + 'k'+str(k)+ '.pdf'
        plt.savefig(filename,dpi=300,bbox_inches="tight")
        plt.close()

        fig,ax = plt.subplots(figsize=(8,8))
        sc.pl.tsne(adata, color=key,
                   title='t-SNE ' + key,
                   show=False,
                   ax=ax)
     
        filename =  args.save_dir + '/scDeepCluster_tSNE_' + 'k'+str(k)+ '.pdf'
        plt.savefig(filename,dpi=300,bbox_inches="tight")                
        plt.close()

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
        df = pd.DataFrame(columns=['cell_id','UMAP_X','UMAP_Y','scDeepCluster_label','Sil_index'])
        df['cell_id'] = adata.obs['cell_id'].values
        df['UMAP_X'] = X_coord
        df['UMAP_Y'] = Y_coord
        df['scDeepCluster_label'] = adata.obs[key].values
        df['Sil_index'] = sil_index[key]

        k = len(set(adata.obs[key]))# number of cluster
        df_file = args.save_dir + "/scDeepCluster_clustering_" + str(avg_sil[key]) + '_k' + str(k) + '.txt'

        df.to_csv(df_file,sep='\t')

    end = time.time()
    Writing_label_runtime = end - start
    logging.info('Writing clustering label done. Runtime: {} seconds'.format(Writing_label_runtime))

    logging.info('scDeepCluster cell clustering is finished.')
  


def silhouette_eval(adata):

    # mean Silhouette coefficient of predicted clustering
    avg_sil = {} # average silhouette width
    sil_index = {} # silhouette index for all cells

    X = adata.obsm['X_scDeepCluster']
    for key in clustering_key:
        pred_label = adata.obs[key].astype('category')

        k = len(set(adata.obs[key]))# number of cluster     

        if (k >= 2) & (k<= len(pred_label)-1): # where Silhouette Coefficient is defined

            avg_sil[key] = np.round(silhouette_score(X,pred_label),3)
            sil_index[key]= silhouette_samples(X,pred_label)

            fig = viz_silhouette(X, pred_label)

            sil_filename = args.save_dir + '/scDeepCluster_Silouette_' + 'k' + str(k) + '.pdf'
            fig.savefig(sil_filename,dpi=300,bbox_inches='tight')

        else:
            avg_sil[key] = None 
            sil_index[key] = [None] * len(pred_label) 

    # Silhouette Plot with true label projected to ReducedDim
    true_label = adata.obs['true_label'] # make sure this column exists
    avg_sil['true_label'] = np.round(silhouette_score(X,true_label),3)
    fig = viz_silhouette(X, true_label)
    sil_filename = args.save_dir + '/scDeepCluster_Silouette_true_label_projection.pdf'
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
        
    eval_metrics = pd.DataFrame(eval_metrics,columns=['res','ARI','AMI','NMI','avg_sil','ASW','AvgBIO'])    
    
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


