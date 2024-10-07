# scMINER_analysis
Codes used for benchmarking scMINER in cell clustering, network inference and hidden driver identification

## network_inference - evaluation of TF and SIG network
network_inference/ATAC_GRN_evaluation - evaluation of TF network based on ATAC-Seq data. 
Required input: expression tables, networks inferred with scMINER, GENIE3, GRNBoost2 and PIDC, list of TFs and SIGs

network_inference/CROP_Seq_GRN_evaluation - evaluation of TF and SIG networks based on CROP-Seq data. 
Required input: expression tables (single-cell and metacell), networks inferred with scMINER, GENIE3, GRNBoost2 and PIDC, list of TFs and SIGs, activity matrix.

## activity_estimation - activity-based analysis with scMINER and SCENIC
activity_estimation/PBMC - analysis of PBMC dataset. 

activity_estimation/TEX - analysis of exhausted T-cells dataset. 

activity_estimation/Tregs - analysis of regulatory T-cells dataset. 

Required input: Expression table, activity matrix computed with scMINER and with SCENIC, clustering results obtained with MICA.
