#!/bin/env bash

script_file='/home/xyao/scDeepCluster/One_shot_benchmark/scDeepCluster_clustering_cpu.py'

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
output_root='/home/xyao/scDeepCluster/One_shot_benchmark'

declare -A n_cluster
n_cluster["Yan"]=8
n_cluster["Goolam"]=5
n_cluster["Buettner"]=3
n_cluster["Pollen"]=11
n_cluster["Chung"]=4
n_cluster["Usoskin"]=4
n_cluster["Kolod"]=3
n_cluster["Klein"]=4
n_cluster["Zeisel"]=7
n_cluster["PBMC14K"]=7
n_cluster["HMC76K"]=20
n_cluster["Covid97K"]=28
n_cluster["Covid650K"]=34



for h5ad_file in $(ls $input_dir);do

	#echo $h5ad_file
	data_type=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $NF}')
	#echo $data_type
	dataset_id=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $1}')
	#echo $dataset_id
	output_dir="${output_root}/${dataset_id}"
	#echo $output_dir

	if [[ ${dataset_id} == "Covid650K" ]]; then
		n_core=40
		memory=25
		queue="large_mem"
	elif [[ ${dataset_id} == "Covid97K" ]]; then
		n_core=40
		memory=10
		queue="large_mem"
	elif [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=40
		memory=6
		queue="large_mem"
	elif [[ ${dataset_id} == "PBMC14K" ]]; then
		n_core=40
		memory=2
		queue="large_mem"
	else
		n_core=40
		memory=2
		queue="large_mem"
	fi

	# echo $n_core
	# echo $memory
	# echo $queue

	# For Covid650K
	#ae_weights='/home/xyao/scDeepCluster/scMINER_13datasets_batch_run_scDeepCluster/Covid650K/AE_weights.pth.tar'

	# Single Run if exceed memory limit, remove `if` for batch_run
	#case ${dataset_id} in "Covid97K" | "Covid650K")

	#if [[ "$dataset_id" != "PBMC14K" && "$dataset_id" != "HMC76K" && "$dataset_id" != "Covid97K" && "$dataset_id" != "Covid650K" ]]; then
	output_dir="${output_root}/${dataset_id}"
	mkdir -p $output_dir

	JOB="${dataset_id}_scDeepCluster_clustering_Kmeans"

	input_file="$input_dir/$h5ad_file"
	if [[ "$dataset_id" == "PBMC14K" ]]; then
		input_file="/home/xyao/Dataset/PBMC14K/micaInput.h5ad"

	PYTHON_CMD="python -W ignore $script_file -i $input_file -o $output_dir -t $data_type \
			--device 'cpu' --n_cluster ${n_cluster[${dataset_id}]}" 
	#echo ${PYTHON_CMD}

	echo ${PYTHON_CMD} > "${output_dir}/${dataset_id}_scDeepCluster.sh"
	chmod +x "${output_dir}/${dataset_id}_scDeepCluster.sh"

	BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${output_dir}/${JOB}.out' -eo '${output_dir}/${JOB}.err' \
				-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} \
				< ${output_dir}/${dataset_id}_scDeepCluster.sh"
	eval $BSUB_CMD

	echo "$JOB submitted to $queue"

	fi
	#;;
	#esac 
	#break
		
done
