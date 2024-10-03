#!/bin/env bash

module load R/4.3.2

script_file='/home/xyao/Seurat/One_shot_benchmark/Seurat_clustering.R'

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
output_root='/home/xyao/Seurat/One_shot_benchmark'


for h5ad_file in $(ls $input_dir);do

	#echo $h5ad_file
	data_type=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $NF}')
	#echo $data_type
	dataset_id=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $1}')
	#echo $dataset_id
	#output_dir="${output_root}/${dataset_id}"
	#echo $output_dir

	if [[ ${dataset_id} == "Covid650K" ]]; then
		n_core=40
		memory=15
		queue="large_mem"
	elif [[ ${dataset_id} == "Covid97K" ]]; then
		n_core=40
		memory=6
		queue="large_mem"
	elif [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=40
		memory=4
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

	# Single Run if exceed memory limit, remove `if` for batch_run
	#case ${dataset_id} in  "PBMC14K" | "HMC76K" | "Covid97K" | "Covid650K")
	#if [[ "$dataset_id" != "PBMC14K" && "$dataset_id" != "HMC76K" && "$dataset_id" != "Covid97K" && "$dataset_id" != "Covid650K" ]]; then
		output_dir="${output_root}/${dataset_id}"
		mkdir -p $output_dir

		JOB="${dataset_id}_Seurat_clustering"

		input_file="$input_dir/$h5ad_file"
	if [[ ${dataset_id} == "PBMC14K" ]]; then
		input_file="/home/xyao/Dataset/PBMC14K/micaInput.h5ad"
	


		R_CMD="Rscript $script_file -i $input_file -o $output_dir -t $data_type" 
		#echo ${PYTHON_CMD}

		echo ${R_CMD} > "${output_dir}/${dataset_id}_Seurat.sh"
		chmod +x "${output_dir}/${dataset_id}_Seurat.sh"

		BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${output_dir}/${JOB}.out' -eo '${output_dir}/${JOB}.err' \
					-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} < ${output_dir}/${dataset_id}_Seurat.sh"
		eval $BSUB_CMD

		echo "$JOB submitted to $queue"
	fi	
	#fi
	#	;;
	#esac
	#break

		
done

