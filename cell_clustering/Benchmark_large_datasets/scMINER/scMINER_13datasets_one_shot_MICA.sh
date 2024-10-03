#!/bin/env bash

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
output_root='/home/xyao/MICA/One_shot_benchmark'


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

	dataset_id=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $1}')
	output_dir="${output_root}/${dataset_id}"
	mkdir -p $output_dir

	if [[ ${dataset_id} == "Covid650K" ]]; then
		n_core=40
		memory=25
		queue="superdome"
		mode="ge"
	elif [[ ${dataset_id} == "Covid97K" ]]; then
		n_core=40
		memory=6
		queue="superdome"
		mode="ge"
	elif [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=40
		memory=4
		queue="superdome"
		mode="ge"
	elif [[ ${dataset_id} == "PBMC14K" ]] || [[ ${dataset_id} == "Zeisel" ]]; then  
		n_core=40
		memory=2
		queue="superdome"
		mode="ge"
	else
		n_core=10
		memory=2
		queue="superdome"
		mode="mds"
	fi

	MICA_input=${input_dir}/$h5ad_file
	if [[ ${dataset_id} == "PBMC14K" ]]; then
		MICA_input="/home/xyao/Dataset/PBMC14K/micaInput.h5ad"
	fi	

	echo $MICA_input
	#continue

	# run MICA
	if [[ $mode == 'mds' ]];then
		MICA_CMD="mica mds -i $MICA_input -o $output_dir -nck ${n_cluster[${dataset_id}]} -nw $n_core"
	elif [[ $mode == 'ge' ]]; then # One shot benchmark for running time, use single resolution
		MICA_CMD="mica ge -i $MICA_input -o $output_dir -nw $n_core "
	else
		echo "Error, wrong mode for MICA"
	fi
	echo ${MICA_CMD} > "${output_dir}/${dataset_id}_MICA.sh"
	chmod +x "${output_dir}/${dataset_id}_MICA.sh"

	JOB="${dataset_id}_MICA"

	BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${output_dir}/${JOB}.out' -eo '${output_dir}/${JOB}.err' \
				-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} < ${output_dir}/${dataset_id}_MICA.sh"
	eval $BSUB_CMD

	echo "$JOB submitted to $queue"

	
done

