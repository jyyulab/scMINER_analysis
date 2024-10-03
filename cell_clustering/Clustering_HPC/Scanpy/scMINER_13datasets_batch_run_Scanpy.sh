#!/bin/env bash
#BSUB -P Scanpy_Batch_true_k
#BSUB -R  "span[hosts=1]" 
#BSUB -R "rusage[mem=2GB]"
#BSUB -oo Scanpy_Batch_true_k.out -eo Scanpy_Batch_true_k.err
#BSUB -J Scanpy_Batch_true_k
#BSUB -q rhel8_superdome
#BSUB -n 1
#BSUB -N xyao@stjude.org

script_file='/home/xyao/Scanpy/scMINER_13datasets_Scanpy/Batch_true_k/Scanpy_clustering.py'

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
output_root="/home/xyao/Scanpy/scMINER_13datasets_Scanpy/Batch_true_k"

HPC_Job_dir="/home/xyao/Scanpy/scMINER_13datasets_Scanpy/Batch_true_k/HPC_Job"
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
		memory=15
		queue="rhel8_large_mem"
	elif [[ ${dataset_id} == "Covid97K" ]]; then
		n_core=40
		memory=5
		queue="rhel8_large_mem"
	elif [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=20
		memory=5
		queue="rhel8_large_mem"
	elif [[ ${dataset_id} == "PBMC14K" ]]; then
		n_core=10
		memory=2
		queue="rhel8_large_mem"
	else
		n_core=5
		memory=2
		queue="rhel8_large_mem"
	fi

	# echo $n_core
	# echo $memory
	# echo $queue
	if [[ ${dataset_id} == "Covid650K" ]]; then
	output_dir="${output_root}/${dataset_id}"

	JOB="${dataset_id}_Scanpy_clustering"

	PYTHON_CMD="python -W ignore $script_file -i "$input_dir/$h5ad_file" -o $output_dir -t $data_type" 
	#echo ${PYTHON_CMD}

	echo ${PYTHON_CMD} > "${HPC_Job_dir}/${dataset_id}_Scanpy.sh"
	chmod +x "${HPC_Job_dir}/${dataset_id}_Scanpy.sh"

	BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${HPC_Job_dir}/${JOB}.out' -eo '${HPC_Job_dir}/${JOB}.err' \
				-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} < ${HPC_Job_dir}/${dataset_id}_Scanpy.sh"
	eval $BSUB_CMD

	echo "$JOB submitted to $queue"

	fi
done

