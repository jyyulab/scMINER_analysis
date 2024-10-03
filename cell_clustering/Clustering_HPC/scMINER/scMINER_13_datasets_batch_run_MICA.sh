#!/bin/env bash
#BSUB -P MICA_Batch_Run
#BSUB -R  "span[hosts=1]" 
#BSUB -R "rusage[mem=2GB]"
#BSUB -oo MICA_Batch_Run.out -eo MICA_Batch_Run.err
#BSUB -J MICA_Batch_Run
#BSUB -q rhel8_large_mem
#BSUB -n 1
#BSUB -N xyao@stjude.org

#conda activate MICA

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
#output_root='/home/xyao/MICA/scMINER_13datasets_batch_run_MICA'

#HPC_Job_dir='/home/xyao/MICA/scMINER_13datasets_batch_run_MICA/HPC_Job'

# for Covid650K only
output_root='/home/xyao/MICA/Covid650K_log1p'
HPC_Job_dir='/home/xyao/MICA/Covid650K_log1p'

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
	#data_type=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $NF}')
	#echo $data_type
	dataset_id=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $1}')
	#echo $dataset_id
	output_dir="${output_root}/${dataset_id}"
	#echo $output_dir

	if [[ ${dataset_id} == "Covid650K" ]]; then
		n_core=40
		memory=15
		queue="rhel8_large_mem"
		mode="ge"
	elif [[ ${dataset_id} == "Covid97K" ]]; then
		n_core=10
		memory=6
		queue="rhel8_large_mem"
		mode="ge"
	elif [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=10
		memory=4
		queue="rhel8_large_mem"
		mode="ge"
	elif [[ ${dataset_id} == "PBMC14K" ]]; then
		n_core=10
		memory=2
		queue="rhel8_large_mem"
		mode="ge"
	elif [[ ${dataset_id} == "Zeisel" ]]; then
		n_core=5
		memory=2
		queue="rhel8_superdome"
		mode="ge"
	else
		n_core=2
		memory=4
		queue="rhel8_superdome"
		mode="mds"
	fi

	
	# after running prepare_MICA, use MICA_eset.log2.txt as MICA clustering input
	if [[ ${dataset_id} == "Covid650K" ]];then
		#MICA_input="${output_dir}/${dataset_id}_MICA.log2.txt"
		MICA_input="/home/xyao/MICA/scMINER_13datasets_batch_run_MICA/Covid650K/Covid650K_MICA.log2.txt"
	#MICA_input="${output_dir}/${dataset_id}_MICA.log2.txt"
	cwl_json="/home/xyao/MICA/scMINER_13datasets_batch_run_MICA/config_cwlexec.json"
	# run MICA
	# Single Run if exceed memory limit, remove `if` for batch_run
	if [[ $mode == 'mds' ]];then
		MICA_CMD="mica mds -i $MICA_input -o $output_dir -pn ${dataset_id} -nc ${n_cluster[${dataset_id}]}"
	elif [[ $mode == 'ge' ]]; then
		MICA_CMD="mica ge -i $MICA_input -o $output_dir -pn ${dataset_id} -ar 4.0 -nw $n_core"
	else
		echo "Error, wrong mode for MICA"
	fi
		
	echo ${MICA_CMD} > "${HPC_Job_dir}/${dataset_id}_MICA.sh"
	chmod +x "${HPC_Job_dir}/${dataset_id}_MICA.sh"

	JOB="${dataset_id}_MICA"

	BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${HPC_Job_dir}/${JOB}.out' -eo '${HPC_Job_dir}/${JOB}.err' \
				-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} < ${HPC_Job_dir}/${dataset_id}_MICA.sh"
	eval $BSUB_CMD

	echo "$JOB submitted to $queue"

	
	#break	
	fi	
done

