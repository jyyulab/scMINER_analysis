#!/bin/env bash

# conda activate mica_bin
# conda activate mica_bin2

input_dir='/home/xyao/Dataset/scMINER_cell_clustering_13datasets_h5ad_input'
output_root='/home/xyao/MICA/mica_100/bin_size_dimension_sweeping/GE_dimension_louvain_res'

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

component_range=(8 12 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96) 

for h5ad_file in $(ls $input_dir);do

	dataset_id=$(echo "${h5ad_file%.*}" |  awk -F '_' '{print $1}')
	#echo $dataset_id

	if [[ ${dataset_id} == "HMC76K" ]]; then
		n_core=20
		memory=3
		queue="rhel8_large_mem"
		mode="ge"

	elif [[ ${dataset_id} == "PBMC14K" ]]; then
		n_core=20
		memory=2
		queue="rhel8_large_mem"
		mode="ge"
	elif [[ ${dataset_id} == "Zeisel" ]]; then
		n_core=10
		memory=2
		queue="rhel8_large_mem"
		mode="ge"
	else
		n_core=4
		memory=4
		queue="rhel8_large_mem"
		mode="ge"
	fi


	# Loop through the array  
	for component in "${component_range[@]}"  
	do  
	   	#echo "$component"  
		#output_dir="${output_root}/${dataset_id}_bin_power_${bin_power}"
		output_dir="${output_root}/${dataset_id}/${dataset_id}_component_${component}"
		#echo $output_dir

		# after running prepare_MICA, use MICA_eset.log2.txt as MICA clustering input
		MICA_input="/home/xyao/MICA/scMINER_13datasets_batch_run_MICA/${dataset_id}/${dataset_id}_MICA.log2.txt"

		if [[ ${dataset_id} == 'PBMC14K' ]]; then
			MICA_input="/home/xyao/Dataset/PBMC14K_filt/PBMC14K_filt_MICAinput.txt"
		fi

		
		# run MICA
		# Single Run if exceed memory limit, remove `if` for batch_run
		if [[ $mode == 'mds' ]];then
			MICA_CMD="mica mds -i $MICA_input -o $output_dir -pn ${dataset_id} -nc ${n_cluster[${dataset_id}]} -dd ${component}"
		elif [[ $mode == 'ge' ]]; then
			MICA_CMD="mica ge -i $MICA_input -o $output_dir -pn ${dataset_id} -nw $n_core -dd ${component}\
			-nnt ann -ir 0.1 -ar 4 -ss 0.2"
		else
			echo "Error, wrong mode for MICA"
		fi

		if [[ ${dataset_id} == 'Zeisel' ]]; then
			mkdir -p $output_dir
			echo ${MICA_CMD} > "${output_dir}/${dataset_id}_MICA.sh"
			chmod +x "${output_dir}/${dataset_id}_MICA.sh"

			JOB="${dataset_id}_MICA"

			BSUB_CMD="bsub -P $JOB  -J $JOB -R \"span[hosts=1]\" -oo '${output_dir}/${JOB}.out' -eo '${output_dir}/${JOB}.err' \
						-n ${n_core} -R \"rusage[mem=${memory}GB]\" -q ${queue} < ${output_dir}/${dataset_id}_MICA.sh"
			eval $BSUB_CMD

			echo "$JOB submitted to $queue"
		
		#break	
		fi	
		#break
	done
done

