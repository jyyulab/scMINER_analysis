#!/bin/env bash
#BSUB -P MICA100
#BSUB -R "rusage[mem=2GB]"
#BSUB -R "span[hosts=1]"
#BSUB -oo MICA.out -eo MICA.err
#BSUB -J MICA100
#BSUB -q superdome
#BSUB -n 40


indir="."
cd $indir


MICA_input="../PBMC14K_act_MICA_input.h5ad"

output_dir="./"
mica louvain -i $MICA_input -o $output_dir -minr 0.5 -maxr 4.0 -ss 0.1 -nw 40 


echo "MICA job completed."
