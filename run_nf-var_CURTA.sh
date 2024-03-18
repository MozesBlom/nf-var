#!/bin/bash -l
#SBATCH -J "nf-var"
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -t 6-12:00:0
#SBATCH --mem=32000
#SBATCH --partition=xxx
#SBATCH --exclude=xxx
#SBATCH --qos standard
#SBATCH	-o nf-var_calling.r1.log

nextflow \
	run \
	/path/to/nf-var/main.nf \
	-profile curta \
	-resume
