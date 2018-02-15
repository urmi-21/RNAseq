#!/bin/bash
#this file starts from fastq-file and finally maps reads to transcriptome
#[1] file_directory_with_SRA
#[2] outdir
#UrMi 17/11/2017

#load all required modules
#module load salmon
#required for bridges
#module load sra-toolkit
file_dir=$1
thisnode=$(/bin/hostname)
thisdir=$(pwd)
#copy index to local
#echo "copying index..."
INDEXDIR="/pylon5/bi5611p/usingh/human_index_final"
#echo "done copying."
#make list of all SRR files in input directory
file_list=($file_dir/*_pass_1.fastq)

for f in "${file_list[@]}"; do
	echo $f
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	this_name=$(awk '{split($0,a,"_"); print a[1]}' | echo "$f")
	echo $this_name
        echo $this_fname
done

