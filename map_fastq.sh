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

#for f in "${file_list[@]}"; do
#	echo $f
#	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
#	this_name=$(echo "$this_fname" | cut -d '_' -f 1)
#	echo $this_name
#       echo $this_fname
#done
cd $file_dir
echo "Running salmon..."
failed_salmon=()
for f in "${file_list[@]}"; do
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	this_name=$(echo "$this_fname" | cut -d '_' -f 1)
	echo $this_name
	echo "salmon quant -i $INDEXDIR -l A -1 "$file_dir/$this_name"_pass_1.fastq -2 "$file_dir/$this_name"_pass_2.fastq -p 28 -o "$file_dir/$this_name"_map -q"
	salmon quant -i $INDEXDIR -l A -1 "$file_dir/$this_name"_pass_1.fastq -2 "$file_dir/$this_name"_pass_2.fastq -p 28 -o "$file_dir/$this_name"_map -q -g /pylon5/bi5611p/usingh/human_genemap.txt

	if [ $? -ne 0 ]; then
                fail_flag=true
                echo "FAILED SALMON FOR " "$this_name"
                echo "$this_name" >> failed_salmon.log
                failed_salmon+=("$this_name")
                continue
        fi

	rm -f "$file_dir/$this_name"_pass_1.fastq
	rm -f "$file_dir/$this_name"_pass_2.fastq
done



echo "finished"

