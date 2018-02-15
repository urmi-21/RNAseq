#!/bin/bash
#this file starts from SRA file and finally maps reads to genome after conversion and filtering
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
file_list=($file_dir/*.sra)
#file_list=(/work/LAS/mash-lab/usingh/lib_urmi/geoscripts/human_rnaseq_srr/all_sra/scripts/new_scripts/testdata_sra/SRR_LISTaa_folder/*.sra)
#run this many fastq dump in parallel
size=15
i=0
for f in "${file_list[@]}"; do
	echo "$f"
	#run fastqdump
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	echo $this_fname
	echo "fastq-dump --readids --split-files --dumpbase --skip-technical --clip --qual-filter-1 --read-filter pass --dumpbase --outdir $file_dir $f && rm -f "$f" &"
	fastq-dump --readids --split-files --dumpbase --skip-technical --clip --qual-filter-1 --read-filter pass --dumpbase --outdir $file_dir $f && rm -f "$f" & 
	#remove .sra
	#&& salmon quant -i $INDEXDIR -l A -1 "$file_dir/$this_fname"_pass_1.fastq -2 "$file_dir/$this_fname"_pass_2.fastq -p 16 -o "$file_dir/$this_fname"_map -q &
	v=$(( $(($i+1)) % $size)) 
	if [ "$v" -eq "0" ]; then
  		echo $i
		echo $v
		echo "waiting..."
		wait
	fi
	i=$(($i+1))
	
done
echo "finally waiting for fastq-dump to finish..."
wait
echo "Running salmon..."
failed_salmon=()
for f in "${file_list[@]}"; do
	this_fname=$(echo "$f" | rev | cut -d"/" -f1 | rev | cut -d"." -f1)
	echo $this_fname
	echo "salmon quant -i $INDEXDIR -l A -1 "$file_dir/$this_fname"_pass_1.fastq -2 "$file_dir/$this_fname"_pass_2.fastq -p 28 -o "$file_dir/$this_fname"_map -q"
	salmon quant -i $INDEXDIR -l A -1 "$file_dir/$this_fname"_pass_1.fastq -2 "$file_dir/$this_fname"_pass_2.fastq -p 28 -o "$file_dir/$this_fname"_map -q -g /pylon5/bi5611p/usingh/human_genemap.txt

	if [ $? -ne 0 ]; then
                fail_flag=true
                echo "FAILED SALMON FOR " "$this_fname"
                echo "$this_fname" >> failed_salmon.log
                failed_salmon+=("$this_fname")
                continue
        fi

	rm -f "$file_dir/$this_fname"_pass_1.fastq
	rm -f "$file_dir/$this_fname"_pass_2.fastq
done



echo "finished"
