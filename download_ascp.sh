#!/bin/bash
#this file downloads SRR files using ascp
#do fastqdump 
#1   SRR ids file 
#[2] outdir (relative path from execution dir)
#this script should be run from head node
#UrMi 16/11/2017

#load all required modules
module load aspera
##get thecurrent node
thisdir=$(pwd)
outdir="$2"
path_to_etc="/opt/packages/aspera/3.6.2/connect/etc"
#scratch=/local/scratch/usingh/
##read list of SRR files
listfile="$1"
sra=($(awk '{print $1}' $listfile))
srp=($(awk '{print $2}' $listfile))
srs=($(awk '{print $3}' $listfile))
srx=($(awk '{print $4}' $listfile))
srr=($(awk '{print $5}' $listfile))
#do everything in outputdir
#full absolute path of outdir
##create outdir if doesnt exist
#make out dir in pylon5
SCRATCH="/pylon5/bi5611p/usingh"
echo "all files will be downloaded into: $SCRATCH/$outdir"
mkdir -p "$SCRATCH/$outdir"
#change to outdir
cd "$SCRATCH/$outdir"

#download each run SRR
fail_flag=false
#declare -a failed_SRR
failed_SRR=()
for ((i=0; i<${#srr[@]};++i));
do
	this_srr=${srr[i]}
	srr_6=$(echo $this_srr | cut -c1-6)
	echo "ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 500m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra ${this_srr}.sra"
	ascp_out=$(ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 700m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra ${this_srr}.sra)
	if [ $? -ne 0 ]; then
		fail_flag=true
		echo "FAILED TO DOWNLOAD " "$this_srr"
		#echo "$this_srr" >> failed_srr.log
		failed_SRR+=("$this_srr")
		continue
	fi

done

if [ ${#failed_SRR[@]} -ne 0 ]; then
	##if some files failed attempt tp download them again 5 times
	a="1"
	max_a=5
	
	while [ $a -lt 4 ]
	do
		echo "Retrying failed downloads. Attemnp $a"
		#wait 5s before attempt
		sleep 20
		toremove=()
		for ((i=0; i<${#failed_SRR[@]};++i));
		do
	        	this_srr=${failed_SRR[i]}
	        	srr_6=$(echo $this_srr | cut -c1-6)
		        echo "ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 500m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra ${this_srr}.sra"
	        	ascp_out=$(ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 700m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra ${this_srr}.sra)
		        if [ $? -eq 0 ]; then
		                #remove from failed list
				echo "Success! $this_srr"
				toremove+=("$this_srr")
		        fi
	
		done
		#remove from failed
		for r in ${toremove[@]}
		do
		   failed_SRR=("${failed_SRR[@]/$r}") #Quotes when working with strings
		done
	
		if [ ${#failed_SRR[@]} -eq 0 ]; then
			a="10"
		fi
	a=$[$a+1]
	done

fi


#if still failed_SRR is not empty write to log
#empty failed_srr.log if exists
[ -e failed_srr.log ] && rm failed_srr.log
if [ ${#failed_SRR[@]} -ne 0 ]; then
	for ((i=0; i<${#failed_SRR[@]};++i));
        	do
                	this_srr=${failed_SRR[i]}
			echo "FAILED to Download $this_srr"
			echo "$this_srr" >> failed_srr.log
		done

fi

#if no files were downloaded exit from here
#file_list=(*sra)
#r=(*sra)
#file_list=("${file_list[@]/$r}")
if ls *sra 1> /dev/null 2>&1; then
	echo "Finished Downloads"

else
	#if [ ${#file_list[@]} -eq 0 ]; then
	echo "NO SRA FILE DOWNLOADED... Exiting..."
	exit 1
fi
#after doenload call mappping script
###create a slurm script in the outdir folder
[ -e "$outdir"_slurm.sh ] && rm "$outdir"_slurm.sh
touch "$outdir"_slurm.sh
#cat $SCRATCH/base.sh > $outdir_slurm.sh

echo "#!/bin/bash" > "$outdir"_slurm.sh
jobname=$(echo "$outdir" | cut -c 9,10,11,12)
echo "#SBATCH -J "$jobname"_runmap" >> "$outdir"_slurm.sh
echo "#SBATCH -p RM" >> "$outdir"_slurm.sh
echo "#SBATCH -N 1" >> "$outdir"_slurm.sh
echo "#SBATCH -n 28" >> "$outdir"_slurm.sh
echo "#SBATCH -t 8:00:00" >> "$outdir"_slurm.sh

echo "time sh $thisdir/map.sh $PWD $PWD" >> "$outdir"_slurm.sh

echo "submitting script to sbatch..."
sbatch "$outdir"_slurm.sh
echo "total size:"
du -sh $PWD
