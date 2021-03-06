#!/bin/bash
#this file downloads SRR files using ascp and saves to the output folder
#do fastqdump 
#[1]   SRR ids file 
#[2] outdir (relative path from execution dir)
#this script should be run from head node
#this file uses ascp command
#UrMi 04/26/2019

thisdir=$(pwd)
outdir="$2"
path_to_etc="/work/LAS/mash-lab/software/ascp_connect/etc"
##read list of SRR files
listfile="$1"
srr=($(awk '{print $1}' $listfile))
#do everything in outputdir
##create outdir if doesnt exist
mkdir -p $outdir
#change to outdir
cd $outdir

#download each run SRR
fail_flag=false
#declare -a failed_SRR
failed_SRR=()
for ((i=0; i<${#srr[@]};++i));
do
        this_srr=${srr[i]}
        srr_6=$(echo $this_srr | cut -c1-6)
        echo "ssh condodtn ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 500m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra $thisdir/$outdir/${this_srr}.sra"
        ascp_out=$(ssh condodtn ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 700m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra $thisdir/$outdir/${this_srr}.sra)
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
                sleep 25
                toremove=()
                for ((i=0; i<${#failed_SRR[@]};++i));
                do
                        this_srr=${failed_SRR[i]}
                        srr_6=$(echo $this_srr | cut -c1-6)
                        echo "ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 500m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra $thisdir/$outdir/${this_srr}.sra"
                        ascp_out=$(ascp -i "$path_to_etc/asperaweb_id_dsa.openssh" -k 2 -T -l 700m anonftp@ftp.ncbi.nlm.nih.gov:/sra/sra-instant/reads/ByRun/sra/SRR/${srr_6}/${this_srr}/${this_srr}.sra $thisdir/$outdir/${this_srr}.sra)
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
                        echo "FAILED to Download $this_srr ..."
                        echo "$this_srr" >> failed_srr.log
                done

fi

cd $thisdir/$outdir
#if no files were downloaded exit from here
if ls *sra 1> /dev/null 2>&1; then
        echo "Finished Downloads"

else
        #if [ ${#file_list[@]} -eq 0 ]; then
        echo "NO SRA FILE DOWNLOADED... Exiting..."
        exit 1
fi

echo "total download size:"
du -sh $PWD

#enter custom script below for downstream processing
