#!/bin/bash

# Runs trim adapt on raw reads which
#      - removes base calls with a Phred score of 20 or lower (assuming Sanger encoding)
#      - removes any signs of the Illumina adapter sequence from the 3' end (AGATCGGAAGAGC)
#      - removes sequences that got shorter than 20 bp
#      - runs fastqc on the resulting reads
#      - stores results in <out_dir>

# Usage ./trim_adapt_condor.sh <in_dir> <out_dir>

dir=$(pwd)
read_dir=$dir/$1
out_dir=$2

big_bin_dir=/mnt/lustre1/users/lazar/bin

# If the out drive doesn't exist, create it
if [ ! -d $out_dir ]; then mkdir $out_dir; fi

# Create HTCondor scripts
for f in $(ls $read_dir/*/*_1.fq.gz)
 do name=$(basename $f | cut -d'.' -f1 | sed 's/_1//')
   f2=$(echo $f | sed 's/_1/_2/')
   echo 'should_transfer_files = IF_NEEDED
   when_to_transfer_output = ON_EXIT
   executable='$big_bin_dir'/trim_galore_zip/trim_galore
   arguments= --paired -o '$out_dir $f $f2' 
   output='$out_dir/$name'.out.txt
   error='$out_dir/$name'.stderr
   log='$out_dir/$name'.log
   request_cpus = 16
   request_memory = 3 GB
   request_disk = 1 MB
   queue 1' > $out_dir/$name.trim_galore.submit
 # Run HTCondor scripts
 condor_submit $out_dir/$name.trim_galore.submit
 done
