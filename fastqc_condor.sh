#!/bin/bash

# nathan dot lazar at gmail dot com

# Generate and run HTCondor script to run fastqc on all reads in a directory

# Usage ./fastqc_raw_condor.sh <in_dir> <out_dir>

dir=$(pwd)
read_dir=$dir/$1
out_dir=$2

big_bin_dir=/mnt/lustre1/users/lazar/bin

# If the out drive doesn't exist, create it
if [ ! -d $out_dir ]; then mkdir $out_dir; fi

# Create HTCondor scripts
for f in $(ls $read_dir/*.fq.gz)
 do name=$(basename $f | cut -d'.' -f1)
   echo 'should_transfer_files = IF_NEEDED
   when_to_transfer_output = ON_EXIT
   executable='$big_bin_dir'/FastQC/fastqc
   arguments = --noextract -o '$out_dir' -t 16' $f'
   output='$out_dir/$name'.out.txt
   error='$out_dir/$name'.stderr
   log='$out_dir/$name'.log
   request_cpus = 16
   request_memory = 3 GB
   request_disk = 1 MB
   queue 1' > $out_dir/$name.raw.fastqc.submit
 # Run HTCondor scripts
 condor_submit $out_dir/$name.raw.fastqc.submit
 done
