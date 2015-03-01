#!/bin/bash

# nathan dot lazar at gmail dot com

# Short script to run seqtk to trim all reads uniformly

big_bin_dir=/mnt/lustre1/users/lazar/bin/
dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/

# Usage ./trim_uniform.sh <in_dir> <out_dir> <5trim> <3trim>

# Example ./trim_uniform.sh ADAPT_TRIM TRIMMED 5 5

in_dir=$1
out_dir=$2
trim5=$3
trim3=$4

#Make out drive if it doesn't exist
if [ ! -d $dir/$out_dir ]
  then mkdir $dir/$out_dir
fi

for f in $(ls $dir/$in_dir/*.fq.gz)
# Make condor submit script for each file (lazy)
  do new=$(basename $f | sed 's/.fq.gz/_trim.fq/')
  echo 'ID=$(Cluster).$(Process)
    should_transfer_files = IF_NEEDED
    when_to_transfer_output = ON_EXIT
    executable='$big_bin_dir'/seqtk/seqtk
    arguments = trimfq -b ' $trim5' -e ' $trim3 $f'
    output='$dir/$out_dir/$new'
    error='$dir/$out_dir/$new.'stderr.$(ID)
    log='$dir/$out_dir/$new.'log.$(ID)
    request_cpus = 1
    request_memory = 12 GB
    request_disk = 1 MB
    queue 1' > $dir/$out_dir/$new.trim.submit

  # Submit scripts to launch children for mapping
  condor_submit $dir/$out_dir/$new.trim.submit
done
