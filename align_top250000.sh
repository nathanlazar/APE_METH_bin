#!/bin/bash

#Runs scripts necessary for aligning the top 250,000 reads 
#using Bsmooth (bowite2) so that M-bias plots can be exampined 
#to determine trimming

#nathan dot lazar at gmail dot com

# Usage:
#   align_top250000.sh \
#     <species> <read_file_1> <read_file_2> <out_dir>

# Example:
#  ./align_top250000.sh \
#     gorilla_gorGor3 \
#     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L7_1_val_1.fq.gz \
#     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L7_2_val_2.fq.gz \
#     TEST_MAP

bin_dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin
dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH

genome=$1
reads_1=$2
reads_2=$3
out_dir=$4

name0=$(basename $reads_1 | sed 's/_1_val//' | sed 's/_R1_val//' | sed 's/_1.fq.gz//')
name1=$(basename $reads_1 | cut -d'.' -f1)
name2=$(basename $reads_2 | cut -d'.' -f1)

#Get top 250,000 reads from each file 
###################################################
get_test_reads() {
  zcat $dir/$reads_1 | \
    head -n 1000000 > $dir/$out_dir/$name1.test.fq
  gzip $dir/$out_dir/$name1.test.fq

  zcat $dir/$reads_2 | \
    head -n 1000000 > $dir/$out_dir/$name2.test.fq
  gzip $dir/$out_dir/$name2.test.fq
}

#Align reads in parallel using Bsmooth:
#######################################
align() {
  $bin_dir/batch_align.sh \
    $genome \
    $out_dir/$name1.test.fq \
    $out_dir/$name2.test.fq \
    $out_dir \
    250000 \
    4
}

#Clean up temp files
#####################
clean_up() {
  rm -r $dir/$out_dir/$name0.ev
  rm -r $dir/$out_dir/$name0.cpg10
  rm -r $dir/$out_dir/$name0.tsv
}

#Main
############################
get_test_reads
align  
#clean_up
