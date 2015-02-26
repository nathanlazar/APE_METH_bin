#!/bin/bash

# Runs scripts necessary to align a small number of reads in parallel
# using Bsmooth (bowite2)

# nathan dot lazar at gmail dot com

# Usage:
#   align_par.sh <output_number> <species> <genome.fa>
#     <read_file1> <read_file2> <out_dir>

# Example:
#   align_par.sh \
#     1 \
#     gorilla_gorGor3/gorilla_gorGor3.fa \
#     TRIMMED/Gorilla_1.fq.gz \
#     TRIMMED/Gorilla_2.fq.gz \
#     MAPPED

#Increase the number of open files allowed
ulimit -n 4096

dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH
cores=12
bsmooth_dir=/mnt/lustre1/users/lazar/GIBBONS/VOK_GENOME/bin/bsmooth-align-0.8.1/bin
process=$1
genome=$(echo $2 | cut -d'.' -f1)
reads_1=$3
reads_2=$4
out_dir=$5
cores=$6

name0=$(echo $reads_1 | sed 's/_1//')

echo $dir
echo $cores
echo $bsmooth_dir
echo $genome
echo $reads_1
echo $reads_2

#Make drive for each set of mapped reads
###################################
if [ ! -d $out_dir/$name0 ]
  then mkdir $out_dir/$name0
fi

#Build Bowtie2 index:
####################
# should already be done

#Align reads using Bsmooth:
#Note: use -bsc to extract measurements from all Cs
#############################################
$bsmooth_dir/bswc_bowtie2_align.pl \
  --bowtie2=/mnt/lustre1/users/lazar/bin/bowtie2-2.2.3/bowtie2 \
  --samtools=/mnt/lustre1/users/lazar/bin/samtools-0.1.19 \
  --out=$out_dir/$name0/$name0.ev \
  --sam=$out_dir/$name0/$name0 \
  --metrics=$out_dir/$name0/$name0.bt2metrics \
  --stderr=$out_dir/$name0/$name0.stderr \
-- $dir/$genome \
-- $dir/$genome.fa \
-- --threads $cores \
-- $out_dir/$reads_1 \
-- $out_dir/$reads_2

echo MAPPING COMPLETE