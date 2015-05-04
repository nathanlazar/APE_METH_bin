#!/bin/bash

# Combine evidence from all of children mappers into one ev directory
# Sort and tabulate this

# nathan dot lazar at gmail dot com

# Usage:
#   combine_par.sh 
#     <name>
#     <genome.fa>
#     <out_dir>

# Example:
#   combine_par.sh \
#     Gorilla
#     gorilla_gorGor3/gorilla_gorGor3.fa \
#     MAPPED

#Increase the number of open files allowed
ulimit -n 4096

dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH
bsmooth_dir=/mnt/lustre1/users/lazar/GIBBONS/VOK_GENOME/bin/bsmooth-align-0.8.1/bin
cores=24
name0=$1
genome=$2
out_dir=$3

# Combine evidence files
mkdir $dir/$out_dir/$name0.ev
for d in $(ls -d $dir/$out_dir/$name0.split_*/$name0.split_*.ev)
  do for f in $(ls $d/*.ev.tsv)
    do chr=$(basename $f | sed 's/.ev.tsv//')
    if [ ! -f $dir/$out_dir/$name0.ev/$chr ]
      then touch $dir/$out_dir/$name0.ev/$chr.ev.tsv
    fi
    cat $f >> $dir/$out_dir/$name0.ev/$chr.ev.tsv
  done
done

#Sort evidence directory
$bsmooth_dir/bsev_sort.pl \
  --ev=$dir/$out_dir/$name0.ev \
  --out=$dir/$out_dir/$name0.tsv \
  --num-threads=$cores

#Tabulate
$bsmooth_dir/bsev_tabulate.pl \
  --cpg=$dir/$out_dir/$name0.cpg10 \
  --mapq-min=10 \
  --num-threads=$cores \
-- $dir/$out_dir/$name0.tsv/ \
-- $dir/$genome

mkdir $dir/$out_dir/$name0

#Create M-bias files
$bsmooth_dir/bsev_mbias.pl \
  --evidence=$dir/$out_dir/$name0.ev \
  --output=$dir/$out_dir/$name0/ \
  --num-threads=$cores

echo COMPLETE