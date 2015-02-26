#!/bin/bash

# Prepare and run scripts for aligning using Bsmooth in parallel

# nathan dot lazar at gmail dot com

# Usage:
#   batch_align.sh <genome.fa>
#     <read_file1> <read_file2> 
#     <out_dir> <chunk.size> <batch.size>

# Example:
#   batch_align.sh \
#     gorilla_gorGor3/gorilla_gorGor3.fa \
#     TRIMMED/Gorilla_1.fq.gz \
#     TRIMMED/Gorilla_2.fq.gz \
#     MAPPED
#     500000
#     50

dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH
bin_dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/gibbon_meth
cores=24
bsmooth_dir=/mnt/lustre1/users/lazar/GIBBONS/VOK_GENOME/bin/bsmooth-align-0.8.1/bin
genome=$(echo $1 | cut -d'.' -f1)
reads_1=$2
reads_2=$3
out_dir=$4
chunk=$5
batch=$6

echo $dir
echo $bin_dir
echo $cores
echo $bsmooth_dir
echo $genome
echo $reads_1
echo $reads_2
echo $out_dir
echo $chunk
echo $batch

name0=$(echo $reads_1 | cut -d'/' -f2 | cut -d'.' -f1 | sed 's/_1//')
name1=$(echo $reads_1 | cut -d'/' -f2 | cut -d'.' -f1)
name2=$(echo $reads_2 | cut -d'/' -f2 | cut -d'.' -f1)

echo $name1
echo $name2

#Make out drive if it doesn't exist
if [ ! -d $dir/$out_dir ]
  then mkdir $dir/$out_dir
fi

#Build Bowtie2 index:
# should already be done

# Split up reads 
zcat $reads_1 | split -d -a5 -l $chunk - $out_dir/$name1.split_
zcat $reads_2 | split -d -a5 -l $chunk - $out_dir/$name2.split_

# Rename split files to remove leading zeros
for f in $(ls $out_dir/$name1.split_*)
  do new=$(echo $f | sed -r 's/_0+/_/g')
  mv $f $new
done
for f in $(ls $out_dir/$name2.split_*)
  do new=$(echo $f | sed -r 's/_0+/_/g')
  mv $f $new
done
mv $out_dir/$name1.split_ $out_dir/$name1.split_0
mv $out_dir/$name2.split_ $out_dir/$name2.split_0

# Count the number of files 
proc=$(ls $out_dir/$name1.split* | wc -l)

# Make HTCondor submit script to spawn children mapping reads
# .. only submit up to batch jobs at a time ..
j=0
while (( j * batch < proc ))
  do
  q=$batch
  if (( (j + 1) * $batch > proc ))
    then q=$(( proc - (j * batch) ))
  fi
  echo 'ID=$(Cluster).$(Process)
    n=$$(['$j' * '$batch' + $(Process)])
    dir='$dir'
    bin_dir='$bin_dir'
    should_transfer_files = IF_NEEDED
    when_to_transfer_output = ON_EXIT
    executable=$(bin_dir)/align_par.sh
    arguments = $(Process)' $genome.fa $name1.split_'$(n)' $name2.split_'$(n)' $out_dir '$$(Cpus)
    output=$(dir)/'$out_dir/$name0.'out.tmp.$(Process).txt
    error=$(dir)/'$out_dir/$name0.'stderr.$(ID)
    log=$(dir)/'$out_dir/$name0.'log.$(ID)
    request_cpus = 16
    request_memory = 12 GB
    request_disk = 1 MB
    queue '$q > $out_dir/$name0.align_par.submit.$j

  # Submit scripts to launch children for mapping
  condor_submit $out_dir/$name0.align_par.submit.$j

  # Wait for each batch to be done before launching more
  # Each job writes out 'MAPPING COMPLETE' when it's done, 
  # we wait for this as the last line of the output file
  n=0
  while (( n < q )) 
    do f=$out_dir/$name0.out.tmp.$n.txt
      while [ "$(tail -n 1 $f)" != "MAPPING COMPLETE" ]
          do sleep 10; done
      m=$(( (j * batch) + n ))
      mv $f $out_dir/$name0.out.$m.txt      
    n=$((n + 1))
    done
  j=$((j + 1))
done

# Make submit script to combine results from children and 
# run sorting of evidence and tabulation of methylation counts
echo 'ID=$(Cluster).$(Process)
  dir='$dir'
  bin_dir='$bin_dir'
  should_transfer_files = IF_NEEDED
  when_to_transfer_output = ON_EXIT
  executable=$(bin_dir)/combine_par.sh
  arguments = '$name0 $genome.fa $out_dir' 
  output=$(dir)/'$out_dir/$name0.'comb.out.$(ID)
  error=$(dir)/'$out_dir/$name0.'comb.stderr.$(ID)
  log=$(dir)/'$out_dir/$name0.'comb.log.$(ID)
  request_cpus = 24
  request_memory = 64 GB
  request_disk = 2 GB
  queue 1' > $out_dir/$name0.comb.submit

# Submit the script
condor_submit $out_dir/$name0.comb.submit

# Wait for everything to complete the output file from the 
# last script will end with 'COMPLETE'
while [ "$(tail -n 1 $out_dir/$name0.comb.out.*)" != "COMPLETE" ]
  do sleep 10
done

# Remove temporary files
mv $out_dir/mbias.tsv $out_dir/$name0.mbias.tsv
#rm -r $out_dir/*split*
#rm -r $out_dir/ev
#rm -r $out_dir/tsv

# TODO:
# dir could be an input or the current directory and 
# the path to the executables be given by flags
# Look at the memory requirements and adjust accordingly.
# gibbon_meth is hard coded as the bin directory
