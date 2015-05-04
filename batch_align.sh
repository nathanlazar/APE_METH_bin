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
bin_dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin
bsmooth_dir=/mnt/lustre1/users/lazar/bin/bsmooth-align-0.8.1/bin
bowtie_dir=/mnt/lustre1/users/lazar/bin/bowtie2-2.2.3
cores=24
genome=$1
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

name0=$(basename $reads_1 | cut -d'.' -f1 | sed 's/_1//')
name1=$(basename $reads_1 | cut -d'.' -f1)
name2=$(basename $reads_2 | cut -d'.' -f1)

echo $name1
echo $name2

#Make out drive if it doesn't exist
if [ ! -d $dir/$genome/ ]
  then mkdir $dir/$out_dir
fi

# Build Bowtie2 index:
# after checking whether it exists
#gen_name=`echo $genome | sed 's/.fa//'`
#if [ ! -d $dir/$gen_name.watson.1.bt2 ]
#  then $bsmooth_dir/bswc_bowtie2_index.pl \
#         --bowtie2-build=$bowtie_dir/bowtie2-build \
#         --name=$dir/$gen_name $dir/$gen_name.fa
#fi

# Split up reads 
zcat $dir/$reads_1 | split -d -a5 -l $chunk - $dir/$out_dir/$name1.split_
zcat $dir/$reads_2 | split -d -a5 -l $chunk - $dir/$out_dir/$name2.split_

# Rename split files to remove leading zeros
for f in $(ls $dir/$out_dir/$name1.split_*)
  do new=$(echo $f | sed -r 's/_0+/_/g')
  mv $f $new
done
for f in $(ls $dir/$out_dir/$name2.split_*)
  do new=$(echo $f | sed -r 's/_0+/_/g')
  mv $f $new
done
mv $dir/$out_dir/$name1.split_ $dir/$out_dir/$name1.split_0
mv $dir/$out_dir/$name2.split_ $dir/$out_dir/$name2.split_0

# Count the number of files 
proc=$(ls $dir/$out_dir/$name1.split* | wc -l)

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
    should_transfer_files = IF_NEEDED
    when_to_transfer_output = ON_EXIT
    executable='$bin_dir'/align_par.sh
    arguments = $(Process)' $genome $name1.split_'$(n)' $name2.split_'$(n)' $out_dir '$$(Cpus)
    output='$dir/$out_dir/$name0.'out.tmp.$(Process).txt
    error='$dir/$out_dir/$name0.'stderr.$(ID)
    log='$dir/$out_dir/$name0.'log.$(ID)
    request_cpus = 16
    request_memory = 15 GB
    request_disk = 1 MB
    queue '$q > $dir/$out_dir/$name0.align_par.submit.$j

  # Submit scripts to launch children for mapping
  condor_submit $dir/$out_dir/$name0.align_par.submit.$j

  # Wait for each batch to be done before launching more
  # Each job writes out 'MAPPING COMPLETE' when it's done, 
  # we wait for this as the last line of the output file
  n=0
  while (( n < q )) 
    do f=$dir/$out_dir/$name0.out.tmp.$n.txt
      while [ "$(tail -n 1 $f)" != "MAPPING COMPLETE" ]
          do sleep 15; done
      m=$(( (j * batch) + n ))
      mv $f $dir/$out_dir/$name0.out.$m.txt      
    n=$((n + 1))
    done
  j=$((j + 1))
  sleep 15 # Wait a bit to make sure everything is done
done

# Make submit script to combine results from children and 
# run sorting of evidence and tabulation of methylation counts
echo 'ID=$(Cluster).$(Process)
  should_transfer_files = IF_NEEDED
  when_to_transfer_output = ON_EXIT
  executable='$bin_dir'/combine_par.sh
  arguments = '$name0 $genome $out_dir' 
  output='$dir/$out_dir/$name0.'comb.out.$(ID)
  error='$dir/$out_dir/$name0.'comb.stderr.$(ID)
  log='$dir/$out_dir/$name0.'comb.log.$(ID)
  request_cpus = 24
  request_memory = 64 GB
  request_disk = 4 GB
  queue 1' > $dir/$out_dir/$name0.comb.submit

# Submit the script
condor_submit $dir/$out_dir/$name0.comb.submit

# Wait for everything to complete the output file from the 
# last script will end with 'COMPLETE'
while [ "$(tail -n 1 $dir/$out_dir/$name0.comb.out.*)" != "COMPLETE" ]
  do sleep 10
done

# Remove temporary files
mv $dir/$out_dir/$name0/mbias.tsv $dir/$out_dir/$name0.mbias.tsv
#rm -r $dir/$out_dir/*split*
#rm -r $dir/$out_dir/ev
#rm -r $dir/$out_dir/tsv
#rm $dir/$out_dir/*align_par.submit*

# TODO:
# dir could be an input or the current directory and 
# the path to the executables be given by flags
# Look at the memory requirements and adjust accordingly.
# gibbon_meth is hard coded as the bin directory
