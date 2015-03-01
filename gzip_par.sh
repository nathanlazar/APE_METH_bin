#!/bin/bash

# nathan dot lazar at gmail dot com

# script to gzip files in parallel on HTCondor

# Usage: ./gzip_par.sh <dir>

dir=$1

for f in $(ls $dir/*.fq)
  do echo 'should_transfer_files = IF_NEEDED
  when_to_transfer_output = ON_EXIT
  executable=/usr/bin/gzip 
  arguments= '$f'
  output= /dev/null
  error= /dev/null
  log=/dev/null
  request_cpus = 1
  request_memory = 1 GB
  request_disk = 1 MB
  queue 1' > $f.gzip.submit
  # Run HTCondor scripts
  condor_submit $f.gzip.submit
done
