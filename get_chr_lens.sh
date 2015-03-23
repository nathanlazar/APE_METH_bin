#!/bin/bash

# nathan dot lazar at gmail dot com

# Get the lengths of the chromosomes given the header in a .sam file

samtools view -S -H $1 | \
  grep "SN" |
  awk 'BEGIN{OFS="\t"} {print $2, $3}' | \
  sed 's/SN://' | \
  sed 's/LN://' 
