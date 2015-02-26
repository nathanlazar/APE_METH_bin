#!/bin/bash

# Human
   ./gibbon_meth/batch_align.sh \
     human_hg19_noIUPAC/human_hg19_noIUPAC.fa \
     TRIMMED/Human_1.fq.gz \
     TRIMMED/Human_2.fq.gz \
     HUMAN_MAPPED \
     500000 \
     50 &> Human_align.out

# Gorilla
   ./gibbon_meth/batch_align.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     TRIMMED/Gorilla_1.fq.gz \
     TRIMMED/Gorilla_2.fq.gz \
     GORILLA_MAPPED \
     500000 \
     50 &> Gorilla_align.out

# Chimp
   ./gibbon_meth/batch_align.sh \
     chimp_panTro4/chimp_panTro4.fa \
     TRIMMED/Chimp_1.fq.gz \
     TRIMMED/Chimp_2.fq.gz \
     CHIMP_MAPPED \
     500000 \
     50 &> Chimp_align.out

# Orangutan
   ./gibbon_meth/batch_align.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     TRIMMED/Orang_1.fq.gz \
     TRIMMED/Orang_2.fq.gz \
     MAPPED \
     ORANG_500000 \
     50 &> Orang_align.out

# Rhesus
   ./gibbon_meth/batch_align.sh \
     rhesus_v7/rhesus_v7.fa \
     TRIMMED/Rhesus_1.fq.gz \
     TRIMMED/Rhesus_2.fq.gz \
     RHESUS_MAPPED \
     500000 \
     50 &> Rhesus_align.out


