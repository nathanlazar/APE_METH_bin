#!/bin/bash

# nathan dot lazar at gmail dot com

bin_dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin

# Gibbon
   $bin_dir/batch_align.sh \
     NomLeu1_0/NomLeu1_0.fa \
     TRIMMMED/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1_trim.fq.gz \
     TRIMMMED/DNA111101LC_62_HSA_normal_NoIndex_L006_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/DNA111101LC_62_HSA_normal_NoIndex_L006.out

   $bin_dir/batch_align.sh \
     NomLeu1_0/NomLeu1_0.fa \
     TRIMMMED/s_3_1_val_1_trim.fq.gz \
     TRIMMMED/s_3_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/s_3.out

   $bin_dir/batch_align.sh \
     NomLeu1_0/NomLeu1_0.fa \
     TRIMMMED/s_8_1_val_1_trim.fq.gz \
     TRIMMMED/s_8_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/s_8.out

# Human
   $bin_dir/batch_align.sh \
     human_hg19_noIUPAC/human_hg19_noIUPAC.fa \
     TRIMMMED/Julia_FCC02FUACXX_L6_1_val_1_trim.fq.gz \
     TRIMMMED/Julia_FCC02FUACXX_L6_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Julia_FCC02FUACXX_L6.out

   $bin_dir/batch_align.sh \
     human_hg19_noIUPAC/human_hg19_noIUPAC.fa \
     TRIMMMED/Julia_FCC02FUACXX_L8_1_val_1_trim.fq.gz \
     TRIMMMED/Julia_FCC02FUACXX_L8_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Julia_FCC02FUACXX_L8.out

# Gorilla
   $bin_dir/batch_align.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     TRIMMMED/Gorilla_FCD0GB2ACXX_L6_1_val_1_trim.fq.gz \
     TRIMMMED/Gorilla_FCD0GB2ACXX_L6_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Gorilla_FCD0GB2ACXX_L6.out

   $bin_dir/batch_align.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     TRIMMMED/Gorilla_FCC02G0ACXX_L7_1_val_1_trim.fq.gz \
     TRIMMMED/Gorilla_FCC02G0ACXX_L7_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Gorilla_FCC02G0ACXX_L7.out

   $bin_dir/batch_align.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     TRIMMMED/Gorilla_FCC02G0ACXX_L8_1_val_1_trim.fq.gz \
     TRIMMMED/Gorilla_FCC02G0ACXX_L8_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Gorilla_FCC02G0ACXX_L8.out

# Chimp
   $bin_dir/batch_align.sh \
     chimp_panTro4/chimp_panTro4.fa \
     TRIMMMED/4369_FCC02G0ACXX_L2_1_val_1_trim.fq.gz \
     TRIMMMED/4369_FCC02G0ACXX_L2_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/4369_FCC02G0ACXX_L2.out

   $bin_dir/batch_align.sh \
     chimp_panTro4/chimp_panTro4.fa \
     TRIMMMED/4369_FCC02G0ACXX_L3_1_val_1_trim.fq.gz \
     TRIMMMED/4369_FCC02G0ACXX_L3_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/4369_FCC02G0ACXX_L3.out

# Orangutan
   $bin_dir/batch_align.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     TRIMMMED/Dunja_FCC0CGEACXX_L4_1_val_1_trim.fq.gz \
     TRIMMMED/Dunja_FCC0CGEACXX_L4_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Dunja_FCC0CGEACXX_L4.out

   $bin_dir/batch_align.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     TRIMMMED/Dunja_FCC0CGUACXX_L4_1_val_1_trim.fq.gz \
     TRIMMMED/Dunja_FCC0CGUACXX_L4_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Dunja_FCC0CGUACXX_L4.out

   $bin_dir/batch_align.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     TRIMMMED/Dunja_FCD0JPEACXX_L1_1_val_1_trim.fq.gz \
     TRIMMMED/Dunja_FCD0JPEACXX_L1_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &> MAPPED/Dunja_FCD0JPEACXX_L1.out

# Rhesus
   $bin_dir/batch_align.sh \
     rhesus_v7/rhesus_v7.fa \
     TRIMMMED/Rhesus_DNA130327BF_25184-0_1_val_1_trim.fq.gz \
     TRIMMMED/Rhesus_DNA130327BF_25184-0_2_val_2_trim.fq.gz \
     MAPPED 500000 50 & &>  MAPPED/Rhesus_DNA130327BF_25184-0.out

