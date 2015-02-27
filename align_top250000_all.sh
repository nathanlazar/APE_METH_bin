#!/bin/bash

bin_dir=/mnt/lustre1/users/lazar/APE_METH/POST_CRASH/APE_METH_bin

# Gibbon
   $bin_dir/align_top250000.sh \
     NomLeu1_0/NomLeu1_0.fa \
     ADAPT_TRIM/DNA111101LC_62_HSA_normal_NoIndex_L006_1_val_1.fq.gz \
     ADAPT_TRIM/DNA111101LC_62_HSA_normal_NoIndex_L006_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/DNA111101LC_62_HSA_normal_NoIndex_L006.out

   $bin_dir/align_top250000.sh \
     NomLeu1_0/NomLeu1_0.fa \
     ADAPT_TRIM/s_3_1_val_1.fq.gz \
     ADAPT_TRIM/s_3_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/s_3.out

   $bin_dir/align_top250000.sh \
     NomLeu1_0/NomLeu1_0.fa \
     ADAPT_TRIM/s_8_1_val_1.fq.gz \
     ADAPT_TRIM/s_8_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/s_8.out

# Human
   $bin_dir/align_top250000.sh 
     human_hg19_noIUPAC/human_hg19_noIUPAC.fa 
     ADAPT_TRIM/Julia_FCC02FUACXX_L6_1_val_1.fq.gz 
     ADAPT_TRIM/Julia_FCC02FUACXX_L6_2_val_2.fq.gz 
     TEST_MAP &> TEST_MAP/Julia_FCC02FUACXX_L6.out

   $bin_dir/align_top250000.sh
     human_hg19_noIUPAC/human_hg19_noIUPAC.fa
     ADAPT_TRIM/Julia_FCC02FUACXX_L8_1_val_1.fq.gz
     ADAPT_TRIM/Julia_FCC02FUACXX_L8_2_val_2.fq.gz
     TEST_MAP &> TEST_MAP/Julia_FCC02FUACXX_L8.out

# Gorilla
   $bin_dir/align_top250000.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     ADAPT_TRIM/Gorilla_FCD0GB2ACXX_L6_1_val_1.fq.gz \
     ADAPT_TRIM/Gorilla_FCD0GB2ACXX_L6_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Gorilla_FCD0GB2ACXX_L6.out

   $bin_dir/align_top250000.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L7_1_val_1.fq.gz \
     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L7_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Gorilla_FCC02G0ACXX_L7.out

   $bin_dir/align_top250000.sh \
     gorilla_gorGor3/gorilla_gorGor3.fa \
     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L8_1_val_1.fq.gz \
     ADAPT_TRIM/Gorilla_FCC02G0ACXX_L8_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Gorilla_FCC02G0ACXX_L8.out

# Chimp
   $bin_dir/align_top250000.sh \
     chimp_panTro4/chimp_panTro4.fa \
     ADAPT_TRIM/4369_FCC02G0ACXX_L2_1_val_1.fq.gz \
     ADAPT_TRIM/4369_FCC02G0ACXX_L2_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/4369_FCC02G0ACXX_L2.out

   $bin_dir/align_top250000.sh \
     chimp_panTro4/chimp_panTro4.fa \
     ADAPT_TRIM/4369_FCC02G0ACXX_L3_1_val_1.fq.gz
     ADAPT_TRIM/4369_FCC02G0ACXX_L3_2_val_2.fq.gz
     TEST_MAP &> TEST_MAP/4369_FCC02G0ACXX_L3.out

# Orangutan
   $bin_dir/align_top250000.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     ADAPT_TRIM/Dunja_FCC0CGEACXX_L4_1_val_1.fq.gz \
     ADAPT_TRIM/Dunja_FCC0CGEACXX_L4_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Dunja_FCC0CGEACXX_L4.out

   $bin_dir/align_top250000.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     ADAPT_TRIM/Dunja_FCC0CGUACXX_L4_1_val_1.fq.gz \
     ADAPT_TRIM/Dunja_FCC0CGUACXX_L4_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Dunja_FCC0CGUACXX_L4.out

   $bin_dir/align_top250000.sh \
     orangutan_ponAbe2/orangutan_ponAbe2.fa \
     ADAPT_TRIM/Dunja_FCD0JPEACXX_L1_1_val_1.fq.gz \
     ADAPT_TRIM/Dunja_FCD0JPEACXX_L1_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Dunja_FCD0JPEACXX_L1.out

# Rhesus
   $bin_dir/align_top250000.sh \
     rhesus_v7/rhesus_v7.fa \
     ADAPT_TRIM/Rhesus_DNA130327BF_25184-0_1_val_1.fq.gz \
     ADAPT_TRIM/Rhesus_DNA130327BF_25184-0_2_val_2.fq.gz \
     TEST_MAP &> TEST_MAP/Rhesus_DNA130327BF_25184-0.out
