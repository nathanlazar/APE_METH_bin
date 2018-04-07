#!/bin/bash

~/APE_METH_bin/LAVA_general.R \
  /u0/dbase/nl/APE_METH/ \
  ~/APE_METH_bin/ \
  LAVA/Human/ \
  human_hg19_noIUPAC/lengths.txt \
  irrelevant \
  RefSeqGenes/Human_hg19_OtherRefSeq.txt \
  LAVA/Gibbon/genes_w_lava.txt > LAVA/Human/permute.out2

~/APE_METH_bin/LAVA_general.R \
  /u0/dbase/nl/APE_METH/ \
  ~/APE_METH_bin/ \
  LAVA/Rhesus/ \
  rhesus_v3/lengths.txt \
  irrelevant \
  RefSeqGenes/Rhesus_rheMac3_OtherRefSeq.txt \
  LAVA/Gibbon/genes_w_lava.txt > LAVA/Rhesus/permute.out2

~/APE_METH_bin/LAVA_general.R \
  /u0/dbase/nl/APE_METH/ \
  ~/APE_METH_bin/ \
  LAVA/Chimp/ \
  chimp_panTro4/lengths.txt \
  irrelevant \
  RefSeqGenes/Chimp_panTro4_OtherRefSeq.txt \
  LAVA/Gibbon/genes_w_lava.txt > LAVA/Chimp/permute.out2

~/APE_METH_bin/LAVA_general.R \
  /u0/dbase/nl/APE_METH/ \
  ~/APE_METH_bin/ \
  LAVA/Orangutan/ \
  orangutan_ponAbe2/lengths.txt \
  irrelevant \
  RefSeqGenes/Orangutan_ponAbe2_OtherRefSeq.txt \
  LAVA/Gibbon/genes_w_lava.txt > LAVA/Orangutan/permute.out2

~/APE_METH_bin/LAVA_general.R \
  /u0/dbase/nl/APE_METH/ \
  ~/APE_METH_bin/ \
  LAVA/Gorilla/ \
  gorilla_gorGor3/lengths.txt \
  irrelevant \
  RefSeqGenes/Gorilla_gorGor3_OtherRefSeq.txt \
  LAVA/Gibbon/genes_w_lava.txt > LAVA/Gorilla/permute.out2
