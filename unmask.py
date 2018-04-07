#!/usr/bin/python

from pyfasta import Fasta
import numpy as np
import sys

fa = Fasta(sys.argv[1]) # softmasked sequence

for seqid, aseq in fa.iteritems():
    aseq = np.array(str(aseq).upper(), dtype="c")
    print ">%s\n%s" % (seqid, aseq.tostring())
