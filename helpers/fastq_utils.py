#!/usr/bin/python

from Bio.SeqIO.QualityIO import FastqGeneralIterator


def count_fastq_reads(fp):
    n = 0
    with open(fp, "rt") as f:
        for record in FastqGeneralIterator(f):
            n += 1
    return n
