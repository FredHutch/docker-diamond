#!/usr/bin/python

import gzip
import logging
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqIO.FastaIO import SimpleFastaParser


def count_fasta_reads(fp):
    n = 0
    if fp.endswith(".gz"):
        with gzip.open(fp, "rt") as f:
            for record in SimpleFastaParser(f):
                n += 1
    else:
        with open(fp, "rt") as f:
            for record in SimpleFastaParser(f):
                n += 1

    return n


def count_fastq_reads(fp):
    n = 0
    if fp.endswith(".gz"):
        with gzip.open(fp, "rt") as f:
            for record in FastqGeneralIterator(f):
                n += 1
    else:
        with open(fp, "rt") as f:
            for record in FastqGeneralIterator(f):
                n += 1

    # If no reads were found, try counting it as a FASTA
    if n == 0:
        logging.info("No FASTQ reads found, trying to read as FASTA")
        n = count_fasta_reads(fp)

    return n
