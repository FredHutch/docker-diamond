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


def clean_fastq_headers(fp_in, fp_out):
    """Read in a FASTQ file and write out a copy with unique headers."""

    # Constraints
    # 1. Headers start with '@'
    # 2. Headers are stripped to the first whitespace
    # 3. Headers are unique
    # 4. Sequence lines are not empty
    # 5. Spacer lines match the header line
    # 6. Quality lines are not empty

    with open(fp_in, "rt") as f_in:
        with open(fp_out, "wt") as f_out:
            # Keep track of the line number
            for ix, line in enumerate(f_in):
                # Get the line position 0-3
                mod = ix % 4

                if mod == 0:
                    # Skip lines that are blank (at the end of the file)
                    if len(line) == 1:
                        continue
                    # 1. Headers start with '@'
                    assert line[0] == '@', "Header lacks '@' ({})".format(line)

                    # 2. Strip to the first whitespace
                    line = line.rstrip("\n").split(" ")[0].split("\t")[0]

                    # 3. Add a unique line number and the newline
                    line = "{}-r{}\n".format(line, 1 + (ix / 4))

                    # Save the header to use for the spacer line
                    header = line[1:]

                elif mod == 1:
                    # 4. Sequence lines are not empty
                    assert len(line) > 1

                elif mod == 2:
                    # 5. Spacer lines start with '+' and match the header
                    assert line[0] == "+"
                    line = "+" + header

                elif mod == 3:
                    # 6. Quality lines are not empty
                    assert len(line) > 1

                # Write out the line
                f_out.write(line)
