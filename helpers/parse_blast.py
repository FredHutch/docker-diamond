#!/usr/bin/python
"""Logic needed to parse a set of BLAST results."""

import json
import argparse
from collections import defaultdict


class BlastParser:
    """Object to parse a set of BLAST results."""

    def __init__(self, blast_fp,
                 qid_ix=0, sid_ix=1, slen_ix=2, sstart_ix=3, send_ix=4, qseq_ix=5,
                 comment_char='@', sep='\t', logging=False):
        """Parse a set of BLAST results."""
        # Store the input variables
        self.blast_fp = blast_fp
        self.qid_ix = qid_ix
        self.sid_ix = sid_ix
        self.slen_ix = slen_ix
        self.sstart_ix = sstart_ix
        self.send_ix = send_ix
        self.qseq_ix = qseq_ix
        self.comment_char = comment_char
        self.sep = sep
        self.logging = logging

        # Store the results of the previous alignment
        self.last_qid = None
        self.last_sid = None
        self.last_qseq = None
        self.last_alen = None
        self.last_pos_range = set([])

        # Store the most recent repeated sequence
        self.last_repeated_qseq = ''

        # Number of reads and bases aligned to each reference
        self.total_reads = defaultdict(int)
        self.unique_reads = defaultdict(int)
        self.total_bases = defaultdict(int)
        self.unique_bases = defaultdict(int)

        # The unique set of positions that are aligned by each reference
        self.total_pos = defaultdict(set)
        self.unique_pos = defaultdict(set)

        # Keep track of the total number of aligned reads
        self.total_aligned_reads = 0

        # The total length of each reference
        self.ref_len = {}

        # Set the highest number of fields that we might have to parse from a line
        self.max_fields = max(sid_ix, slen_ix, sstart_ix, send_ix, qseq_ix) + 1

    def parse(self):
        """Parse the file."""

        with open(self.blast_fp) as f:
            ix = 0
            for line in f:
                if ix % 100000 == 0:
                    if self.logging:
                        self.logging.info("Processed {} alignments".format(ix))
                # Skip lines starting with '@', by default
                if line[0] == self.comment_char:
                    continue

                qid, sid, qseq, sstart, send, slen = self.parse_line(line)

                # If this is a new query, increment the counter
                if qid != self.last_qid:
                    self.total_aligned_reads += 1
                    self.last_qid = qid

                # Calculate the alignment length
                alen = send - sstart

                # If this is the first time we've seen this subject, add the length
                if sid not in self.ref_len:
                    self.ref_len[sid] = float(slen)

                # For the first alignment, just set the "last_*" variables and move on
                if ix > 0:
                    # Check to see if this alignment is repeated the same query sequence
                    if qseq == self.last_qseq:
                        self.last_repeated_qseq = qseq
                    # If the last alignment was unique, add it as such
                    if self.last_qseq != self.last_repeated_qseq:
                        # The previous alignment was unique
                        self.unique_bases[self.last_sid] += self.last_alen
                        self.unique_reads[self.last_sid] += 1
                        # Add the alignment positions
                        self.unique_pos[self.last_sid] |= self.last_pos_range
                    # No matter what, add the total reads and bases for the previous sequence
                    self.total_bases[self.last_sid] += self.last_alen
                    self.total_reads[self.last_sid] += 1
                    self.total_pos[self.last_sid] |= self.last_pos_range
                self.last_qseq = qseq
                self.last_sid = sid
                self.last_alen = alen
                self.last_pos_range = set(range(sstart, send))

                ix += 1

        # Add the final line
        if self.last_qseq != self.last_repeated_qseq:
            self.unique_bases[self.last_sid] += self.last_alen
            self.unique_reads[self.last_sid] += 1
            self.unique_pos[self.last_sid] |= self.last_pos_range
        self.total_bases[self.last_sid] += self.last_alen
        self.total_reads[self.last_sid] += 1
        self.total_pos[self.last_sid] |= self.last_pos_range

        msg = "Processed {} alignments, {} total aligned reads".format(ix, self.total_aligned_reads)
        if self.logging:
            self.logging.info(msg)

    def parse_line(self, line):
        """Parse one line."""
        line = line.strip('\n').split('\t', self.max_fields)

        # Subject ID
        sid = line[self.sid_ix]
        if sid == '*':
            # Read is not aligned
            return None, None, None, None, None, None

        # Query ID
        qid = line[self.qid_ix]

        # Aligned query sequence
        qseq = line[self.qseq_ix]
        if qseq == '*':
            # Read is not aligned
            return None, None, None, None, None, None

        # Alignment positions
        sstart = int(line[self.sstart_ix])  # Position of alignment start on subject
        send = int(line[self.send_ix])  # Position of alignment end on subject
        # Orient the alignment positions so that sstart < send
        if send < sstart:
            send, sstart = sstart, send

        # Convert to 0-index
        sstart = sstart - 1

        # Total subject (reference) length
        slen = int(line[self.slen_ix])

        return qid, sid, qseq, sstart, send, slen

    def make_summary(self):
        """Make the final output."""

        # Make the output object as a list
        out = []

        for k, v in self.total_bases.items():
            # Information for a single item
            d = {'id': k}
            # Reference length (stored as a float)
            rl = self.ref_len[k]
            d['length'] = int(rl)
            # Depth = total aligned bases / reference length
            d['total_depth'] = round(v / rl, 4)
            d['unique_depth'] = round(self.unique_bases.get(k, 0) / rl, 4)
            # Coverage = number of positions covered / reference length
            d['total_coverage'] = round(len(self.total_pos[k]) / rl, 4)
            d['unique_coverage'] = round(len(self.unique_pos[k]) / rl, 4)
            # RPKM = aligned reads / kilobase of reference / million aligned reads
            d['total_rpkm'] = round(self.rpkm(self.total_reads[k], rl, self.total_aligned_reads), 6)
            d['unique_rpkm'] = round(self.rpkm(self.unique_reads[k], rl, self.total_aligned_reads), 6)

            out.append(d)

        return out

    def rpkm(self, reads, ref_len, total_reads, amino_acid_ref=True):
        """Calculate RPKM (reads per kilobase per million reads)."""
        # For amino acid alignments, coordinates indicate 1aa == 3nt
        if amino_acid_ref:
            return reads / ((3 * ref_len / 1000.) * (total_reads / 1000000.))
        else:
            return reads / ((ref_len / 1000.) * (total_reads / 1000000.))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Parse a set of BLAST results.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location of BLAST results file.""")
    parser.add_argument("--out",
                        type=str,
                        help="""Path to write results out to.""")
    parser.add_argument("--qseqid", default=0, type=int,
                        help=("Index position (0-based) for column with query ID"))
    parser.add_argument("--sseqid", default=1, type=int,
                        help=("Index position (0-based) for column with subject ID"))
    parser.add_argument("--slen", default=2, type=int,
                        help=("Index position (0-based) for column with subject length"))
    parser.add_argument("--sstart", default=3, type=int,
                        help=("Index position (0-based) for column with subject start position"))
    parser.add_argument("--send", default=4, type=int,
                        help=("Index position (0-based) for column with subject end position"))
    parser.add_argument("--qseq", default=5, type=int,
                        help=("Index position (0-based) for column with aligned query sequence"))

    args = parser.parse_args()

    blast_parser = BlastParser(args.input, args.qseqid, args.sseqid,
                               args.slen, args.sstart, args.send,
                               args.qseq)
    blast_parser.parse()
    results = blast_parser.make_summary()
    with open(args.out, 'wt') as fo:
        json.dump(results, fo)
