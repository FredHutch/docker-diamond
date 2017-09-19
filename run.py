#!/usr/bin/python
"""Align a set of reads against a reference database with DIAMOND, and save the results."""

import os
import json
import logging
import argparse
import subprocess
from helpers.parse_blast import BlastParser


def run_cmds(commands):
    """Run a set of commands and write out the log, combining the STDOUT and STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    p = subprocess.Popen(commands,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        logging.info(stdout)
    if stderr:
        logging.info("Standard error of subprocess:")
        logging.info(stderr)

    # Check the exit code
    assert exitcode == 0, "Exit code {}".format(exitcode)


def calc_abund(input_str,
               ref_db,
               output_folder,
               evalue=0.00001,
               blocks=1,
               query_gencode=11,
               threads=16,
               temp_folder='/mnt/temp'):
    """Align a set of reads against a reference database."""

    # Define the location of temporary file used for DIAMOND output
    blast_fp = os.path.join(temp_folder, 'temp.blast')

    # Make sure that the temporary file does not already exist
    if os.path.exists(blast_fp):
        os.unlink(blast_fp)

    # Get the reads
    read_fp, read_prefix = get_reads_from_url(input_str, temp_folder)

    # Get the reference database
    db_fp = get_reference_database(ref_db, temp_folder)
    logging.info("Reference database: " + db_fp)

    # Align the reads against the reference database
    logging.info("Aligning reads")
    align_reads(read_fp,
                db_fp,
                blast_fp,
                threads=threads,
                evalue=evalue,
                blocks=blocks,
                query_gencode=query_gencode)

    # Parse the alignment to get the abundance summary statistics
    logging.info("Parsing the output")
    parser = BlastParser(blast_fp)
    parser.parse()
    abund_summary = parser.make_summary()

    os.unlink(blast_fp)

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()

    # Make an object with all of the results
    out = {
        "input_path": input_str,
        "input": read_prefix,
        "output_folder": output_folder,
        "logs": logs,
        "ref_db": ref_db,
        "results": abund_summary
    }

    # Write out the final results as a JSON object and write them to the output folder
    return_results(out, read_prefix, output_folder, temp_folder)


def get_reads_from_url(input_str, temp_folder):
    """Get a set of reads from a URL -- return the downloaded filepath and file prefix."""
    logging.info("Getting reads from {}".format(input_str))
    error_msg = "{} must start with s3://, sra://, or ftp://".format(input_str)
    assert input_str.startswith(('s3://', 'sra://', 'ftp://')), error_msg

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    # Get files from AWS S3
    if input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds(['aws', 's3', 'cp', '--sse', 'AES256', input_str, temp_folder])
        return local_path, filename

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])
        return local_path, filename

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = input_str[6:]
        logging.info("Getting reads from SRA: " + accession)
        local_path = os.path.join(temp_folder, accession + ".fastq")
        # Download from NCBI
        run_cmds(["fastq-dump",
                  "--skip-technical",
                  "--readids",
                  "--read-filter",
                  "pass",
                  "--dumpbase",
                  "--clip",
                  "--outdir",
                  temp_folder,
                  accession])

        # Rename the file (which automatically has '_pass' included)
        run_cmds(["mv",
                  os.path.join(temp_folder, accession + "_pass.fastq"),
                  local_path])
        return local_path, accession

    else:
        raise Exception("Did not recognize prefix to fetch reads: " + input_str)


def get_reference_database(ref_db, temp_folder):
    """Get a reference database."""
    assert ref_db.endswith('.dmnd'), "Ref DB must be *.dmnd ({})".format(ref_db)
    # Get files from AWS S3
    if ref_db.startswith('s3://'):
        logging.info("Getting reference database from S3: " + ref_db)
        run_cmds(['aws', 's3', 'cp', '--sse', 'AES256', ref_db, temp_folder])

        # Return the prefix for the reference database
        database_prefix = ref_db.split('/')[-1][:-5]
        return os.path.join(temp_folder, database_prefix)

    else:
        # Treat the input as a local path
        logging.info("Getting reference database from local path: " + ref_db)
        assert os.path.exists(ref_db)
        return ref_db


def align_reads(read_fp,
                db_fp,
                blast_fp,
                threads=16,
                evalue=0.00001,
                blocks=1,
                query_gencode=11):
    """Align the reads against the reference database."""
    run_cmds(["diamond",
              "blastx",
              "--threads",
              str(threads),
              "--query",
              read_fp,
              "--db",
              db_fp,
              "--outfmt",
              "6",
              "qseqid",
              "sseqid",
              "slen",
              "sstart",
              "send",
              "qseq",
              "--out",
              blast_fp,
              "--top",
              "0",
              "--evalue",
              str(evalue),
              "-b",
              str(blocks),
              "--query-gencode",
              str(query_gencode)])


def return_results(out, read_prefix, output_folder, temp_folder):
    """Write out the final results as a JSON object and write them to the output folder."""
    # Make a temporary file
    temp_fp = os.path.join(temp_folder, read_prefix + '.json')
    with open(temp_fp, 'wt') as fo:
        json.dump(out, fo)
    # Compress the output
    run_cmds(['gzip', temp_fp])
    temp_fp = temp_fp + '.gz'

    if output_folder.startswith('s3://'):
        # Copy to S3
        run_cmds(['aws', 's3', 'cp', '--sse', 'AES256', temp_fp, output_folder])
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_folder])


def make_scratch_space(scratch_size, temp_folder):
    """Create scratch space using a ramdisk."""
    run_cmds(['mount', '-t', 'tmpfs', '-o', 'size={}g'.format(scratch_size),
              'tmpfs', temp_folder])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Align a set of reads against a reference database with DIAMOND, and save the results.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location for input file(s). Comma-separated.
                                (Supported: sra://, ftp://).""")
    parser.add_argument("--ref-db",
                        type=str,
                        help="""Folder containing reference database.
                                (Supported: s3://, ftp://).""")
    parser.add_argument("--output-folder",
                        type=str,
                        help="""Folder to place results.
                                (Supported: s3://, ftp://).""")
    parser.add_argument("--scratch-size",
                        type=int,
                        default=20,
                        help="Size of scratch space created with ramdisk (Gb).")
    parser.add_argument("--evalue",
                        type=float,
                        default=0.00001,
                        help="E-value used to filter alignments.")
    parser.add_argument("--blocks",
                        type=int,
                        default=5,
                        help="""Number of blocks used when aligning.
                              Value relates to the amount of memory used.""")
    parser.add_argument("--query-gencode",
                        type=int,
                        default=11,
                        help="Genetic code used to translate nucleotide reads.")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use aligning.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/mnt/temp',
                        help="Folder used to mount ramdisk used for temporary files.")

    args = parser.parse_args()

    # Set up logging
    log_fp = 'log.txt'
    logFormatter = logging.Formatter('%(asctime)s %(levelname)-8s [run.py] %(message)s')
    rootLogger = logging.getLogger()
    rootLogger.setLevel(logging.INFO)

    # Write to file
    fileHandler = logging.FileHandler(log_fp)
    fileHandler.setFormatter(logFormatter)
    rootLogger.addHandler(fileHandler)
    # Also write to STDOUT
    consoleHandler = logging.StreamHandler()
    consoleHandler.setFormatter(logFormatter)
    rootLogger.addHandler(consoleHandler)

    # Set up the scratch space
    logging.info("Setting up scratch space ({}Gb)".format(args.scratch_size))
    make_scratch_space(args.scratch_size, args.temp_folder)

    # Align each of the inputs and calculate the overall abundance
    for input_str in args.input.split(','):
        logging.info("Processing input argument: " + input_str)
        calc_abund(input_str,
                   args.ref_db,
                   args.output_folder,
                   evalue=args.evalue,
                   blocks=args.blocks,
                   query_gencode=args.query_gencode,
                   threads=args.threads,
                   temp_folder=args.temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()
