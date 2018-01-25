#!/usr/bin/python
"""Align a set of reads against a reference database with DIAMOND, and save the results."""

import os
import sys
import time
import json
import uuid
import boto3
import logging
import argparse
import traceback
import subprocess
from helpers.parse_blast import BlastParser
from helpers.fastq_utils import count_fastq_reads
from helpers.fastq_utils import clean_fastq_headers


def run_cmds(commands, retry=0, catchExcept=False):
    """Run commands and write out the log, combining STDOUT & STDERR."""
    logging.info("Commands:")
    logging.info(' '.join(commands))
    p = subprocess.Popen(commands,
                         stdout=subprocess.PIPE,
                         stderr=subprocess.STDOUT)
    stdout, stderr = p.communicate()
    exitcode = p.wait()
    if stdout:
        logging.info("Standard output of subprocess:")
        for line in stdout.split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.split('\n'):
            logging.info(line)

    # Check the exit code
    if exitcode != 0 and retry > 0:
        msg = "Exit code {}, retrying {} more times".format(exitcode, retry)
        logging.info(msg)
        run_cmds(commands, retry=retry - 1)
    elif exitcode != 0 and catchExcept:
        msg = "Exit code was {}, but we will continue anyway"
        logging.info(msg.format(exitcode))
    else:
        assert exitcode == 0, "Exit code {}".format(exitcode)


def calc_abund(input_str,
               db_fp,
               db_url,
               output_folder,
               evalue=0.00001,
               blocks=1,
               query_gencode=11,
               threads=16,
               temp_folder='/mnt/temp',
               random_string=uuid.uuid4(),
               overwrite=False,
               align_mode="blastx"):
    """Align a set of reads against a reference database."""

    # Record the start time
    start_time = time.time()

    # Use the read prefix to name the output and temporary files
    read_prefix = input_str.split('/')[-1]

    # Define the location of temporary file used for DIAMOND output
    blast_fp = os.path.join(
        temp_folder,
        '{}-{}.blast'.format(random_string, read_prefix))

    # Make sure that the temporary file does not already exist
    assert os.path.exists(blast_fp) is False, "Alignment file already exists"

    # Check to see if the output already exists, if so, skip this sample
    output_fp = output_folder.rstrip('/') + '/' + read_prefix + '.json.gz'
    if output_fp.startswith('s3://'):
        # Check S3
        logging.info("Making sure that the output path doesn't already exist on S3")
        bucket = output_fp[5:].split('/')[0]
        prefix = '/'.join(output_fp[5:].split('/')[1:])
        client = boto3.client('s3')
        results = client.list_objects(Bucket=bucket, Prefix=prefix)
        if 'Contents' in results:
            if overwrite:
                logging.info("Overwriting output ({})".format(output_fp))
            else:
                logging.info("Output already exists, skipping ({})".format(output_fp))
                return
    else:
        # Check local filesystem
        if os.path.exists(output_fp):
            if overwrite:
                logging.info("Overwriting output ({})".format(output_fp))
            else:
                logging.info("Output already exists, skipping ({})".format(output_fp))
                return

    # Get the reads
    read_fp = get_reads_from_url(
        input_str,
        temp_folder,
        random_string=random_string
    )

    # Align the reads against the reference database
    logging.info("Aligning reads")
    align_reads(read_fp,
                db_fp,
                blast_fp,
                threads=threads,
                evalue=evalue,
                blocks=blocks,
                query_gencode=query_gencode,
                align_mode=align_mode)

    # Parse the alignment to get the abundance summary statistics
    logging.info("Parsing the output")
    parser = BlastParser(blast_fp, logging=logging)
    parser.parse()
    aligned_reads, abund_summary = parser.make_summary()

    # Count the total number of reads
    logging.info("Counting the total number of reads")
    n_reads = count_fastq_reads(read_fp)
    logging.info("Reads in input file: {}".format(n_reads))

    os.unlink(blast_fp)
    os.unlink(read_fp)

    # Read in the logs
    logging.info("Reading in the logs")
    logs = open(log_fp, 'rt').readlines()

    # Make an object with all of the results
    out = {
        "input_path": input_str,
        "input": read_prefix,
        "output_folder": output_folder,
        "logs": logs,
        "ref_db": db_fp,
        "ref_db_url": db_url,
        "results": abund_summary,
        "aligned_reads": aligned_reads,
        "total_reads": n_reads,
        "time_elapsed": time.time() - start_time
    }

    # Write out the final results as a JSON object and write them to the output folder
    return_results(out, read_prefix, output_folder, temp_folder)

    # Delete any temporary files that might be hanging around
    for fp in os.listdir(temp_folder):
        if fp.startswith(read_prefix):
            if fp.endswith(".json.gz"):
                continue
            fp = os.path.join(temp_folder, fp)
            logging.info("Removing temporary file: {}".format(fp))
            os.unlink(fp)


def get_sra(accession, temp_folder):
    """Get the FASTQ for an SRA accession via ENA."""
    local_path = os.path.join(temp_folder, accession + ".fastq")
    # Download from ENA via FTP
    # See https://www.ebi.ac.uk/ena/browse/read-download for URL format
    url = "ftp://ftp.sra.ebi.ac.uk/vol1/fastq"
    folder1 = accession[:6]
    url = "{}/{}".format(url, folder1)
    if len(accession) > 9:
        if len(accession) == 10:
            folder2 = "00" + accession[-1]
        elif len(accession) == 11:
            folder2 = "0" + accession[-2:]
        elif len(accession) == 12:
            folder2 = accession[-3:]
        else:
            logging.info("This accession is too long: " + accession)
            assert len(accession) <= 12
        url = "{}/{}".format(url, folder2)
    # Add the accession to the URL
    url = "{}/{}/{}".format(url, accession, accession)
    logging.info("Base info for downloading from ENA: " + url)
    # There are three possible file endings
    file_endings = ["_1.fastq.gz", "_2.fastq.gz", ".fastq.gz"]
    # Try to download each file
    for end in file_endings:
        run_cmds(["curl",
                  "-o", os.path.join(temp_folder, accession + end),
                  url + end], catchExcept=True)
    # If none of those URLs downloaded, fall back to trying NCBI
    if any([os.path.exists("{}/{}{}".format(temp_folder, accession, end))
            for end in file_endings]):
        # Combine them all into a single file
        logging.info("Combining into a single FASTQ file")
        with open(local_path, "wt") as fo:
            cmd = "gunzip -c {}/{}*fastq.gz".format(temp_folder, accession)
            gunzip = subprocess.Popen(cmd, shell=True, stdout=fo)
            gunzip.wait()

        # Clean up the temporary files
        logging.info("Cleaning up temporary files")
        for end in file_endings:
            fp = "{}/{}{}".format(temp_folder, accession, end)
            if os.path.exists(fp):
                os.unlink(fp)
    else:
        logging.info("No files found on ENA, trying SRA")
        local_sra_path = os.path.join(temp_folder, accession + ".sra")
        run_cmds([
            "curl", "-o",
            local_sra_path,
            sra_url(accession)])
        run_cmds([
            "fastq-dump", "--split-files", "--outdir", temp_folder, local_sra_path
        ])
        # Combine any multiple files that were found
        with open(local_path + ".temp", "wt") as fo:
            cmd = "cat {}/{}*fastq".format(temp_folder, accession)
            cat = subprocess.Popen(cmd, shell=True, stdout=fo)
            cat.wait()

        if os.path.exists(local_path + ".temp"):
            run_cmds(["mv", local_path + ".temp", local_path])

        # Check to see if the file was downloaded
        msg = "File could not be downloaded from SRA: {}".format(accession)
        assert os.path.exists(local_path), msg

        # Remove the temporary SRA file
        run_cmds(["rm", local_sra_path])

    # Return the path to the file
    logging.info("Done fetching " + accession)
    return local_path


def sra_url(accession):
    """Path to the SRA file accessible via FTP."""
    base = "ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra"
    first_three = accession[:3]
    first_six = accession[:6]
    url = "{}/{}/{}/{}/{}.sra"
    url = url.format(
        base,
        first_three,
        first_six,
        accession,
        accession
        )
    return url


def get_reads_from_url(input_str, temp_folder, random_string=uuid.uuid4()):
    """Get a set of reads from a URL -- return the downloaded filepath."""
    logging.info("Getting reads from {}".format(input_str))

    filename = input_str.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if not input_str.startswith(('s3://', 'sra://', 'ftp://')):
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(input_str)
        assert os.path.exists(input_str), msg
        logging.info("Copying to temporary folder, cleaning up headers")

        # Add a random string to the filename
        local_path = local_path.split('/')
        local_path[-1] = "{}-{}".format(random_string, local_path[-1])
        local_path = '/'.join(local_path)

        # Make the FASTQ headers unique
        clean_fastq_headers(input_str, local_path)

        return local_path

    # Get files from AWS S3
    if input_str.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', input_str, temp_folder
            ])

    # Get files from an FTP server
    elif input_str.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, input_str])

    # Get files from SRA
    elif input_str.startswith('sra://'):
        accession = filename
        logging.info("Getting reads from SRA: " + accession)
        local_path = get_sra(accession, temp_folder)

    else:
        raise Exception("Did not recognize prefix to fetch reads: " + input_str)

    # Add a random string to the filename
    new_path = local_path.split('/')
    new_path[-1] = "{}-{}".format(random_string, new_path[-1])
    new_path = '/'.join(new_path)
    logging.info(
        "Copying {} to {}, cleaning up FASTQ headers".format(
            local_path, new_path
            )
        )
    clean_fastq_headers(local_path, new_path)
    logging.info("Deleting old file: {}".format(local_path))
    os.unlink(local_path)
    return new_path


def get_reference_database(ref_db, temp_folder, random_string=uuid.uuid4()):
    """Get a reference database."""
    assert ref_db.endswith('.dmnd'), "Ref DB must be *.dmnd ({})".format(ref_db)
    # Get files from AWS S3
    if ref_db.startswith('s3://'):
        logging.info("Getting reference database from S3: " + ref_db)

        # Save the database to a local path with a random string prefix, to avoid collision
        local_fp = os.path.join(temp_folder, "{}.{}".format(random_string, ref_db.split('/')[-1]))

        assert os.path.exists(local_fp) is False

        logging.info("Saving database to " + local_fp)
        run_cmds(['aws', 's3', 'cp', '--quiet', '--sse', 'AES256', ref_db, local_fp])

        return local_fp[:-5]

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
                query_gencode=11,
                align_mode="blastx"):
    """Align the reads against the reference database."""
    if align_mode == "blastx":
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
    elif align_mode == "blastp":
        run_cmds(["diamond",
                  "blastp",
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
                  str(blocks)])


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
        run_cmds(['aws', 's3', 'cp', '--quiet', '--sse', 'AES256', temp_fp, output_folder])
        os.unlink(temp_fp)
    else:
        # Copy to local folder
        run_cmds(['mv', temp_fp, output_folder])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Align a set of reads against a reference database with DIAMOND, and save the results.
    """)

    parser.add_argument("--input",
                        type=str,
                        help="""Location for input file(s). Comma-separated.
                                (Supported: sra://, s3://, or ftp://).""")
    parser.add_argument("--ref-db",
                        type=str,
                        help="""Folder containing reference database.
                                (Supported: s3://, ftp://, or local path).""")
    parser.add_argument("--output-folder",
                        type=str,
                        help="""Folder to place results.
                                (Supported: s3://, or local path).""")
    parser.add_argument("--overwrite",
                        action="store_true",
                        help="""Overwrite output files. Off by default.""")
    parser.add_argument("--evalue",
                        type=float,
                        default=0.00001,
                        help="E-value used to filter alignments.")
    parser.add_argument("--blocks",
                        type=int,
                        default=5,
                        help="""Number of blocks used when aligning.
                              Value relates to the amount of memory used.""")
    parser.add_argument("--align-mode",
                        type=str,
                        default="blastx",
                        help="Input is nucleotide (blastx) or protein (blastp).")
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
                        default='/share',
                        help="Folder used for temporary files (and ramdisk, if specified).")

    args = parser.parse_args()

    # Make sure that the align mode is either blastx or blastp
    assert args.align_mode in ["blastx", "blastp"]

    # Set a random string, which will be appended to all temporary files
    random_string = uuid.uuid4()

    # Set up logging
    log_fp = '{}-log.txt'.format(random_string)
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

    # Get the reference database
    db_fp = get_reference_database(
        args.ref_db,
        args.temp_folder,
        random_string=random_string
    )
    logging.info("Reference database: " + db_fp)

    # Align each of the inputs and calculate the overall abundance
    for input_str in args.input.split(','):
        logging.info("Processing input argument: " + input_str)
        # Capture in a try statement
        try:
            calc_abund(input_str,             # ID for single sample to process
                       db_fp,                 # Local path to DB
                       args.ref_db,           # URL of ref DB, used for logging
                       args.output_folder,    # Place to put results
                       evalue=args.evalue,
                       blocks=args.blocks,
                       query_gencode=args.query_gencode,
                       threads=args.threads,
                       temp_folder=args.temp_folder,
                       random_string=random_string,
                       overwrite=args.overwrite,
                       align_mode=args.align_mode)
        except:
            # There was some error
            # Capture the traceback
            logging.info("There was an unexpected failure")
            exc_type, exc_value, exc_traceback = sys.exc_info()
            for line in traceback.format_tb(exc_traceback):
                logging.info(line)

            # Delete any files that were created in this process
            for fp in os.listdir(args.temp_folder):
                if fp.startswith(str(random_string)) or \
                   fp.startswith(input_str.split('/')[-1]):
                    logging.info("Deleting temporary file {}".format(fp))
                    os.unlink(os.path.join(args.temp_folder, fp))

            # Exit
            logging.info("Exit type: {}".format(exc_type))
            logging.info("Exit code: {}".format(exc_value))
            sys.exit(exc_value)

    # Delete any files that were created in this process
    for fp in os.listdir(args.temp_folder):
        if fp.startswith(str(random_string)):
            logging.info("Deleting temporary file {}".format(fp))
            os.unlink(os.path.join(args.temp_folder, fp))

    # Stop logging
    logging.info("Done")
    logging.shutdown()
