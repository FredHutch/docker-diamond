#!/usr/bin/env python3
"""Run DIAMOND and return the BLAST alignment file for two FASTA files."""

import os
import sys
import time
import json
import uuid
import shutil
import logging
import argparse
import traceback
import subprocess


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
        for line in stdout.decode("utf-8").split('\n'):
            logging.info(line)
    if stderr:
        logging.info("Standard error of subprocess:")
        for line in stderr.decode("utf-8").split('\n'):
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


def get_file_from_url(url_path, temp_folder):
    """Get a file from a URL -- return the downloaded filepath."""
    logging.info("Getting file from {}".format(url_path))

    filename = url_path.split('/')[-1]
    local_path = os.path.join(temp_folder, filename)

    if temp_folder.endswith("/") is False:
        temp_folder = temp_folder + "/"

    logging.info("Filename: " + filename)
    logging.info("Local path: " + local_path)

    if "://" in url_path is False:
        logging.info("Treating as local path")
        msg = "Input file does not exist ({})".format(url_path)
        assert os.path.exists(url_path), msg
        logging.info("Making symbolic link in temporary folder")
        os.symlink(url_path, local_path)

    # Get files from AWS S3
    elif url_path.startswith('s3://'):
        logging.info("Getting reads from S3")
        run_cmds([
            'aws', 's3', 'cp', '--quiet', '--sse',
            'AES256', url_path, temp_folder
        ])

    # Get files from an FTP server
    elif url_path.startswith('ftp://'):
        logging.info("Getting reads from FTP")
        run_cmds(['wget', '-P', temp_folder, url_path])

    else:
        raise Exception(
            "Did not recognize prefix to fetch reads: " + url_path)

    return local_path


def exit_and_clean_up(temp_folder):
    """Log the error messages and delete the temporary folder."""
    # Capture the traceback
    logging.info("There was an unexpected failure")
    exc_type, exc_value, exc_traceback = sys.exc_info()
    for line in traceback.format_tb(exc_traceback):
        logging.info(line)

    # Delete any files that were created for this sample
    logging.info("Removing temporary folder: " + temp_folder)
    shutil.rmtree(temp_folder)

    # Exit
    logging.info("Exit type: {}".format(exc_type))
    logging.info("Exit code: {}".format(exc_value))
    sys.exit(exc_value)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""
    Run BLAST (DIAMOND) on a pair of FASTA files.
    """)

    parser.add_argument("--query",
                        type=str,
                        required=True,
                        help="""Location for 'query' input file.
                                (Supported: local path, s3://, or ftp://).""")
    parser.add_argument("--subject",
                        type=str,
                        required=True,
                        help="""Location for 'subject' input file.
                                (Supported: local path, s3://, or ftp://).""")
    parser.add_argument("--output-aln",
                        type=str,
                        required=True,
                        help="""URL for output alignment file (local path, or S3://).""")
    parser.add_argument("--outfmt",
                        type=str,
                        default="6",
                        help="""Output format for DIAMOND -- default of 6 is tabular BLAST-format.""")
    parser.add_argument("--perc-identity",
                        type=float,
                        default=50.,
                        help="""Minimum percent identity for alignments.""")
    parser.add_argument("--query-cover",
                        type=float,
                        default=50.,
                        help="""Percent query coverage.""")
    parser.add_argument("--subject-cover",
                        type=float,
                        default=50.,
                        help="""Percent subject coverage.""")
    parser.add_argument("--max-target-seqs",
                        type=int,
                        default=500,
                        help="""Maximum number of alignments to report per query.""")
    parser.add_argument("--output-log",
                        type=str,
                        required=True,
                        help="""URL for output log file (local path, or S3://).""")
    parser.add_argument("--blast-type",
                        type=str,
                        default="blastp",
                        help="""Type of BLAST to run (blastx or blastp).""")
    parser.add_argument("--query-gencode",
                        type=int,
                        default=11,
                        help="""Genetic code to use in case of blastx.""")
    parser.add_argument("--threads",
                        type=int,
                        default=16,
                        help="Number of threads to use assembling.")
    parser.add_argument("--temp-folder",
                        type=str,
                        default='/share',
                        help="Folder used for temporary files.")

    args = parser.parse_args()

    # Make sure the type of BLAST is allowed
    msg = "Blast type '{}' not recognized".format(args.blast_type)
    assert args.blast_type in ["blastp", "blastx"], msg

    # Check that the temporary folder exists
    assert os.path.exists(args.temp_folder)

    # Set a random string, which will be appended to all temporary files
    random_string = str(uuid.uuid4())[:8]

    # Make a temporary folder within the --temp-folder with the random string
    temp_folder = os.path.join(args.temp_folder, str(random_string))
    # Make sure it doesn't already exist
    msg = "Collision, {} already exists".format(temp_folder)
    assert os.path.exists(temp_folder) is False, msg
    # Make the directory
    os.mkdir(temp_folder)
    
    # Make folders for the query, subject, and output
    temp_folders = {}
    for n in ["query", "subject", "output", "db"]:
        temp_folders[n] = os.path.join(temp_folder, n)
        os.mkdir(temp_folders[n])

    # Set up logging
    log_fp = '{}/log.txt'.format(temp_folder)
    logFormatter = logging.Formatter(
        '%(asctime)s %(levelname)-8s [run_blast.py] %(message)s')
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

    # Get the subject and query files
    try:
        query_fp = get_file_from_url(args.query, temp_folders["query"])
    except:
        exit_and_clean_up(temp_folder)
    try:
        subject_fp = get_file_from_url(args.subject, temp_folders["subject"])
    except:
        exit_and_clean_up(temp_folder)

    logging.info("Query: " + query_fp)
    logging.info("Subject: " + subject_fp)

    # Make a reference database from the subject
    logging.info("Making reference database")
    db_prefix = os.path.join(temp_folders["db"], "db")
    try:
        run_cmds([
            "diamond",
            "makedb",
            "--in", subject_fp,
            "-d", db_prefix
        ])
    except:
        exit_and_clean_up(temp_folder)
    logging.info("Done making reference database")

    # Set up the arguments
    output_fp = os.path.join(temp_folders["output"], "output.aln")
    arg_list = [
        "diamond",
        args.blast_type,
        "--db", db_prefix,
        "--query", query_fp,
        "--out", output_fp,
        "--outfmt", args.outfmt,
        "--id", str(args.perc_identity),
        "--max_target_seqs", str(args.max_target_seqs),
        "--query-cover", str(args.query_cover),
        "--subject-cover", str(args.subject_cover),
        "--threads", str(args.threads)
    ]
    if args.blast_type == "blastx":
        # For BLASTX, specify the genetic code
        arg_list = arg_list + [
            "--query-gencode", str(args.query_gencode),
        ]
    # Gzip if specified
    if args.output_aln.endswith(".gz"):
        arg_list = arg_list + [
            "--compress", "1"
        ]

    # Run the command
    try:
            run_cmds(arg_list)
    except:
        exit_and_clean_up(temp_folder)

    # Return the results
    for local_path, remote_path in [
        (output_fp, args.output_aln),
        (log_fp, args.output_log)
    ]:
        if remote_path.startswith("s3://"):
            try:
                run_cmds([
                    "aws",
                    "s3",
                    "cp",
                    local_path,
                    remote_path
                ])
            except:
                exit_and_clean_up(temp_folder)
        else:
            try:
                run_cmds([
                    "mv",
                    local_path,
                    remote_path
                ])
            except:
                exit_and_clean_up(temp_folder)

    # Delete any files that were created in this process
    logging.info("Deleting temporary folder: {}".format(temp_folder))
    shutil.rmtree(temp_folder)

    # Stop logging
    logging.info("Done")
    logging.shutdown()
