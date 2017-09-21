# docker-diamond
Docker image running DIAMOND


[![Docker Repository on Quay](https://quay.io/repository/fhcrc-microbiome/docker-diamond/status "Docker Repository on Quay")](https://quay.io/repository/fhcrc-microbiome/docker-diamond)


This repository provides a Docker image for running DIAMOND that is compatible with automated analysis on AWS Batch. Specifically, this image includes a wrapper script that:


	1. Creates a working directory in memory as a ramdisk

    2. Downloads reference databases

    3. Downloads input data

    4. Aligns reads with DIAMOND

    5. Calculates summary statistics for each reference (coverage and depth)

    6. Saves the outputs to stable file storage


In order to be compatible with AWS Batch, all of these steps are parameterizable and are run with a single command.


### Wrapper script options

#### --input

Specifies the set of FASTQ reads that will be aligned. Supports files from SRA, S3, or FTP. Use the file prefix to specify the source (`s3://`, `sra://`, or `ftp://`). Note that for SRA, just provide the accession (e.g. `sra://SRR123456`).

#### --ref-db

Path to the DIAMOND reference database (file ending in .dmnd). Supports `s3://`, `ftp://`, or a local path.

#### --output-folder

Folder to place the output in, supporting either `s3://` or a local path. Output files will take the form of `<prefix>.json.gz`, where `<prefix>` is the SRA accession (if specified), or otherwise the prefix of the input file from S3 or ftp. 

#### --scratch-size

Size of the temporary folder that will be created to store the reference database, input data, temporary files, and output (in Gb). The folder will be created in memory, and so it will directly reduce the amount of memory that can be used for alignment.

#### --evalue

The evalue used to filter alignments. Defaults to 0.00001.

#### --blocks

Number of 'blocks' used by DIAMOND when loading the reference database in for alignment. According to the DIAMOND manual, the amount of memory used will be roughly 6X the number of blocks (in Gb). So setting `--blocks` to 5 would result in ~30Gb of memory being used during the alignment. 

#### --query-gencode

Genetic code used for six-frame translation, defaults to 11 for microbial.

#### --threads

Number of threads used by DIAMOND during alignment, defaults to 16.

#### --temp-folder

The path to the folder used to create the temporary ramdisk in for storage. There is no obvious reason why a user would need to change this setting.


### Output format

The output of this analysis will summarize the abundance of each individual protein in the reference database. The output is in JSON format, with the following fields:

```
{
  "input_path": <PATH_TO_INPUT_DATA>,
  "input": <INPUT_DATA_PREFIX>,
  "output_folder": <OUTPUT_FOLDER_PATH>,
  "logs": <ANALYSIS_LOGS>,
  "ref_db": <PATH_TO_REF_DB>,
  "results": [
    {
      "id": "gene1",
      "length": 340,
      "total_depth": 1.5,
      "total_coverage": 0.8,
      "total_rpkm": 234.156,
      "unique_depth": 0.75,
      "unique_coverage": 0.3,
      "unique_rpkm": 24.561
    },
    {
      "id": "gene2",
      "length": 120,
      "total_depth": 4.6,
      "total_coverage": 0.98,
      "total_rpkm": 534.156,
      "unique_depth": 1.21,
      "unique_coverage": 0.36,
      "unique_rpkm": 55.751
    },
    ...
  ]
}
```
