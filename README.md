# docker-diamond
Docker image running DIAMOND


This repository provides a Docker image for running DIAMOND that is compatible with automated analysis on AWS Batch. Specifically, this image includes a wrapper script that:


    1. Downloads reference databases

    2. Downloads input data

    3. Aligns reads with DIAMOND

    4. Calculates summary statistics for each reference (coverage and depth)

    5. Saves the outputs to stable file storage


In order to be compatible with AWS Batch, all of these steps are parameterizable and are run with a single command.

