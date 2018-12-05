FROM ubuntu:16.04
MAINTAINER sminot@fredhutch.org

# Install prerequisites
RUN apt update && \
	apt-get install -y build-essential wget unzip python2.7 python-dev git python-pip bats awscli curl

# Use /share as the working directory
RUN mkdir /share
WORKDIR /share

# Add files
RUN mkdir /usr/diamond
ADD requirements.txt /usr/diamond

# Install python requirements
RUN pip install -r /usr/diamond/requirements.txt && rm /usr/diamond/requirements.txt


# Install DIAMOND v0.9.23
RUN cd /usr/diamond && \
	wget https://github.com/bbuchfink/diamond/releases/download/v0.9.23/diamond-linux64.tar.gz && \
	tar xzvf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz


# Install the SRA toolkit
RUN cd /usr/local/bin && \
	wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.2/sratoolkit.2.8.2-ubuntu64.tar.gz && \
	tar xzvf sratoolkit.2.8.2-ubuntu64.tar.gz && \
	ln -s /usr/local/bin/sratoolkit.2.8.2-ubuntu64/bin/* /usr/local/bin/ && \
	rm sratoolkit.2.8.2-ubuntu64.tar.gz


# Add the run script to the PATH
ADD run.py /usr/diamond
ADD helpers /usr/diamond/helpers
RUN cd /usr/diamond && \
	ln -s /usr/diamond/run.py /usr/bin/
ADD run_blast.py /usr/bin


# Run tests and then remove the folder
ADD tests /usr/diamond/tests
RUN bats /usr/diamond/tests/ && rm -r /usr/diamond/tests/
