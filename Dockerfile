FROM ubuntu:20.04
MAINTAINER sminot@fredhutch.org

# Install wget, curl, and Python3
RUN apt update && \
	DEBIAN_FRONTEND="noninteractive" \
	apt-get install -y wget curl build-essential \
	python3-dev python3-pip python3

# Install BioPython
RUN pip3 install biopython==1.70

# Install DIAMOND v2.0.6
RUN mkdir /usr/diamond && cd /usr/diamond && \
	wget https://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz && \
	tar xzvf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz
