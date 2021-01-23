FROM ubuntu:20.04
MAINTAINER sminot@fredhutch.org

# Install wget and curl
RUN apt-get update && \
	apt-get install -y wget curl

# Install DIAMOND v2.0.6
RUN mkdir /usr/diamond && cd /usr/diamond && \
	wget https://github.com/bbuchfink/diamond/releases/download/v2.0.6/diamond-linux64.tar.gz && \
	tar xzvf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz
