FROM ubuntu:20.04
MAINTAINER sminot@fredhutch.org

# Install DIAMOND v0.9.31
RUN mkdir /usr/diamond && cd /usr/diamond && \
	wget https://github.com/bbuchfink/diamond/releases/download/v0.9.31/diamond-linux64.tar.gz && \
	tar xzvf diamond-linux64.tar.gz && \
	mv diamond /usr/bin/ && \
	rm diamond-linux64.tar.gz
