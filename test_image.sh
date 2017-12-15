#!/bin/bash

set -e

test_image(){
	
	img_tag=$1

	[[ "${#img_tag}" == "0" ]] && echo "Please specify image" && return

	[[ "$(docker run $img_tag echo True)" != "True" ]] && echo "Tag not found ($img_tag)" && return

	docker run \
		-v $PWD/temp:/share \
		-v $HOME:/root \
		--rm \
		$img_tag \
			run.py \
				--input sra://SRR2241114 \
				--ref-db s3://fh-pi-fredricks-d/lab/Sam_Minot/dbs/genbank_viral_diamond_20170911/genbank_viral.dmnd \
				--output-folder s3://fh-pi-fredricks-d/lab/Sam_Minot/data/test/ \
				--temp-folder /share \
				--overwrite
}

test_image $1

