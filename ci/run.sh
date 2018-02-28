#!/bin/bash
set -e
set -u

# checksums
MD5SUM="19e50b22e5ee88314ad6a6630d202277"
MD5SUM_ANNO="da002cc4c9c4f2c77e4401c97564be94"

# run ggsashimi without annotation
docker run --rm -w $PWD -v $PWD:$PWD guigolab/ggsashimi -b examples/input_bams.tsv -c chr10:27040584-27048100 -o ci/sashimi
[[ $(grep -avE 'CreationDate|ModDate' ci/sashimi.pdf | md5sum | awk '{$0=$1}1') == $MD5SUM ]]

# run ggsashimi with annotation
docker run --rm -w $PWD -v $PWD:$PWD guigolab/ggsashimi -g examples/annotation.gtf -b examples/input_bams.tsv -c chr10:27040584-27048100 -o ci/sashimi-anno
[[ $(grep -avE 'CreationDate|ModDate' ci/sashimi-anno.pdf | md5sum | awk '{$0=$1}1') == $MD5SUM_ANNO ]]