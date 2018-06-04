#!/bin/bash
set -e
set -u

# checksums
sashimi_md5="86e5924ecf8ce1272635ff43b244b32e"
sashimi_anno_md5="216c07785889074f69cb94cc2af7cb00"

pdfmd5() {
    grep -avE 'CreationDate|ModDate' $1 | md5sum | awk '{$0=$1}1'
}

fail() {
    echo ${1-""} >&2 && exit 1
}

files=( sashimi sashimi_anno )

anno=""
for f in ${files[@]}; do
    [[ $f == "sashimi_anno" ]] && anno="-g examples/annotation.gtf"
    docker run --rm -w $PWD -v $PWD:$PWD guigolab/ggsashimi $anno -b examples/input_bams.tsv -c chr10:27040584-27048100 -o ci/$f
    md5=$(pdfmd5 ci/$f.pdf)
    [[ $md5 == $(eval 'echo $'$f'_md5') ]] || fail "== Wrong checksum for $f.pdf: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0