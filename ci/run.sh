#!/bin/bash
set -e
set -u

# checksums
sashimi_md5="d11062f49b37fc2f54faa36382a5380f"
sashimi_anno_md5="c6f8411283bcc6eb80f68da139af3505"

pdfmd5() {
    grep -avE 'CreationDate|ModDate|Producer' $1 | md5sum | awk '{$0=$1}1'
}

fail() {
    echo ${1-""} >&2 && exit 1
}

files=( sashimi sashimi_anno )

anno=""
for f in ${files[@]}; do
    [[ $f == "sashimi_anno" ]] && anno="-g examples/annotation.gtf"
    ./sashimi-plot.py $anno -b examples/input_bams.tsv -c chr10:27040584-27048100 -o ci/$f
    md5=$(pdfmd5 ci/$f.pdf)
    [[ $md5 == $(eval 'echo $'$f'_md5') ]] || fail "== Wrong checksum for $f.pdf: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0

