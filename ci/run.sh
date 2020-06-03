#!/bin/bash
set -e
set -u

# checksums
sashimi_md5="ff7d755c26673831db3fd0dd224b8bc0"
sashimi_anno_md5="e817cbf59608d090cf8292290dadc1cb"

fail() {
    echo ${1-""} >&2 && exit 1
}

export GGSASHIMI_DEBUG=1

modes=( sashimi sashimi_anno )

anno=""
for m in ${modes[@]}; do
    [[ $m == "sashimi_anno" ]] && anno="-g examples/annotation.gtf"
    ./sashimi-plot.py $anno -b examples/input_bams.tsv -c chr10:27040584-27048100
    md5=$(md5sum R_script | awk '$0=$1')
    [[ $md5 == $(eval 'echo $'$m'_md5') ]] || fail "== Wrong checksum for $m mode: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0

