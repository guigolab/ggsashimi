#!/bin/bash
set -e
set -u

# checksums
sashimi_md5="b9eae97d095ab04724c521a5ba1edd5b"
sashimi_anno_md5="38a99516895dc3c1ab58aafbea845a4c"
sashimi_color_md5="b03de03e5465707041ea90583c475cb0"

fail() {
    echo ${1-""} >&2 && exit 1
}

# export GGSASHIMI_DEBUG=1

modes=( sashimi sashimi_anno sashimi_color )

for m in ${modes[@]}; do
    [[ $m == "sashimi_anno" ]] && anno="-g examples/annotation.gtf" || anno=""
    [[ $m == "sashimi_color" ]] && color="-C 3" || color=""
    ./sashimi-plot.py $anno -b examples/input_bams.tsv -c chr10:27040584-27048100 $color
    md5=$(sed '/^\/\(.\+Date\|Producer\)/d' sashimi.pdf | md5sum | awk '$0=$1')
    [[ $md5 == $(eval 'echo $'$m'_md5') ]] || fail "== Wrong checksum for $m mode: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0

