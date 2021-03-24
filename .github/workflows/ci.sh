#!/bin/bash
set -e
set -u

fail() {
    echo ${1-""} >&2 && exit 1
}

# define global variables
bams="examples/input_bams.tsv"
region="chr10:27040584-27048100"

modes=(
    sashimi
    sashimi_anno
    sashimi_color
    sashimi_aggr
)

for mode in ${modes[@]}; do
    anno=""
    color=""
    aggr=""
    case $mode in
    sashimi_anno)
        sashimi_md5=(
            "c824e33480b828b0258eac90e1420991" # R 3.3.2 - ggplot2 2.2.1
            "f1183e65b5b89af995af3ba2e6c306c3" # R 3.6.3 - ggplot2 3.3.3
            "b11dec12a89ecdd2039d475b747e82fa" # R 4.0.3 - ggplot2 3.3.3
        )
        anno="-g examples/annotation.gtf"
        ;;
    sashimi_color)
        sashimi_md5=(
            "e62d13587f71d74dc1df86bfd0aed2b7" # R 3.3.2 - ggplot2 2.2.1
            "532581f0e2287a32535e2c5b848eb59c" # R 3.6.3 - ggplot2 3.3.3
            "38681c67b6db83206dc879b69ce2be2f" # R 4.0.3 - ggplot2 3.3.3
        )
        color="-C 3"
        ;;
    sashimi_aggr)
        sashimi_md5=(
            "ee3c894ffc734d957013b6cfefdb2d01" # R 3.3.2 - ggplot2 2.2.1
            "88afa86dbedea4a335b3fc07f94b6cd6" # R 3.6.3 - ggplot2 3.3.3
            "4049559a67d8ab4b5e17e49ed285febc" # R 4.0.3 - ggplot2 3.3.3
        )
        aggr="-C 3 -O 3 -A mean_j"
        ;;
    *)
        sashimi_md5=(
            "8be57782f81217393db59786e16049c7" # R 3.3.2 - ggplot2 2.2.1
            "b5b9b33bfc4b51944236baea43096b06" # R 3.6.3 - ggplot2 3.3.3
            "15c609eda4319f5395b1db8ab747d6e9" # R 4.0.3 - ggplot2 3.3.3
        )
        ;;
    esac
    ./ggsashimi.py $anno -b $bams -c $region $color $aggr
    md5=$(sed '/^\/\(.\+Date\|Producer\)/d' sashimi.pdf | md5sum | awk '$0=$1')
    [[ " ${sashimi_md5[@]} " =~ " ${md5} " ]] || fail "== Wrong checksum for $mode mode: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0

