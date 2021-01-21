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
            "131d3ad1bb27c91693c342fbc0667d48" # R 3.3.2 - ggplot2 2.2.1
            "e7f72a5d373c93b4b0cbfa4976e9bd15" # R 3.6.3 - ggplot2 3.3.3
            "5617a83e8b3ba6d7478c40fcff99f5d3" # R 4.0.3 - ggplot2 3.3.3
        )
        anno="-g examples/annotation.gtf"
        ;;
    sashimi_color)
        sashimi_md5=(
            "41fee9e7e51a13a09c6d068cc6cb93cf" # R 3.3.2 - ggplot2 2.2.1
            "bb186bcd74674e7c8027e0dcfe2eb1fe" # R 3.6.3 - ggplot2 3.3.3
            "2606cf532430534c15cfc3a5d0db3620" # R 4.0.3 - ggplot2 3.3.3
        )
        color="-C 3"
        ;;
    sashimi_aggr)
        sashimi_md5=(
            "50674d8812194bbb8f8f8ddab9238b3d" # R 3.3.2 - ggplot2 2.2.1
            "5a0342a6093284ec0b9f85b01acb72e4" # R 3.6.3 - ggplot2 3.3.3
            "69b2a3f5028b1d9d0675ca98967cf353" # R 4.0.3 - ggplot2 3.3.3
        )
        aggr="-C 3 -O 3 -A mean_j"
        ;;
    *)
        sashimi_md5=(
            "b8ee3a415a80030ce92de15f40054a40" # R 3.3.2 - ggplot2 2.2.1
            "567d245fcd2d3d02288976e18bdee6b0" # R 3.6.3 - ggplot2 3.3.3
            "62d039688abd9f082b44abf98f313ae2" # R 4.0.3 - ggplot2 3.3.3
        )
        ;;
    esac
    ./sashimi-plot.py $anno -b $bams -c $region $color $aggr
    md5=$(sed '/^\/\(.\+Date\|Producer\)/d' sashimi.pdf | md5sum | awk '$0=$1')
    [[ " ${sashimi_md5[@]} " =~ " ${md5} " ]] || fail "== Wrong checksum for $mode mode: $md5"
done

echo "== All checksums match"
echo "== DONE"

exit 0

