## ====================
##  ggsashimi examples
## ====================

## Example #1. Overlay, intron shrinkage, gene annotation, PDF output, custom size and colors
../sashimi-plot.py -b input_bams.tsv -c chr10:27040584-27048100 -g annotation.gtf -M 10 -C 3 -O 3 --shrink --alpha 0.25 --base-size=20 --ann-height=4 --height=3 --width=18 -P palette.txt

## Example #2. Median coverage and number of reads supporting inclusion and exclusion, no gene annotation, TIFF output (350 PPI), custom size, default colors
../sashimi-plot.py -b input_bams.tsv -c chr10:27040584-27048100 -M 10 -C 3 -O 3 -A median --alpha 1 -F tiff -R 350 --base-size=16 --height=3 --width=18


