#!/usr/bin/env python

# Import modules
from argparse import ArgumentParser
#from subprocess import Popen, PIPE
import subprocess as sp
import sys
import re






def parse_coordinates(c):
	chr = c.split(":")[0]
	start, end = c.split(":")[1].split("-")
	# Convert to 0-based 
	start, end = int(start) - 1, int(end)
	return chr, start, end



def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):

	# Match
	if CIGAR_op == "M":
		for i in range(pos, pos + CIGAR_len):
			if i < start or i >= end:
				continue
			ind = i - start
			a[ind] += 1

	# Insertion
	if CIGAR_op == "I":
		return pos

	# Deletion or Soft-clip
	if CIGAR_op == "D" or CIGAR_op == "S":
		pass

	# Junction
	if CIGAR_op == "N":
		don = pos
		acc = pos + CIGAR_len
		if don > start and acc < end:
			junctions[(don,acc)] = junctions.setdefault((don,acc), 0) + 1

	pos = pos + CIGAR_len

	return pos



def read_bam(f, c):

	chr, start, end = parse_coordinates(c)

	# Initialize coverage array 
	a = [0] * (end - start)

	# Junction dictionary
	junctions = dict()

	p = sp.Popen("samtools view %s %s " %(f, c), shell=True, stdout=sp.PIPE)
	for line in p.communicate()[0].strip().split("\n"):

		line_sp = line.strip().split("\t")
		samflag, read_start, CIGAR = line_sp[1], int(line_sp[3]), line_sp[5]

		# Ignore reads with more exotic CIGAR operators
		if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])): 
			continue

		CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
		CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

		pos = read_start

		for n, CIGAR_op in enumerate(CIGAR_ops):
			CIGAR_len = int(CIGAR_lens[n])
			pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions)

	p.stdout.close()
	
	return a, junctions


def prepare_for_R(a, junctions, c, m):

	chr, start, end = parse_coordinates(args.coordinates)

	# Convert the array index to genomic coordinates
	x = list(i+start for i in range(len(a)))
	y = a

	# Arrays for R
	dons, accs, yd, ya, counts = [], [], [], [], []

	# Prepare arrays for junctions (which will be the arcs)
	for (don, acc), n in junctions.iteritems():

		# Do not add junctions with less than defined coverage
		if n < m:
			continue

		dons.append(don)
		accs.append(acc)
		counts.append(n)

		yd.append( a[ don - start ])
		ya.append( a[ acc - start ])

	return x, y, dons, accs, yd, ya, counts


if __name__ == "__main__":

	# Argument parsing
	
	parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
	parser.add_argument("-b", "--bam", type=str, help="Bam file")
	parser.add_argument("-c", "--coordinates", type=str, help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
	parser.add_argument("-m", "--min_coverage", type=int, default=0, 
		help="Minimum number of reads supporting a junction to be drawn [default=0]")
#	parser.add_argument("-s", "--smooth", action="store_true", default=False, help="Smooth the signal histogram")
	args = parser.parse_args()
	
	args.coordinates = "chrX:9609491-9612406"
	args.bam = "/nfs/no_backup/rg/epalumbo/projects/tg/work/8b/8b0ac8705f37fd772a06ab7db89f6b/2A_m4_n10_toGenome.bam"
	
	a, junctions = read_bam(args.bam, args.coordinates)
	x, y, dons, accs, yd, ya, counts = prepare_for_R(a, junctions, args.coordinates, args.min_coverage)
	

	print """
	library(ggplot2)
	library(grid)
	d = data.frame(x=c(%(x)s), y=c(%(y)s))
	d[['y']] = as.double(smooth(d[['y']]))

	junctions = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s))

	gp = ggplot(d) + geom_bar(aes(x, y), position='identity', stat='identity')

	curve = function(angle, color, lwd, curvature=-0.3) {
		curveGrob(0, 1, 1, 0, default.units = "npc",
		curvature = curvature, ncp = 50, shape = -1,
		square = T, squareShape = 1,
		inflect = F, open = TRUE, 
		name = NULL, vp = NULL, debug = FALSE, 
		gp = gpar(col=color, lwd=lwd), angle = angle)
	}

	for (i in 1:nrow(junctions)) {
		j = as.numeric(junctions[i,])
		# Find intron midpoint 
		xmid = round(mean(j[1:2]), 1)
		ymid = round(mean(j[3:4]), 1)
		
		# Left curve
		curve = curveGrob(0, 0, 1, 1, default.units = "npc",
			curvature = -0.4, ncp = 100, shape = -1,
			square = T, squareShape = 0,
			inflect = F, open = TRUE, 
			name = NULL, vp = NULL, debug = FALSE, 
			gp = gpar(col="grey", lwd=1), angle = 120
		)
		gp = gp + annotation_custom(grob = curve, j[1], xmid, j[3], 20)

		# Right curve
		curve = curveGrob(1, 0, 0, 1, default.units = "npc",
			curvature = +0.4, ncp = 100, shape = -1,
			square = T, squareShape = 0,
			inflect = F, open = TRUE, 
			name = NULL, vp = NULL, debug = FALSE, 
			gp = gpar(col="grey", lwd=1), angle = 40
		)
		gp = gp + annotation_custom(grob = curve, xmid, j[2], 20, j[4])

#		my_grob = curve(90, "black", 1, -0.3)
#		gp = gp + annotation_custom(grob = my_grob, j[1], j[2], j[3], j[4])

	}


	ggsave('%(out)s', h=4, w=10)
	""" %({
		'x' : ",".join(map(str, x)),
		'y' : ",".join(map(str, y)),
		'dons' : ",".join(map(str, dons)),
		'accs' : ",".join(map(str, accs)),
		'yd' : ",".join(map(str, yd)),
		'ya' : ",".join(map(str, ya)),
		'out' : "tmp.pdf",
		})
	exit()






