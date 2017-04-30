#!/usr/bin/env python

# Import modules
from argparse import ArgumentParser
#from subprocess import Popen, PIPE
import subprocess as sp
import sys, re
import operator


def define_options():
	# Argument parsing
	parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
	parser.add_argument("-b", "--bam", type=str, 
		help="Individual bam file or file with a list of bam files and ids")
	parser.add_argument("-c", "--coordinates", type=str,
		help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
	parser.add_argument("-M", "--min_coverage", type=int, default=1, 
		help="Minimum number of reads supporting a junction to be drawn [default=1]")
	parser.add_argument("-g", "--gtf", 
		help="Gtf file with annotation (only exons is enough)")
#	parser.add_argument("-s", "--smooth", action="store_true", default=False, help="Smooth the signal histogram")
	return parser



def parse_coordinates(c):
	chr = c.split(":")[0]
	start, end = c.split(":")[1].split("-")
	# Convert to 0-based 
	start, end = int(start) - 1, int(end)
	return chr, start, end



def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions, line):

	# Match
	if CIGAR_op == "M":
		for i in range(pos, pos + CIGAR_len):
			if i < start or i >= end:
				continue
			ind = i - start
			a[ind] += 1

	# Insertion or Soft-clip
	if CIGAR_op == "I" or CIGAR_op == "S":
		return pos

	# Deletion 
	if CIGAR_op == "D":
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
			pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions, line=line)

	p.stdout.close()
	
	return a, junctions


def read_bam_input(f):
	if f.endswith(".bam"):
		bn = f.strip().split("/")[-1].strip(".bam")
		return [(bn, f)]
	with open(f) as openf:
		return list(line.strip().split("\t") for line in openf)


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

		yd.append( a[ don - start -1 ])
		ya.append( a[ acc - start +1 ])

	return x, y, dons, accs, yd, ya, counts


def read_gtf(f, c):
	exons = {}
	chr, start, end = parse_coordinates(c)
	with open(f) as openf:
		for line in openf:
			el_chr, ann, el, el_start, el_end, score1, strand, score2, tags = line.strip().split("\t")
			if el_chr != chr:
				continue
			if el != "exon":
				continue
			exon_start, exon_end = int(el_start), int(el_end)
			# Ignore elements not included in the region
			if not (start < exon_start < end or start < exon_end < end):
				continue

			d = dict(kv.split(" ") for kv in tags.strip(";").split("; "))
			transcript_id = d["transcript_id"]
			strand = '"' + strand + '"'
			exons.setdefault(transcript_id, []).append((max(exon_start, start), min(end, exon_end), strand))
	return exons


if __name__ == "__main__":

	parser = define_options()
	args = parser.parse_args()
	
	args.coordinates = "chrX:9609491-9612406"
#	args.coordinates = "chrX:9609491-9610000"
#	args.bam = "/nfs/no_backup/rg/epalumbo/projects/tg/work/8b/8b0ac8705f37fd772a06ab7db89f6b/2A_m4_n10_toGenome.bam"

	if args.gtf:
		annotation = read_gtf(args.gtf, args.coordinates)

	bam_dict = {}
	for id, bam in read_bam_input(args.bam):
		a, junctions = read_bam(bam, args.coordinates)
	 	bam_dict[id] = prepare_for_R(a, junctions, args.coordinates, args.min_coverage)
	
	print """
	library(ggplot2)
	library(grid)

	density_list = list()
	junction_list = list()
	"""
	
	if args.gtf:
		print """
		annotation = data.frame(
			tx = rep(c(%(tx)s), c(%(n_exons)s)), 
			exon_start = c(%(exon_start)s),
			exon_end = c(%(exon_end)s),
			strand = c(%(strand)s)
		)
		print(annotation)
		""" %({
			"tx": ",".join(annotation.iterkeys()),
			"n_exons": ",".join(map(str, map(len, annotation.itervalues()))),
			"exon_start" : ",".join(map(str, (v[0] for vs in annotation.itervalues() for v in vs))),
			"exon_end" : ",".join(map(str, (v[1] for vs in annotation.itervalues() for v in vs))),
			"strand" : ",".join(map(str, (v[2] for vs in annotation.itervalues() for v in vs))),
		})
		
	for k, v in bam_dict.iteritems():
		x, y, dons, accs, yd, ya, counts = v
		
		print """
		density_list$%(id)s = data.frame(x=c(%(x)s), y=c(%(y)s))
		junction_list$%(id)s = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
		""" %({
			"id": k,
			'x' : ",".join(map(str, x)),
			'y' : ",".join(map(str, y)),
			'dons' : ",".join(map(str, dons)),
			'accs' : ",".join(map(str, accs)),
			'yd' : ",".join(map(str, yd)),
			'ya' : ",".join(map(str, ya)),
			'counts' : ",".join(map(str, counts)),
		})

	print """	
	height = 4
	base_size = 14
	theme_set(theme_bw(base_size=base_size))
	theme_update(
		panel.grid = element_blank()
	)

	pdf('%s', h=height, w=10)
	grid.newpage()
	pushViewport(viewport(layout = grid.layout(length(density_list), 1)))

	for (bam_index in 1:length(density_list)) {
	
		id = names(density_list)[bam_index]
		d = density_list[[id]]
		junctions = junction_list[[id]]
	
		maxheight = max(d[['y']])
	
		scale_lwd = function(r) {
			lmin = 0.1
			lmax = 4
			return( r*(lmax-lmin)+lmin )
		}
	
	
		# Density plot
		gp = ggplot(d) + geom_bar(aes(x, y), position='identity', stat='identity')
	
		gp = gp + labs(title=id)

		j_tot_counts = sum(junctions[['count']])
	
		for (i in 1:nrow(junctions)) {
	
			j = as.numeric(junctions[i,])
	
			# Find intron midpoint 
			xmid = round(mean(j[1:2]), 1)
			ymid = max(j[3:4]) * 1.1
	
			# Thickness of the arch
			lwd = scale_lwd(j[5]/j_tot_counts)
	
			curve_par = gpar(lwd=lwd)
	
			# Choose position of the arch (top or bottom)
			nss = length(match(j[1], junctions[,1]))
			if (nss%%%%2 == 0) {  #bottom
				ymid = -0.4 * maxheight
				# Draw the archs
				# Left
				curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
				gp = gp + annotation_custom(grob = curve, j[1], xmid, 0, ymid)
				# Right
				curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
				gp = gp + annotation_custom(grob = curve, xmid, j[2], 0, ymid)
			} 
	
			if (nss%%%%2 != 0) {  #top
				# Draw the archs
				# Left
				curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
				gp = gp + annotation_custom(grob = curve, j[1], xmid, j[3], ymid)
				# Right
				curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(0, 1, 1, 1), shape=1, gp=curve_par)
				gp = gp + annotation_custom(grob = curve, xmid, j[2], j[4], ymid)
		
				gp = gp + annotate("label", x = xmid, y = ymid, label = j[5], 
					vjust=0.5, hjust=0.5, label.padding=unit(0.01, "lines"), 
					label.size=NA, size=(base_size*0.352777778)*0.6
				)
			}
	

	#		gp = gp + annotation_custom(grob = rectGrob(x=0, y=0, gp=gpar(col="red"), just=c("left","bottom")), xmid, j[2], j[4], ymid)
	#		gp = gp + annotation_custom(grob = rectGrob(x=0, y=0, gp=gpar(col="green"), just=c("left","bottom")), j[1], xmid, j[3], ymid)
	
	
		}
		print(gp, vp=viewport(layout.pos.row = bam_index, layout.pos.col = 1))

	}

	

	dev.off()
	""" %("tmp.pdf")

	exit()






