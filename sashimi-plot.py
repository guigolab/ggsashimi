#!/usr/bin/env python

# Import modules
from argparse import ArgumentParser
#from subprocess import Popen, PIPE
import subprocess as sp
import sys, re, copy


def define_options():
	# Argument parsing
	parser = ArgumentParser(description='Create sashimi plot for a given genomic region')
	parser.add_argument("-b", "--bam", type=str,
		help="Individual bam file or file with a list of bam files and ids")
	parser.add_argument("-c", "--coordinates", type=str,
		help="Genomic region. Format: chr:start-end. Remember that bam coordinates are 0-based")
	parser.add_argument("-o", "--out-prefix", type=str, dest="out_prefix", default="sashimi",
		help="Prefix for plot file name [default=%(default)s]")
	parser.add_argument("-S", "--out-strand", type=str, dest="out_strand", default="both",
		help="Only for --strand other than 'NONE'. Choose which signal strand to plot: <both> <plus> <minus> [default=%(default)s]")
	parser.add_argument("-M", "--min-coverage", type=int, default=1, dest="min_coverage",
		help="Minimum number of reads supporting a junction to be drawn [default=1]")
	parser.add_argument("-g", "--gtf",
		help="Gtf file with annotation (only exons is enough)")
	parser.add_argument("-s", "--strand", default="NONE", type=str,
		help="Strand specificity: <NONE> <SENSE> <ANTISENSE> <MATE1_SENSE> <MATE2_SENSE> [default=%(default)s]")
	parser.add_argument("--shrink", action="store_true",
		help="Shrink the junctions by a factor for nicer display [default=%(default)s]")
	parser.add_argument("-O", "--overlay", type=int,
		help="Index of column with overlay levels (1-based)")
	parser.add_argument("-C", "--color-factor", type=int, dest="color_factor",
		help="Index of column with color levels (1-based)")
	parser.add_argument("-L", "--labels", type=int, dest="labels", default=1,
		help="Index of column with labels (1-based) [default=%(default)s]")
	parser.add_argument("--height", type=int, default=6,
		help="Height of the individual signal plot in inches [default=%(default)s]")
	parser.add_argument("--ann-height", type=float, default=1.5, dest="ann_height",
		help="Height of annotation plot in inches [default=%(default)s]")
	parser.add_argument("--width", type=int, default=10,
		help="Width of the plot in inches [default=%(default)s]")
	parser.add_argument("--base-size", type=int, default=14, dest="base_size",
		help="Base character size of the plot in pch [default=%(default)s]")
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


def flip_read(s, samflag):
	if s == "NONE" or s == "SENSE":
		return 0
	if s == "ANTISENSE":
		return 1
	if s == "MATE1_SENSE":
		if int(samflag) & 64:
			return 0
		if int(samflag) & 128:
			return 1
	if s == "MATE2_SENSE":
		if int(samflag) & 64:
			return 1
		if int(samflag) & 128:
			return 0


def read_bam(f, c, s):

	chr, start, end = parse_coordinates(c)

	# Initialize coverage array and junction dict
	a = {"+" : [0] * (end - start)}
	junctions = {"+": {}}
	if s != "NONE":
		a["-"] = [0] * (end - start)
		junctions["-"] = {}

	p = sp.Popen("samtools view %s %s " %(f, c), shell=True, stdout=sp.PIPE)
	for line in p.communicate()[0].strip().split("\n"):

		if line == "":
			continue

		line_sp = line.strip().split("\t")
		samflag, read_start, CIGAR = line_sp[1], int(line_sp[3]), line_sp[5]

		# Ignore reads with more exotic CIGAR operators
		if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
			continue

		read_strand = ["+", "-"][flip_read(s, samflag) ^ bool(int(samflag) & 16)]
		if s == "NONE": read_strand = "+"

		CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
		CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

		pos = read_start

		for n, CIGAR_op in enumerate(CIGAR_ops):
			CIGAR_len = int(CIGAR_lens[n])
			pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a[read_strand], junctions[read_strand], line=line)

	p.stdout.close()

	if a == {"+" : [0] * (end - start)}:
		print "WARNING: No reads in the specified area."

	return a, junctions


def read_bam_input(f, overlay, color, label):
	if f.endswith(".bam"):
		bn = f.strip().split("/")[-1].strip(".bam")
		yield bn, f, None, None, None
		return
	with open(f) as openf:
		for line in openf:
			line_sp = line.strip().split("\t")
			overlay_level = line_sp[overlay-1] if overlay else None
			color_level = line_sp[color-1] if color else None
			label_text = line_sp[label-1] if label else None
			yield line_sp[0], line_sp[1], '"%s"' %(overlay_level), '"%s"' %(color_level), label_text


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


def intersect_introns(data):
	data = sorted(data)
	it = iter(data)
	a, b = next(it)
	for c, d in it:
		if b > c:  # Use `if b > c` if you want (1,2), (2,3) not to be
			        # treated as intersection.
			b = min(b, d)
			a = max(a, c)
		else:
			yield a, b
			a, b = c, d
	yield a, b


def shrink_annotation(ann):
	return


def shrink_density(x, y, introns):
	new_x, new_y = [], []
	shift = 0
	start = 0
	# introns are already sorted by coordinates
	for a,b in introns:
		end = x.index(a)+1
		new_x += [int(i-shift) for i in x[start:end]]
		new_y += y[start:end]
		start = x.index(b)
		l = (b-a)
		shift += l-l**0.7
	new_x += [i-shift for i in x[start:]]
	new_y += y[start:]
	return new_x, new_y

def shrink_junctions(dons, accs, introns):
	new_dons, new_accs = [0]*len(dons), [0]*len(accs)
	shift_acc = 0
	shift_don = 0
	s = set()
	junctions = zip(dons, accs)
	for a,b in introns:
		l = b - a
		shift_acc += l-int(l**0.7)
		for i, (don, acc) in enumerate(junctions):
			if a >= don and b <= acc:
				if (don,acc) not in s:
					new_dons[i] = don - shift_don
					new_accs[i] = acc - shift_acc
				else:
					new_accs[i] = acc - shift_acc
				s.add((don,acc))
		shift_don = shift_acc
	return new_dons, new_accs

def read_gtf(f, c):
	exons = {}
	transcripts = {}
	chr, start, end = parse_coordinates(c)
	end = end -1
	with open(f) as openf:
		for line in openf:
			if line.startswith("#"):
				continue
			el_chr, ann, el, el_start, el_end, score1, strand, score2, tags = line.strip().split("\t")
			if el_chr != chr:
				continue
			d = dict(kv.strip().split(" ") for kv in tags.strip(";").split("; "))
			transcript_id = d["transcript_id"]
			el_start, el_end = int(el_start) -1, int(el_end)
			strand = '"' + strand + '"'
			if el == "transcript":
				if (el_end > start and el_start < end):
					transcripts[transcript_id] = max(start, el_start), min(end, el_end), strand
				continue
			if el == "exon":
				if (start < el_start < end or start < el_end < end):
					exons.setdefault(transcript_id, []).append((max(el_start, start), min(end, el_end), strand))

	return transcripts, exons


def make_introns(transcripts, exons, intersected_introns=None):
	new_transcripts = copy.deepcopy(transcripts)
	new_exons = copy.deepcopy(exons)
	introns = {}
	if intersected_introns:
		for tx, (tx_start,tx_end,strand) in new_transcripts.iteritems():
			total_shift = 0
			for a,b in intersected_introns:
				l = b - a
				shift = l - int(l**0.7)
				total_shift += shift
				for i, (exon_start,exon_end,strand) in enumerate(exons[tx]):
					new_exon_start, new_exon_end = new_exons[tx][i][:2]
					if a < exon_start:
						if b > exon_end:
							if i ==  len(exons[tx])-1:
								total_shift = total_shift - shift + (exon_start - a)*(1-int(l**-0.3))
							shift = (exon_start - a)*(1-int(l**-0.3))
							new_exon_end = new_exons[tx][i][1] - shift
						new_exon_start = new_exons[tx][i][0] - shift
					if b <= exon_end:
						new_exon_end = new_exons[tx][i][1] - shift
					new_exons[tx][i] = (new_exon_start,new_exon_end,strand)
			tx_start = min(tx_start, sorted(new_exons[tx])[0][0])
			new_transcripts[tx] = (tx_start, tx_end - total_shift, strand)

	for tx, (tx_start,tx_end,strand) in new_transcripts.iteritems():
		intron_start = tx_start
		ex_end = 0
		for ex_start, ex_end, strand in sorted(new_exons.get(tx, [])):
			intron_end = ex_start
			if tx_start < ex_start:
				introns.setdefault(tx, []).append((intron_start, intron_end, strand))
			intron_start = ex_end
		if tx_end > ex_end:
			introns.setdefault(tx, []).append((intron_start, tx_end, strand))
	d = {'transcripts': new_transcripts, 'exons': new_exons, 'introns': introns}
	return d


def gtf_for_ggplot(annotation, start, end, arrow_bins):
	arrow_space = (end - start)/arrow_bins
	s = """

	# data table with exons
	ann_list = list(
		"exons" = data.table(),
		"introns" = data.table()
	)
	"""

	if annotation["exons"]:
	
		s += """
		ann_list[['exons']] = data.table(
			tx = rep(c(%(tx_exons)s), c(%(n_exons)s)),
			start = c(%(exon_start)s),
			end = c(%(exon_end)s),
			strand = c(%(strand)s)
		)
		""" %({
		"tx_exons": ",".join(annotation["exons"].keys()),
		"n_exons": ",".join(map(str, map(len, annotation["exons"].itervalues()))),
		"exon_start" : ",".join(map(str, (v[0] for vs in annotation["exons"].itervalues() for v in vs))),
		"exon_end" : ",".join(map(str, (v[1] for vs in annotation["exons"].itervalues() for v in vs))),
		"strand" : ",".join(map(str, (v[2] for vs in annotation["exons"].itervalues() for v in vs))),
		})

	if annotation["introns"]:
	
		s += """
		ann_list[['introns']] = data.table(
			tx = rep(c(%(tx_introns)s), c(%(n_introns)s)),
			start = c(%(intron_start)s),
			end = c(%(intron_end)s),
			strand = c(%(strand)s)
		)
		# Create data table for strand arrows
		txarrows = data.table()
		introns = ann_list[['introns']]
		# Add right-pointing arrows for plus strand
		if ("+" %%in%% introns$strand) {
			txarrows = rbind(
				txarrows,
				introns[strand=="+" & end-start>5, list(
					seq(start+4,end,by=%(arrow_space)s)-1,
					seq(start+4,end,by=%(arrow_space)s)
					), by=.(tx,start,end)
				]
			)
		}
		# Add left-pointing arrows for minus strand
		if ("-" %%in%% introns$strand) {
			txarrows = rbind(
				txarrows,
				introns[strand=="-" & end-start>5, list(
					seq(start,max(start+1, end-4), by=%(arrow_space)s),
					seq(start,max(start+1, end-4), by=%(arrow_space)s)-1
					), by=.(tx,start,end)
				]
			)
		}
		""" %({
			"tx_introns": ",".join(annotation["introns"].keys()),
			"n_introns": ",".join(map(str, map(len, annotation["introns"].itervalues()))),
			"intron_start" : ",".join(map(str, (v[0] for vs in annotation["introns"].itervalues() for v in vs))),
			"intron_end" : ",".join(map(str, (v[1] for vs in annotation["introns"].itervalues() for v in vs))),
			"strand" : ",".join(map(str, (v[2] for vs in annotation["introns"].itervalues() for v in vs))),
			"arrow_space" : arrow_space,
		})

	s += """

	gtfp = ggplot()
	if (length(ann_list[['introns']]) > 0) {
		gtfp = gtfp + geom_segment(data=ann_list[['introns']], aes(x=start, xend=end, y=tx, yend=tx), size=0.3)
		gtfp = gtfp + geom_segment(data=txarrows, aes(x=V1,xend=V2,y=tx,yend=tx), arrow=arrow(length=unit(0.02,"npc")))
	}
	if (length(ann_list[['exons']]) > 0) {
		gtfp = gtfp + geom_segment(data=ann_list[['exons']], aes(x=start, xend=end, y=tx, yend=tx), size=5, alpha=1)
	}
	gtfp = gtfp + scale_y_discrete(expand=c(0,0.5))
	gtfp = gtfp + scale_x_continuous(expand=c(0,0.25), limits = c(%s,%s))
#	gtfp = gtfp + labs(y="m") + theme(axis.title.y=element_text(colour="white"))
	gtfp = gtfp + labs(y=NULL)
	gtfp = gtfp + theme(plot.margin=unit(c(0,0,0,0), "npc"))
	""" %(start, end)
	
	return s


def setup_R_script(h, w, b, label_dict):
	s = """
	library(ggplot2)
	library(grid)
	library(data.table)

	scale_lwd = function(r) {
		lmin = 0.1
		lmax = 4
		return( r*(lmax-lmin)+lmin )
	}

	height = %(h)s
	width = %(w)s
	base_size = %(b)s
	theme_set(theme_bw(base_size=base_size))
	theme_update(
		panel.grid = element_blank(),
		plot.margin = unit.c(unit(0,"npc"), unit(0,"npc"), unit(0,"npc"), unit(1,"strwidth","FbGN00000")),
		axis.line = element_line(size=0.5),
#		axis.text.y = element_blank(),
		axis.title.x = element_blank()
	)

	labels = list(%(labels)s)

	density_list = list()
	junction_list = list()

	""" %({
		'h': h,
		'w': w,
		'b': b,
		'labels': ",".join(('"%s"="%s"' %(id,lab) for id,lab in label_dict.iteritems())),
	})
	return s

def density_overlay(d, R_list):
#	lapply(names(l), function(x) cbind(l[[`x`]], id=x))
#	setNames(lapply(levels(as.factor(names(v))), function(y) {rbindlist(lapply(v[which(names(v)==y)], function(x) d[[as.character(x)]]))}), levels(as.factor(names(v))))
	s = """
	f = data.frame(id=c(%(id)s), fac=rep(c(%(levels)s), c(%(length)s)))
	%(R_list)s = setNames(
		lapply(
			levels(f$fac), function(y) {
				rbindlist(lapply(subset(f, fac==y)$id, function(x) %(R_list)s[[as.character(x)]]))
			}
		),
		levels(f$fac)
	)
	""" %({
		"levels": ",".join(d.keys()),
		"id": ",".join(map(str, ('"%s"' %(v) for vs in d.itervalues() for v in vs))),
		"length": ",".join(map(str, map(len, d.values()))),
		"R_list": R_list,
	})
	return s


def plot(R_script):
	p = sp.Popen("R --vanilla --slave", shell=True, stdin=sp.PIPE)
	p.communicate(input=R_script)
	p.stdin.close()
	p.wait()
	return


def colorize(d, p, color_factor):
	levels = sorted(set(d.itervalues()))
	n = len(levels)
	if n > len(p):
		p = (p*n)[:n]
	if color_factor:
		s = "color_list = list(%s)\n" %( ",".join('%s="%s"' %(k, p[levels.index(v)]) for k,v in d.iteritems()) )
	else:
		s = "color_list = list(%s)\n" %( ",".join('%s="%s"' %(k, "grey") for k,v in d.iteritems()) )
	return s



if __name__ == "__main__":

	parser = define_options()
	if len(sys.argv)==1:
	    parser.print_help()
	    sys.exit(1)
	args = parser.parse_args()

#	args.coordinates = "chrX:9609491-9612406"
	args.coordinates = "chrX:9609491-9610000"
#	args.bam = "/nfs/no_backup/rg/epalumbo/projects/tg/work/8b/8b0ac8705f37fd772a06ab7db89f6b/2A_m4_n10_toGenome.bam"


	strand_dict = {"plus": "+", "minus": "-"}
	bam_dict, overlay_dict, color_dict, id_list, label_dict = {"+":{}}, {}, {}, [], {}
	if args.strand != "NONE": bam_dict["-"] = {}

	
	for id, bam, overlay_level, color_level, label_text in read_bam_input(args.bam, args.overlay, args.color_factor, args.labels):
		id_list.append(id)
		label_dict[id] = label_text
		a, junctions = read_bam(bam, args.coordinates, args.strand)
		for strand in a:
		 	bam_dict[strand][id] = prepare_for_R(a[strand], junctions[strand], args.coordinates, args.min_coverage)
		if color_level is None:
			color_dict.setdefault(id, id)
		if overlay_level != '"None"':
			overlay_dict.setdefault(overlay_level, []).append(id)
			if color_level:
				color_dict.setdefault(overlay_level, overlay_level)
		if overlay_level == '"None"':
			color_dict.setdefault(id, color_level)


	if args.gtf:
		transcripts, exons = read_gtf(args.gtf, args.coordinates)


	for strand in bam_dict:

		# Output file name		
		out_prefix = args.out_prefix + "_" + strand
		if args.strand == "NONE":
			out_prefix = args.out_prefix
		else:
			if args.out_strand != "both" and strand != strand_dict[args.out_strand]:
				continue
				
		intersected_introns = None

		if args.shrink:
			introns = (v for vs in bam_dict[strand].values() for v in zip(vs[2], vs[3]))
			intersected_introns = list(intersect_introns(introns))

		if args.gtf:
		    annotation = make_introns(transcripts, exons, intersected_introns)

		# Define plot height
		bam_height = args.height * len(id_list)
		if args.overlay:
			bam_height = args.height * len(overlay_dict)
		if args.gtf:
			bam_height += args.ann_height
		R_script = setup_R_script(bam_height, args.width, args.base_size, label_dict)
		palette = "#ff0000", "#000000", "#00ff00"

		R_script += colorize(color_dict, palette, args.color_factor)


		for i, k in enumerate(id_list):
			v = bam_dict[strand][k]
			x, y, dons, accs, yd, ya, counts = v
			if args.shrink:
				x, y = shrink_density(x, y, intersected_introns)
				dons, accs = shrink_junctions(dons, accs, intersected_introns)
#				dons, accs, yd, ya, counts = [], [], [], [], []

			if i == 0:
				arrow_bins = 50
				if args.gtf:
					R_script += gtf_for_ggplot(annotation, x[0], x[-1], arrow_bins)


			R_script += """
			
			density_list[["%(id)s"]] = data.frame(x=c(%(x)s), y=c(%(y)s))
			junction_list[["%(id)s"]] = data.frame(x=c(%(dons)s), xend=c(%(accs)s), y=c(%(yd)s), yend=c(%(ya)s), count=c(%(counts)s))
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

		if args.overlay:
			R_script += density_overlay(overlay_dict, "density_list")
			R_script += density_overlay(overlay_dict, "junction_list")

		R_script += """

		pdf("%(out)s", h=height, w=width)
		grid.newpage()
		pushViewport(viewport(
			layout = grid.layout(
				nrow=length(density_list)+%(args.gtf)s, 
				ncol=1, 
				heights=unit(c(rep(%(signal_height)s,length(density_list)), %(ann_height)s*%(args.gtf)s), "in")
				)
			)
		)

		for (bam_index in 1:length(density_list)) {

			id = names(density_list)[bam_index]
			d = density_list[[id]]
			junctions = junction_list[[id]]

			maxheight = max(d[['y']])

			# Density plot
			gp = ggplot(d) + geom_bar(aes(x, y), position='identity', stat='identity', fill=color_list[[id]], alpha=1/2)
			gp = gp + labs(y=labels[[id]])
			gp = gp + theme(axis.text.y=element_blank())
			gp = gp + scale_x_continuous(expand=c(0,0.2))
#			gp = gp + scale_y_continuous(expand=c(0,0))
			if (bam_index != length(density_list)) {
				gp = gp + theme(axis.text.x = element_blank())
			}



			if (nrow(junctions)>0) {row_i = 1:nrow(junctions)} else {row_i = c()}


			for (i in row_i) {

				j_tot_counts = sum(junctions[['count']])

				j = as.numeric(junctions[i,])

				# Find intron midpoint
				xmid = round(mean(j[1:2]), 1)
				ymid = max(j[3:4]) * 1.1

				# Thickness of the arch
				lwd = scale_lwd(j[5]/j_tot_counts)

				curve_par = gpar(lwd=lwd, col=color_list[[id]])

				# Choose position of the arch (top or bottom)
#				nss = sum(junctions[,1] %%in%% j[1])
#				nss = i
				nss = 1
				if (nss%%%%2 == 0) {  #bottom
					ymid = -0.4 * maxheight
					# Draw the arcs
					# Left
					curve = xsplineGrob(x=c(0, 0, 1, 1), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
					gp = gp + annotation_custom(grob = curve, j[1], xmid, 0, ymid)
					# Right
					curve = xsplineGrob(x=c(1, 1, 0, 0), y=c(1, 0, 0, 0), shape=1, gp=curve_par)
					gp = gp + annotation_custom(grob = curve, xmid, j[2], 0, ymid)
				}

				if (nss%%%%2 != 0) {  #top
					# Draw the arcs
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

		if (%(args.gtf)s == 1) {
			print(gtfp, vp=viewport(layout.pos.row = bam_index+1, layout.pos.col = 1))
		}


		dev.off()

		""" %({
			"out": "%s.pdf" %out_prefix,
			"args.gtf": float(bool(args.gtf)),
			"signal_height": args.height,
			"ann_height": args.ann_height,
			})


		plot(R_script)
	exit()






