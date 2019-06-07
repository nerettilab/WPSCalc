#!/usr/bin/env python

# This script calculates the Window Protection Score (WPS) for an aligned sequencing data file
# WPS of a window of size K is defined as the number of molecules spanning the window minus those with an end point within the window. WPS is assigned to the coordinate in the center of the window.
# Author: nicholas_skvir@brown.edu

#Script takes .BED or .BEDPE file as input
#Note that bed files are zero-indexed and the last base is exclusive (For example, the first 100 bases of a chromosome are defined as chromStart=0, chromEnd=100, and span the bases numbered 0-99.)

import sys, os
import argparse
import gzip
import math
import csv
import io
from os.path import basename

parser = argparse.ArgumentParser(description='Script requires a bedfile of the alignments')
parser.add_argument('--version', action='version', version='%(prog)s 0.1')
parser.add_argument('alignment_file', action= 'store', metavar='alignment_file', help='Enter bedfile to use for calculating window protection score. Example /files/myfile.bed')
parser.add_argument('genome', action= 'store', metavar='genome', default='hg38', help='Specify the genome to which the bamfile was aligned')
parser.add_argument('chromosome', action= 'store', metavar='chromosome', help='Specify chromosome')
parser.add_argument('start', action= 'store', metavar='start', help='Specify start coordinate')
parser.add_argument('end', action= 'store', metavar='end', help='Specify end coordinate')
parser.add_argument('window_size', action= 'store', metavar='window_size', default=50, help='Specify the size of the sliding window used to calculate window protection score.')
parser.add_argument('--pairedend', action= 'store', dest='pairedend', default= False, help='Designate this option for paired-end sequencing. Default False, change to \'True\'')
args = parser.parse_args()

# parameters
bedfile = args.alignment_file
genome = args.genome
chromosome = args.chromosome
start = args.start
end = args.end
window_size = int(args.window_size)
paired_end = args.pairedend

filename = basename(bedfile)
outfile = str(filename) + "_" + str(chromosome) + "_" + str(start) + "_" + str(end) + "_WPS.txt"

print("Performing run on " + str(filename))
##Get appropriate chromosome lengths from reference file based on the annotation used for aligning the bam file (hg19 or hg38)
if genome=='hg19':
	with open('./hg19_chr_lengths.txt') as c:
		chromosomes=[]
		for line in c:
			line=line.split()
			if line:
				line=[str(i) for i in line]
				chromosomes.append(line)
elif genome=='hg38':
	with open('./hg38_chr_lengths.txt') as c:
		chromosomes=[]
		for line in c:
			line=line.split()
			if line:
				line=[str(i) for i in line]
				chromosomes.append(line)
else:
	sys.exit('must specify either \'hg19\' or \'hg38\' as reference genome')

#Create dictionary with keys for each chromosome and values of length
chr_lengths={}
for a in chromosomes:
	chr_lengths[a[0]]=a[1]

print("Calculating WPS for " + str(chromosome) + " from " + str(start) + " to " + str(end) + "...")
global_interval=[0,int(chr_lengths.get(chromosome))]
print "Chromosome length of " + str(chromosome) + ": " + str(global_interval[1])

#Check whether query range is valid
if (int(start)>= int(global_interval[1]) or int(end) > int(global_interval[1])):
	sys.exit('Specified \'Start\' and \'End\' coordinates must fall within the length of ' + str(chromosome))

#Read in file
with open(bedfile, 'rb') as file:
    csvin = csv.reader(file, delimiter='\t')
    lines = list(csvin)

#Add reads to list if they belong to the specified chromosome
data = []
for line in lines:
	if str(line[0]).lower() == str(chromosome).lower():
		data.append(line)

print(str(len(data)) + " reads found that map to " + str(chromosome))

#Split input data into lists of chromosomes, start, and end coordinates... add reads to lists if they fall within the specified start/end range
chrs = []
starts = []
ends = []

#Check whether reading in single or paired end bed file and add reads accordingly
if paired_end:
	print("Reading paired end data...")
	for read in data:
		if (int(read[1]) >= int(start)-int(window_size) and int(read[1]) <= int(end)-int(window_size)):
			chrs.append(read[0])
			starts.append(read[1])
			ends.append(read[5])

elif not paired_end:
	print("Reading unpaired data...")
	for read in data:
		if (int(read[1]) >= int(start)-int(window_size) and int(read[1]) <= int(end)-int(window_size)):
			chrs.append(read[0])
			starts.append(read[1])
			ends.append(read[2])

print((str(len(starts))) + " reads found within query range...")

#Create output file
with io.open(outfile, 'w', encoding='utf-8') as f:
       f.write(u"Chromosome" + u"\t" + u"Position" + u"\t" + u"Score" + u"\n")

#Define initial window midpoint, start point, and end point based on window size and global range over which the data points can occur
#Midpoint will start several positions in, depending on overall window size
#If window size is an even number, the midpoint will be the first position plus half the window size (if the window size is 6, and the first point is 0, the 'midpoint' will be at 3) and the window will be 1 longer on the left side than the right
#If the window size is an odd number, the midpoint will be the actual midpoint
if window_size%2==0:
	start_point=int(start)
	lower_bound=start_point-(window_size/2)
	upper_bound=start_point+(window_size/2)-1
	end_point=int(end)
else:
	start_point=int(start)
	lower_bound=start_point-(window_size//float(2))
	upper_bound=start_point+(window_size//float(2))
	end_point=int(end)


#print "Window Size: " + str(window_size)
#print "Starting Midpoint: "+ str(start_point)
#print "Starting Lower Bound: "+ str(lower_bound)
#print "Starting Upper Bound: "+ str(upper_bound)
#print "Final Midpoint: "+ str(end_point)


##FOR each position p in range of the midpoints to be checked...
WPS={}
for p in range(start_point, (end_point)):
	p_score=0

##FOR each interval on the current chromosome, IF the start point occurs before the start of the window, AND the end point occurs after the end of the window
	for v in range(0,len(starts)):
		if int(starts[v]) <= int(lower_bound) and int(ends[v]) > int(upper_bound):
##Increment the score at position p by +1
			p_score=p_score+1
			#print("ADDING 1 point at " + str(p) + " for data point " + str(v) + " on " + str(keys))
			#print("because " + str(v[0]) + " is less than or equal to " + str(lower_bound) + " and " + str(v[1]) + " is greater than " + str(upper_bound))

##If any start points or any end points occur between start and end of window
		elif (int(starts[v]) >= int(lower_bound) and int(starts[v]) <= int(upper_bound)) or (int(ends[v]) > int(lower_bound) and int(ends[v]) <= int(upper_bound)):
##Decrease the score at position p by -1
			p_score=p_score-1
			#print("REMOVING 1 point at " + str(p) + " for data point " + str(v) + " on " + str(keys))
			#print("because " + str(v[0]) + " is greater than " + str(lower_bound) + " and less than or equal to " + str(upper_bound) + " OR " + str(v[1]) + " is greater than " + str(lower_bound) + " and less than or equal to " + str(upper_bound))

##Increment the lower and upper bounds of the window by 1 and record the window protection score at in a dictionary with key corresponding to point p
	lower_bound=lower_bound+1
	upper_bound=upper_bound+1
	WPS[p]=p_score
	with io.open(outfile, 'a', encoding='utf-8') as out:
		out.write((str(chromosome)) + u" " + str(p) + u" " + str(WPS[p]) + u"\n")
