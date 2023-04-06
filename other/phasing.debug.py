#!/usr/bin/env python
import sys

def read_phase_file(read_name, infile):
	result= {}
	with open(infile) as IN:
		for line in IN:
			lines = line.strip().split()
			if lines[0] not in read_name:
				continue
			result[read_name[lines[0]]] = lines[1]
	return result

def read_names(infile):
	result = {}
	with open(infile) as IN:
		for line in IN:
			lines = line.strip().split()
			result[lines[0]] = lines[1]
	return result

def read_paf(read_name, infile):
	result = {}
	with open(infile) as IN:
		for line in IN:
			lines = line.strip().split()
			if int(lines[11]) <= 0:
				continue
			if lines[0] not in read_name:
				continue
			result[read_name[lines[0]]] = [lines[7], lines[8]]
	return result

r'''
read_file
m64011_190830_220126/85068381/ccs 1
m64011_190901_095311/177079360/ccs 2

paf file
read_map_to_ref

reads phasing file
m64011_190830_220126/3/ccs          0           
m64011_190830_220126/12/ccs         0           
m64011_190830_220126/14/ccs         m           
m64011_190830_220126/16/ccs         p           
m64011_190830_220126/18/ccs         p           
m64011_190830_220126/21/ccs         m
m64011_190830_220126/24/ccs         p

test.log
9350 true 178632 614 {9216, 8759, 9783, 8491, 9704...
'''

if len(sys.argv) != 5:
	print("python phasing.debug.py read_name paf_files phase_reads conmunicate_log")
	sys.exit(1)

read_name = read_names(sys.argv[1])
paf_files = read_paf(read_name, sys.argv[2])
phase_reads = read_phase_file(read_name, sys.argv[3])

with open(sys.argv[4]) as IN:
	for line in IN:
		lines = line.strip().split()
		phases = {}
		s = 1000000000
		e = 0

		new_reads = []
		for read in lines[6:]:
			read = read.strip('{').strip('}').strip(',')
			phase = phase_reads[read] if read in phase_reads else 'non'
			if phase not in phases:
				phases[phase] = 0
			phases[phase] += 1
			if int(read) < s:
				s = int(read)
			if int(read) > e:
				e = int(read)
			new_reads.append("%s.%s" % (read, phase))
		print("%s :%s" % (s, str(paf_files[str(s)])), "%s :%s" % (e, paf_files[str(e)]), list(sorted(phases.items(), reverse=True)), 
			lines[:6], sorted(new_reads, key = lambda x: int(x.split('.')[0])))
