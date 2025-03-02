#!/usr/bin/env python
# -*- coding: utf-8 -*-
# CREATED BY: Alexander Garcia-Zea
# jagarcia@eafit.edu.co

# Converts genbank to mochiview formats:
# Converts a genbank formatted genome to a MochiView format. However, it adds an AA sequence column

import sys
import os
########################################################################################################################
def collect_all_info(file_name, GET_FASTA):
	"""Uses a list of line numbers to identify lines corresponding to a given gene's information. 
	Collects gene location, name, and description as a list of lists like [[gene1, start, stop, description], [gene2..."""

	# Information value is an error by default unless correctly identified by this function.
	f = open(file_name, "r")
	lines = f.readlines()
	index = len(lines)
	annotation_switch = "no"
	features_switch = "no"
	genome_accession = file_name.split(".")[0]  # Use the filename without the extension
	annotation_start_line = 0

	genes = []
	line_number = -1

	for i in range(index):
		line_number += 1
		try:
			current_line = lines[line_number]
		except IndexError:
			continue

		if line_number == 0:
			genome_name = genome_accession  # Replace genome name with filename without extension
		else:
			if annotation_start_line == 0:  # Haven't reached annotations yet
				if annotation_switch == "no":  # Before "features" header, no gene info
					data = current_line[:21].strip(" ")
					if data == "FEATURES":  # Next actual feature name = start of gene info
						annotation_switch = "yes"
						features_switch = "yes"
					elif data == "LOCUS":  # "Locus" first line of doc, has systematic genome name
						splt = current_line.split(" ")
						clean_splt = []
						for i in splt:
							if i != '':
								clean_splt.append(i)

				else:
					if features_switch == "yes":
						tester = current_line[:21].strip(" ")
						if tester == "":
							continue
						else:
							if tester == "source":  # Skip first "source" entry
								continue
							else:
								annotation_start_line = line_number
			else:
				if current_line[:6] == "ORIGIN":  # FASTA DNA sequence started - toggle switch. Collect/print?
					last_data_line = line_number + 2
					if GET_FASTA == "Y":
						clean_fasta(last_data_line, genome_accession, lines)  # Write fasta data to outfile.
					break

	# From all data lines, return a list of features, each of which is a list of items		
	all_feature_info = collect_feature_data(lines, annotation_start_line, last_data_line)

	# Identify items in each field needed for output
	outdata = []
	for feature in all_feature_info:
		if len(feature) < 2:
			continue
		else:
			result = identify_data(feature)
			outdata.append(result)
	return genome_name, outdata

########################################################################################################################
# The rest of the script remains unchanged until create_data_lines function.
def collect_feature_data(all_lines, start_pos, end_pos):
	'''
	- input is all data lines, and gene/feature annotations start line number (zero indexed)
	- start_pos is for the whole file (start position of annotations)
	- Returns all gene/feature annotation data
	'''
	feature_data = []
	current_feature = []
	current_line = ""
	for k in range(start_pos, end_pos):
		clean_line = all_lines[k][21:].strip("\n").strip('"')
		current_line = current_line + clean_line	# Start with current line
		try:
			next_line = all_lines[k+1]
		except IndexError:							# End of data edge case
			current_feature.append(current_line)
			feature_data.append(current_feature)
		
		new_feature = next_line[:21].strip(" ")		# Next line is new feature?
		new_item = next_line[21:22]					# Next line is new item?
		if len(new_feature) == 0:					# Should just be if new_feature == "": But it didn't work for some reason.
			if new_item != "/":
				continue
			else:
				current_feature.append(current_line)
				current_line = ""
		else:
			current_feature.append(current_line)
			feature_data.append(current_feature)
			current_line = ""
			current_feature = []
	return feature_data

####################################################################################################
def identify_data(clean_data_list):
	'''
	For one feature, which is a list of items: 
	- Parses and identifies discrete data fields individual features list [data field 1, data field 2,...]
	- Returns identified data
	'''

	description = "None"
	systematic_name = "Error"
	pseudo_gene = "N"
	gene_name = "Unknown"
	alias = "None"
	translation = "No Data"
	product = "None"
	start = "No Data"
	end = "No Data"
	strand = "No Data"

	for l in clean_data_list:
		slash_detector = l[:1]
		if slash_detector == "/":
			# Type of data is evident from first few letters of the line
			splt = l.strip("/").split("=")
			data_ID = splt[0]
			try:
				data_actual = splt[1]
				if data_ID == "product":
					description = data_actual.strip("\n").strip('"')
				elif data_ID == 'locus_tag':
					systematic_name = data_actual.strip("\n").strip('"')
				elif data_ID == "pseudo" or data_ID == "pseudogene":
					pseudo_gene = "Y"
				elif data_ID == "mobile_element":
					description = data_actual.strip("\n").strip('"')
				elif data_ID == 'gene':
					gene_name = data_actual.strip("\n").strip('"')
				elif data_ID == "protein_id":
					alias = data_actual.strip("\n").strip('"')					
				elif data_ID == "translation":
					translation = l[14:].strip('"')
				else: 
					pass
			except IndexError:
				continue

		elif slash_detector != "/":
			l_fix1 = l.replace(">", "")
			l_clean = l_fix1.replace("<", "")
			if l_clean[:10] == "complement":
				strand = "-"
				comp_a = l_clean.strip("complement(").strip(")")
				if comp_a[:4] == "join":
					start, end, gene_name = join_cleaner(comp_a)
					break
				elif comp_a[:5] == "order":
					start, end, gene_name = order_cleaner(comp_a)
					break
				else:
					splt3 = comp_a.split("..")
					start = splt3[0]
					end = splt3[1]
			elif l_clean[:5] == "order":
				start, end, gene_name = order_cleaner(l_clean)
				break
			elif l_clean[:4] == "join":
				start, end, gene_name = join_cleaner(l_clean)
				break
			else:
				y = l_clean.split("..")
				strand = "+"
				try:
					test = int(y[0])			# Ensures that the result doesn't have additional formatting
					test2 = int(y[1])			# Same test, but for end number
					start = y[0]
					end = y[1]
				except ValueError:
					print("ERROR - gene coordinate formatting.")
					print("Line: ", l)
					
	if gene_name == "Unknown":
		gene_name = systematic_name

	organized_data = [start, end, strand, description, systematic_name, gene_name, pseudo_gene, alias, translation]
	return organized_data

########################################################################################################################
def join_cleaner(joined_coordinates):
	"""Certain features are annotated as 'join(123..543, 678..987)' and these need to be cleaned up. Returns
	start and stop coordinates, but for the whole thing. i.e. it joins all the individual pieces together,
	which throws away annotations indicating individual pieces. This process could be improved in future versions."""

	strp_a = joined_coordinates.strip("join(").strip(")")
	splt_a = strp_a.split(",")
	final_coords = []
	for i in splt_a:
		coordinates_ab = i.split("..")
		for j in coordinates_ab:
			final_coords.append(j)
	start_position = final_coords[0]
	stop_position = final_coords[-1]

	# Could add in some sort of annotation for individual joined regions within the name in the future
	gene_nm = ""

	return start_position, stop_position, gene_nm

########################################################################################################################
def order_cleaner(ordered_coordinates):
	"""Certain features are annotated as 'order(123..543, 678..987)' and these need to be cleaned up. Returns
	start and stop coordinates, but for the whole thing. i.e. it joins all the individual pieces together,
	which throws away annotations indicating individual pieces. This process could be improved in future versions."""

	strp_y = ordered_coordinates.strip("order(").strip(")")
	splt_1 = strp_y.split(",")
	final_coords_ordered = []
	for i in splt_1:
		coordinates_abc = i.split("..")
		for j in coordinates_abc:
			final_coords_ordered.append(j)
	start_position = final_coords_ordered[0]
	stop_position = final_coords_ordered[-1].strip(")")

	# Could add in some sort of annotation for individual joined regions within the name in the future
	gene_nm2 = ""

	return start_position, stop_position, gene_nm2

########################################################################################################################
def clean_fasta(firstlinenumber, accession, datalines):
	'''Infile = genbank data file with genome sequence. Using first and last line of genome sequence data, spits out 
	a cleaned fasta file.'''

	fasta = open(accession + ".fasta", "w")
	genome_seq_dirty = ""
	genome_seq_clean = ""
	bad_characters = ["0","1","2","3","4","5","6","7","8","9"," ", "/"]

	#write the genome name
	fasta.write(">" + accession + "\n")
	lastline = len(datalines)

	# collect genome data with numbers
	fasta_line = 0
	for line in datalines:
		fasta_line += 1
		if fasta_line >= firstlinenumber:
			l = line.replace(" ","")
			if l[:2] == "//":
				break
			else:
				genome_seq_dirty = genome_seq_dirty + l
	
	# clean the spaces and numbers
	for i in genome_seq_dirty:
		if i in bad_characters:
			continue
		else:
			genome_seq_clean = genome_seq_clean + i
	fasta.write(genome_seq_clean)
	return None

########################################################################################################################
def remove_duplicate_lines(semifinal_data):
	"""
	- Tests final data, then refines the list to remove duplicate lines (almost all are duplicated).
	- Preference is given to the SECOND entry, which usually has more info than the first
	- Otherwise duplicate lines with AA sequences are kept over those without AA sequences
	"""
	finaldata_list = []
	previous_line = []		# keeps track of previous line. It will be appended to the final data list or not depending
							# upon what the next line is. If its a duplicate, the lines will be compared and the good
							# one will be kept while the other is tossed. 
	previous_start = ""		# keeps track of duplicate start coords which can occur in successive lines
	previous_AA = ""		# keep strack of previous line's presence/absence of an AA sequence
	systematic_list = []
	datarange = len(semifinal_data)

	for x in range(datarange):
		data_point = semifinal_data[x]

		if x == datarange-1:
			if data_point[4] not in systematic_list:
				finaldata_list.append(data_point)				# Edge case, last line
		elif x == 0:
			previous_line = data_point					# Current line is now the "previous" line
			previous_start = data_point[0]
			previous_end = data_point[1]
			previous_AA = data_point[8]
		else:
			if data_point[0] == "No Data":						# Toss data if there is not start coordinate
				continue
			elif data_point[4] == "Error":						# Toss data point if there is no systematic name
				continue
			else:
				recorded = "no"
				if previous_line == []:							# Starting point edge case - no previous line, so just toggle to next line
					previous_line = data_point					# Current line is now the "previous" line
					previous_start = data_point[0]
					previous_end = data_point[1]
					previous_AA = data_point[8]
					recorded = "yes"
				elif previous_line != []:
					currentstart = data_point[0]
					currentend = data_point[1]
					
					if previous_start == currentstart or previous_end == currentend:	# Duplicate line found. Compare/Record only one.
						aa_entry = data_point[8]					# Current AA entry
						if aa_entry != "" and previous_AA == "":	# Second line has AA sequence, first does not. 99% of entries
							
							if data_point[4] not in systematic_list:
								finaldata_list.append(data_point)		# Record second line
								systematic_list.append(data_point[4])
							recorded = "yes"
							previous_line = data_point						# Reset
							previous_start = currentstart
							previous_AA = aa_entry
						else:
							if aa_entry == "" and previous_AA != "":		# Record first line. (Only time when first should be recorded)
								if data_point[4] not in systematic_list:
									finaldata_list.append(previous_line)
									systematic_list.append(data_point[4])
								recorded = "yes"
								previous_line = []							# Reset
								previous_start = ""
								previous_AA = ""
							else: 											# Record second line
								if data_point[4] not in systematic_list:
									finaldata_list.append(data_point)
									systematic_list.append(data_point[4])
								recorded = "yes"
								previous_line = []							# Reset
								previous_start = ""
								previous_AA = ""
						
					else: #previous_start != currentstart:		# Previous line is not duplicated. Record previous line, toggle current. 
						if data_point[4] not in systematic_list:
							finaldata_list.append(previous_line)
							systematic_list.append(data_point[4])
						previous_line = data_point				# Current line is now the "previous" line
						previous_start = data_point[0]
						previous_AA = data_point[8]
				if recorded == "no":							# Just record one line. Still not sure why this is needed [Rolls eyes]
					if data_point[4] not in systematic_list:
						finaldata_list.append(data_point)
						systematic_list.append(data_point[4])
	return finaldata_list

def create_data_lines(finaldata, infile_name, genome_nm, wantproteinYN):
	"""Creates data line in MochiView format, returns the line as a list with a single string entry."""

	out = open(infile_name.split(".")[0] + ".mochi", "w")
	if wantproteinYN == "Y":
		out.write("SEQ_NAME\tSTART\tEND\tSTRAND\tFEATURE_NAME\tGENE_NAME\tALIASES\tDESCRIPTION\tTXN_START\tTXN_END\tCDS_START\tCDS_END\tEXON_COUNT\tEXON_STARTS\tEXON_ENDS\tIS_PRIMARY\tAAseq\n")
	elif wantproteinYN == "N":
		out.write("SEQ_NAME\tSTART\tEND\tSTRAND\tFEATURE_NAME\tGENE_NAME\tALIASES\tDESCRIPTION\tTXN_START\tTXN_END\tCDS_START\tCDS_END\tEXON_COUNT\tEXON_STARTS\tEXON_ENDS\tIS_PRIMARY\n")
	# Order: sys (0), start(1), stop(2), strand(3), desc(4), name(5), alias(6), "N"(7), translation(8)

	seq_name = infile_name.split(".")[0]  # Extract filename without extension for SEQ_NAME

	for dp in finaldata:
		try:
			record_dp(dp, out, seq_name, wantproteinYN)
		except IndexError:
			try:
				record_dp(dp[0], out, seq_name, wantproteinYN)
			except IndexError:
				print("Index Error:")
				print(dp)
				print()
				continue

	return None

########################################################################################################################
def record_dp(datapt, outfilename, seq_name, proteinYN):
	'''Writes data from a list to the outfile'''

	# [start, end, strand, description, systematic_name, gene_name, pseudo_gene, alias, translation]
	outfilename.write(seq_name + "\t" 
	+ str(datapt[0]) + "\t"
	+ str(datapt[1]) + "\t"
	+ str(datapt[2]) + "\t"
	+ str(datapt[4]) + "\t"
	+ str(datapt[5]) + "\t"
	+ str(datapt[7]) + "\t"
	+ str(datapt[3]) + "\t"
	+ str(datapt[0]) + "\t"
	+ str(datapt[1]) + "\t"
	+ str(datapt[0]) + "\t"
	+ str(datapt[1]) + "\t"
	+ "1\t"	
	+ str(datapt[0]) + "\t"
	+ str(datapt[1]) + "\t"
	+ "y")

	if proteinYN == "Y":
		outfilename.write("\t" + datapt[8] + "\n")
	else:
		outfilename.write("\n")

	return None

########################################################################################################################
# Main Execution
########################################################################################################################

all_files = (os.listdir())
all_gbk_files = []

wantfasta = input("Output fasta file as well? (Y/N): ").strip().upper()  # Normalize input to uppercase
wantprotein = input("Output protein sequence in the MochiView file? (Y/N): ").strip().upper()

for i in all_files:
	if i.endswith(".gbff") or i.endswith(".refseq"):
		all_gbk_files.append(i)

# For every genbank file, convert to MochiView format
for f in all_gbk_files:
	print("Processing File: ", f)
	GenomeName, all_out_data = collect_all_info(f, wantfasta)
	polished_data = remove_duplicate_lines(all_out_data)
	create_data_lines(polished_data, f, GenomeName, wantprotein)

print("Work complete!")
print()
