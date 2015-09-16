_author__ = 'Evanthia Kaimaklioti'
#requirements are sqlite3 
# there are several print statements, the most important one the print json


import sqlite3
import re
import itertools
import collections
from operator import itemgetter
import operator
import argparse
import json
import json as simplejson
from collections import defaultdict
from itertools import izip


parser = argparse.ArgumentParser()
parser.add_argument('inputfile', nargs = '?', default=argparse.SUPPRESS)
parser.add_argument('-in', dest='inputfile', default = 'NO FILE SUPPLIED')          # e.g. -in testing.blastx
	
parser.add_argument('blast_score', nargs = '?', default=argparse.SUPPRESS)
parser.add_argument('-b', dest='blast_score', default = 50)                         # -b 75 

parser.add_argument('outputfile', nargs = '?', action = 'store', default=argparse.SUPPRESS)
parser.add_argument('-out', dest='outputfile', default = 'NO OUTPUT FILE NAME SUPPLIED')  # -out *output_filename.txt

parser.add_argument('blast_threshold', nargs = '?',  default=argparse.SUPPRESS)
parser.add_argument('-bt', dest='blast_threshold', default = 0.9)    # -bt 0.8 

filename = parser.parse_args().inputfile
blast_score_arg = parser.parse_args().blast_score
out_filename = parser.parse_args().outputfile
blastThreshold = parser.parse_args().blast_threshold

def main():	

	taxonomy_database = make_db()
	##all_taxids_defaultdict, all_hits = parse_blast_file_open()    ## if you wish to call this function instead see instructions below
	all_taxids_defaultdict = parse_blast_file_open1()
	print "Contig and the list of taxids of the successful hits", all_taxids_defaultdict
	##print "the list of successful hits for this contig is", all_hits   
	
	
	##count_hits = ""
	##for contig1, element in all_hits.iteritems():
	##	count_hits = len(element)
	##	print "hits count for this contig is", contig1, count_hits
	
	
	overall_data = collections.defaultdict(list)
	LCA = ""    #declare the LCA variable
	for contig in all_taxids_defaultdict:
		my_contig_taxonomies = [] 
		 
		selected_taxids = all_taxids_defaultdict[contig]  #captures the taxids of the selected hits
		for taxid in selected_taxids:
			if taxid == 0:   # i.e. this contig, did not have significant hits (above the defined minimum blast score)
				LCA = 0
			else:
				taxonomy1 = create_taxonomy_list3(taxid)
				my_contig_taxonomies.append(taxonomy1)
				LCA = find_LCA(my_contig_taxonomies)
		overall_data[LCA].append(contig)
	print "overall_data", overall_data

	LCA_list = list(overall_data.keys())
	print "the LCA list is", LCA_list

	
	
	
	# write the information to an excel file with the columns corresponding to the following in the order from left to right
	# LCA's taxid, organisms name, all of the contigs that correspond to that LCA, sum (overall number of contigs) that are under that LCA
	with open (out_filename, 'w') as file:
		update_taxonomyDB_OK = update_taxonomyDb()
		commonKeys=set(overall_data.keys()) & set(update_taxonomyDB_OK.keys())  #change this to the variable, name_on_fly_dict
		list_keys = list(commonKeys)
		
		for key in list_keys: 
			length_of_list = len(overall_data[key])

			file.writelines(str(key))
			file.writelines('\t')
			file.writelines(str(update_taxonomyDB_OK.get(key)))
			file.writelines('\t')
			file.writelines(str(overall_data.get(key)))
			file.writelines('\t')
			file.writelines(str(length_of_list))
			file.writelines('\n')
			
	get_everything_dict = create_lineage(LCA_list)
	names_dict = create_names_dict(get_everything_dict)

	names_lineage_dict_prior_json = create_structure_before_json(names_dict, 'root')
	print "names_lineage_dict_prior_json",  names_lineage_dict_prior_json
	print (json.dumps(names_lineage_dict_prior_json, indent = 1))
	#eventually you can code for this structure to be printed in a separate file, as *filename.json to be quickly 
	# visualising the dendrogram 

# input: blast tabulated format file 
# function: finds the highest blast bit score (BBS) for each contig, sets the BBS threshold window, selects the successful hits
# for the contig 
# note: selects on the first TAXID from the hits that have more than one taxid


def parse_blast_file_open1():
	highest_score = 0
	threshold = 0
	info = []
	blast_score = 0
	taxid = ""
	data_dict = collections.defaultdict(list)
	with open (filename) as blast_file:
		for line in blast_file:
			for field in line: 
				field = line.strip().split("\t")
				contig = field[0]
				blast_score = field[11].replace(" ", "")
				taxid = field[12]
				# remove the ";" and select only the first taxid
				new_taxid = taxid.split(";") # new_taxid[0] is the first taxid
			info.append((contig, float(blast_score), new_taxid[0]))     # is the list of lists 
			sorted_info = sorted(info, key=itemgetter(1), reverse = True)
			d = collections.defaultdict(list)
			for k, v, c in sorted_info: 
				d[k].append(v)
				d[k].append(c)
		for key, value in d.items():         # value1 is ('Contig104', ['185', '4072', '185', '4072', '185', '4113', '183'])
			i = 0 
			blast_score = value[0::2]    #get every other item, starting with the first
			highest_score = float(blast_score[0])
			threshold = float(blastThreshold) * float(highest_score)
		
			for items in blast_score:
				taxid = value[1::2]
				if highest_score < float(blast_score_arg):
					data_dict[key].append(0)
				elif items <= threshold:
					pass
				else:
					data_dict[key].append(int(taxid[i]))
				i += 1
	blast_file.close()
	#print "data_dict", data_dict
	return data_dict
	
	
	
	
# input is the blast output tab-delimited file
# also parses only the blast hits that are within the user defined threshold for the highest blast score
# selects only the first taxid for any hits that have more than one taxid, i.e. as MEGAN does. 
# output: the defaultdict of lists (data_dict): key is the Contig name and value is the list of the taxids of the successful hits
# output: defaultdict of lists (hits_dict): key is contig name and the value list of the subject name of the successful hits 
# if you select to work with this function you need to uncomment out lines the lines with the ## in front
# i.e. lines : 	40, 43, 46-49. and comment out line: 41
def parse_blast_file_open():
	highest_score = 0
	threshold = 0
	info = []
	blast_score = 0
	taxid = ""
	data_dict = collections.defaultdict(list)
	hits_dict = collections.defaultdict(list)
	hits_list = []

	for line in open(filename):
		if ';' in line:
			pass
		else:
			for field in line: 
				field = line.strip().split("\t")
				contig = field[0]
				hit_info = field[1]
				blast_score = field[11].replace(" ", "") #change here depending where the Blast bit score is. 
				taxid = field[12]   #change here depending where the TAXID is. 
			info.append((contig, float(blast_score), int(taxid), hit_info))   
			# sort the list according to the blast_score in descending order
			sorted_info = sorted(info, key=operator.itemgetter(1), reverse = True)
			d = collections.defaultdict(list)
			for k, v, c, d1 in sorted_info: 
				d[k].append(v)
				d[k].append(c)
				d[k].append(d1)
	for key, value in d.items(): # value1 is ('Contig104', ['185', '4072', '185', '4072', '185', '4113', '183'])
		i = 0 
		blast_score = value[0::3]    #get every 3 items, starting with the first
		highest_score = float(blast_score[0])
		print "the highest blast bit for this contig is", key, highest_score
		threshold = float(blastThreshold) * float(highest_score)
		print "the threshold for this contig is,", key, threshold
		contigs_list = value[2::3]
		for items in blast_score:
			taxid = value[1::3]
			if highest_score < float(blast_score_arg):
				data_dict[key].append(0)
				hits_dict[key].append(0)
			elif items <= threshold:
				pass
			else:
				hits_dict[key].append(contigs_list[i])
				data_dict[key].append(int(taxid[i]))		
			i += 1
	return (data_dict, hits_dict)




# its input is the taxonomies to be tested
# outputs the LCA (the first common taxid in all taxonomies)
# eventually this function will have as input the selected taxids and wiill go to the Db to bring the taxonomy
# exception handling: if a taxid is not in the Db then it should be excluded from the selected hits. (Make a note of this in the documentation)
# input dictionary holding the taxonomy for each contig
# output: the LCA out of all taxonomies
def find_LCA(taxonomy):
    for taxid in taxonomy[0]:   
	if all([taxid in taxidB for taxidB in taxonomy[1:]]):
		return taxid      # return the first taxid common
    return False             # if there is no taxid in common, don't return anything



# input: nodes.dmp & names.dmp
# output: a four column SQL database
# this function will be utilised to produce the taxonomy for the records for which the LCA will be calculated
# output: builds the four column table that will contain the taxid, parent, rank and name of the NCBI hit record
def make_db():
	con = sqlite3.connect('taxonomy_db.db')
	cur = con.cursor()
	cur.execute('DROP TABLE IF EXISTS Taxonomy')      # if the Taxonomy table exists it drops it, to avoid the shell error 
	cur.execute('CREATE TABLE Taxonomy(taxid INT, parentID INT, rank TEXT, name TEXT, PRIMARY KEY (taxid))')
	data = []   #build the data as a list so that the executemany() and the commit() won't have to be executed everytime there is a new record
	for line in open('nodes.dmp', 'rU'):
		fields = re.split('\t\|\t', line)
		taxid = fields[0]
		parent = fields[1]
		rank = fields[2]
		data.append((taxid, parent, rank, 'NULL'))
	cur.executemany("INSERT INTO Taxonomy(taxid, parentID, rank, name) VALUES (?, ?, ?, ?);", data)
	con.commit()
	return 'done'
	
	
# updates the Taxonomy table by putting the scientific name in the NULL column.
# output: names_dict is build as to be used when returning the csv file instead of parsing the MySQL db
def update_taxonomyDb():
	con = sqlite3.connect('taxonomy_db.db')
	print "Opened Database successfully"
	cur = con.cursor()
	#cur.execute('DROP TABLE IF EXISTS Taxonomy')      # if the Taxonomy table exists it drops it, to avoid the shell error 
	scientific_name = ""
	names_fly_dict = {}
	names_fly_dict[0] = 'Not assigned'
	taxid_names = ""
	for line in open('names.dmp', 'r'):
		field = line.strip().split("\t")
		taxid_names = field[0]  #captures all taxids, irrespective if they correspond to 'scientific name' or not
		notion = field[6]
		if notion == 'scientific name':
			scientific_name = field[2]      #captures the scientific name
		#error catching : need to strip the scientific name from the apostrophies
		cur.execute('UPDATE Taxonomy SET name = ? WHERE taxid=?',(scientific_name, taxid_names)) # it changed the taxid '2' with Hedyotis-Oldenlandia complex
		names_fly_dict[int(taxid_names)] = scientific_name
	con.commit()
	con.close()
	return names_fly_dict

	

# input a taxid 
# this function goes through the taxonomy_db.db and searches with the taxid PRIMARY KEY to bring the taxid append it
# into a taxonomy list and then continue appending the parent id
# then this parent id will be the taxid.
# querying a database to select rows and retrieving row entries
# output: brings the taxonomy lineage until the 'root' for each taxid 
def create_taxonomy_list3(tax_id):   
	conn = sqlite3.connect('taxonomy_db.db')
	c = conn.cursor()
	
	taxonomy = []
	taxonomy.append(tax_id)
	while tax_id !=1:
		c.execute("SELECT (parentID) FROM Taxonomy WHERE taxid=:some_id",{"some_id": tax_id})
		all_rows = c.fetchone()
		if all_rows:
			parentID = all_rows[0]
			taxonomy.append(parentID)
			tax_id = parentID
		else:
			break
	#c.close()
	return taxonomy


#input: LCA_list
#function & output: creates a defaultdict of lists structure that holds the lineage relationship of the taxids in the LCA_list
def create_lineage(lca_list):  
	conn = sqlite3.connect('taxonomy_db.db')
	c = conn.cursor()
	lineage_dict = collections.defaultdict(list)
	for tax_id in lca_list:
		taxonomy = []
		taxonomy.append(tax_id)
		# put exception loop: 
			#if the tax_id is 0 or 1 then append Root. AND THUS IN THE DICTIONARY PUT .... defaultdict['Root'].append((1,0))
		if tax_id == 0:    # 0 is the 'Not assigned'
			lineage_dict['root'].append(tax_id)
		else:
			while tax_id != 1:
				c.execute("SELECT (parentID) FROM Taxonomy WHERE taxid=:some_id",{"some_id": tax_id})
				all_rows = c.fetchone()
				if all_rows:
					parentID = all_rows[0]
					taxonomy.append(parentID)
					if tax_id in lineage_dict[parentID]:  # if the taxid is a child node pass
						pass
					else:
						lineage_dict[parentID].append(tax_id)  # else put as key the parent of that taxid and make that taxid as the value (i.e. child)
					tax_id = parentID
				else:
					break
		#c.close()
	return lineage_dict

	
# function & output: takes the lineage_dict and returns a defaultdict of list with the same relationships in the nodes
# with their scientific names
# it's best if this is eventually merged with the create_lineage() function	
def create_names_dict(d):
	conn = sqlite3.connect('taxonomy_db.db')
	c = conn.cursor()
	new_name_dict = collections.defaultdict(list)
	for key, values in d.iteritems():
		c.execute("SELECT (name) FROM Taxonomy WHERE taxid=:some_id",{"some_id": key})
		all_rows = c.fetchone()
		if all_rows:
			name = str(all_rows[0])	
		for element in values: 
			if element == 0:
				new_name_dict['root'].append('Not assigned')
			else:
				c.execute("SELECT (name) FROM Taxonomy WHERE taxid=:some_id",{"some_id": element})
				all_rows1 = c.fetchone()
				if all_rows1:
					name1 = str(all_rows1[0])
					new_name_dict[name].append(name1)
	print "new name dict", new_name_dict
	return new_name_dict


# recursion function to build the dictionary in the form as the flare.json, currently without the value variable
def create_structure_before_json(d, node):
	new_node = str(node)
	return(
		{'name':new_node, 'children':[create_structure_before_json(d, child_node) for child_node in d[new_node]]}
	)
	
	
if __name__ == "__main__":
	main()





