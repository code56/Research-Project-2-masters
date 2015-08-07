__author__ = 'evanthia'
#comments : description of the software
#uses this version of biopython
#its requirements are sqlite3 

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

'''
Escherichia albertii B090 taxid: 550692   550692 208962 561 543 91347 1236 1224 2 131567

Streptomyces abikoensis taxid: 97398      97398 1883 2062 85011 1760 201174 2 131567

Streptomyces abietis taxid: 1227734       1227734 1883 2062 85011 1760 201174 2 131567

Trypanosoma vivax 	5699 47571 5690 5654 5653 33682 2759 131567 

Helicobacter brantae 	375927 209 72293 213849 29547 68525 1224 2 131567

Helicobacter cetorum 	138563 209 72293 213849 29547 68525 1224 2 131567 

Helicobacter macacae 	398626 209 72293 213849 29547 68525 1224 2 131567 

Helicobacter cholecystus 	45498 209 72293 213849 29547 68525 1224 2 131567

Olea ambrensis 	495749 4145 426106 4144 4143 91888 71274 1437201 91827 71240 1437183 3398 58024 78536 58023 3193 131221 35493 33090 2759 131567 

Daucus carota 	4039 4038 241799 241789 241778 4037 364270 4036 91882 71274 1437201 91827 71240 1437183 3398 58024 78536 58023 3193 131221 35493 33090 2759 131567 

'''

parser = argparse.ArgumentParser()
parser.add_argument('inputfile', nargs = '?', default=argparse.SUPPRESS)
parser.add_argument('-in', dest='inputfile', default = 'NO FILE SUPPLIED')
	
parser.add_argument('blast_score', nargs = '?', default=argparse.SUPPRESS)
parser.add_argument('-b', dest='blast_score', default = 50)

filename = parser.parse_args().inputfile
blast_score_arg = parser.parse_args().blast_score

LCA_list = [0, 1, 3398, 2759, 71240, 58024, 131567, 12176, 1437201, 12178, 91827, 131221, 3193, 1437183]

def main():	
	# Thrips. Correct LCA taxonomy ID = 45057
	taxonomy_test1 = {
		'670479' : ['670479','45057','153976','45053','45049','38130','30262','33342','33340','7496','85512','50557','6960','197562','197563','6656','88770','1206794','33317','33213','6072','33208','33154','2759','131567','1'],
		'161013' : ['161013','45057','153976','45053','45049','38130','30262','33342','33340','7496','85512','50557','6960','197562','197563','6656','88770','1206794','33317','33213','6072','33208','33154','2759','131567','1'],
		'1291313' : ['1291313','45057','153976','45053','45049','38130','30262','33342','33340','7496','85512','50557','6960','197562','197563','6656','88770','1206794','33317','33213','6072','33208','33154','2759','131567','1']
	}

	list_of_hits = {'550692': [550692,208962, 561, 543, 91347, 1236, 1224, 2, 131567],
                '97398' : [97398, 1883, 2062, 85011, 1760, 201174, 2, 131567],
               '1227734': [1227734, 1883, 2062, 85011, 1760, 201174, 2, 131567],
               '5699' : [5699, 47571, 5690, 5654, 5653, 33682, 2759, 131567],
				'375927' : [375927, 209, 72293, 213849, 29547, 68525, 1224, 2, 131567],
				'138563' : [138563, 209, 72293, 213849, 29547, 68525, 1224, 2, 131567],
				'398626' : [398626, 209, 72293, 213849, 29547, 68525, 1224, 2, 131567],
				 '45498' : [45498, 209, 72293, 213849, 29547, 68525, 1224, 2, 131567],
				'495749' : [495749, 4145, 426106, 4144, 4143, 91888, 71274, 1437201, 91827, 71240, 1437183, 3398, 58024, 78536, 58023, 3193, 131221, 35493, 33090, 2759, 131567],
				'4039' : [4039, 4038, 241799, 241789, 241778, 4037, 364270, 4036, 91882, 71274, 1437201, 91827, 71240, 1437183, 3398, 58024, 78536, 58023, 3193, 131221, 35493, 33090, 2759, 131567]
	}		

	#list = list_of_hits.values()   #this gives just the taxonomies
	#LCA = find_LCA(list)
	#test1 = find_LCA(taxonomy_test1.values())
	#print LCA
	#print test1,test2,test3,test4
	#taxonomy_database = make_db()
	#update_taxonomyDB = update_taxonomyDb()
	all_taxids_defaultdict = parse_blast_file_open()
	


	overall_data = collections.defaultdict(list)   
	
	for contig in all_taxids_defaultdict:
		my_contig_taxonomies = []  # or a dict
		#print "contig" + " " + str(contig) 
		selected_taxids = all_taxids_defaultdict[contig]
		#print "selected_taxids" + " " + str(selected_taxids)
		for taxid in selected_taxids:
			taxonomy1 = create_taxonomy_list3(taxid)
			my_contig_taxonomies.append(taxonomy1)
			#list = my_contig_taxonomies.values()   #this gives just the taxonomies
			LCA = find_LCA(my_contig_taxonomies)
			
		overall_data[LCA].append(contig)
			


# input is the blast output tab-delimited file
# also parses only the blast hits that are within the 90% of the highest blast score
# eventually the input blast file with be given from an argparse function
# output: the best hits will be written in an excel file
# the contig name will be the key
# all the lines in a list
def parse_blast_file_open():
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
				description = field[1:]
				blast_score = field[11].replace(" ", "")
				taxid = field[12]
				# remove the ";" and select only the first taxid
				new_taxid = taxid.split(";") # new_taxid[0] is the first taxid

			info.append((contig, blast_score, new_taxid[0]))     # is the list of lists 
			# sort the list according to the blast_score in descending order
			sorted_info = sorted(info, key=itemgetter(1), reverse = True)
			d = collections.defaultdict(list)
			for k, v, c in sorted_info: 
				d[k].append(v)
				d[k].append(c)
		for key, value in d.items():         # value1 is ('Contig104', ['185', '4072', '185', '4072', '185', '4113', '183'])
			i = 0 
			blast_score = value[0::2]    #get every other item, starting with the first
			highest_score = float(blast_score[0])
			threshold = 0.9 * float(highest_score)
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
	return data_dict
	
	
# finds the LCA out of all taxonomies
# its input is the taxonomies to be tested
# outputs the LCA (the first common taxid in all taxonomies)
# eventually this function will have as input the selected taxids and wiill go to the Db to bring the taxonomy
# exception handling: if a taxid is not in the Db then it should be excluded from the selected hits. (Make a note of this in the documentation)
# input dictionary holding the taxonomy for each contig

def find_LCA(taxonomy):
    for taxid in taxonomy[0]:   #taxonomy[0] is the first taxonomy in the list of taxonomies e.g. [ [97398, 1883, 2062, 85011, 1760, 201174, 2, 131567]]
	if all([taxid in taxidB for taxidB in taxonomy[1:]]):
		return taxid      # return the first taxid common
    return False             # if there is no taxid in common, don't return anything


# builds the four column table that will contain the taxid, parent, rank and name of the NCBI hit record
# input: nodes.dmp & names.dmp
# output: a four column SQL database
# this function will be utilised to produce the taxonomy for the records for which the LCA will be calculated
# there will be another function where the name column will be populated by searching the 

def make_db():
	con = sqlite3.connect('taxonomy_db.db')
	cur = con.cursor()
	cur.execute('DROP TABLE IF EXISTS Taxonomy') # if the Taxonomy table exists it drops it, to avoid the shell error 
	cur.execute('CREATE TABLE Taxonomy(taxid INT, parentID INT, rank TEXT, name TEXT, PRIMARY KEY (taxid))')
	
	data = []  
	for line in open('nodes.dmp', 'rU'):
		fields = re.split('\t\|\t', line)
		taxid = fields[0]
		parent = fields[1]
		rank = fields[2]
		data.append((taxid, parent, rank, 'NULL'))
	cur.executemany("INSERT INTO Taxonomy(taxid, parentID, rank, name) VALUES (?, ?, ?, ?);", data)
	con.commit()

	return 'done'


def update_taxonomyDb():
	con = sqlite3.connect('taxonomy_db.db')
	#print "Opened Database successfully"
	cur = con.cursor()
	scientific_name = ""
	dict_names = {}

	for line in open('names.dmp', 'r'):
		field = line.strip().split("\t")
		taxid_names = field[0]  #captures all taxids, irrespective if they correspond to 'scientific name' or not
 		notion = field[6]
		if notion == 'scientific name':
			scientific_name = field[2]      #captures only the scientific name from the file
		#error catching : need to strip the scientific name from the apostrophies
		# here I can also be producing the names dictionary. So that when I will be looking for the names I won't have to be going row by row in the SQL database. 
		
		cur.execute('UPDATE Taxonomy SET name = ? WHERE taxid=?',(scientific_name, taxid_names)) # it changed the taxid '2' with Hedyotis-Oldenlandia complex
		#dict_names[taxid_names] = scientific_name
	con.commit()
	con.close()
	print "I created the dictionary of names while updating the taxonomy database", dict_names
	return dict_names


# input a taxid 
# this function goes through the taxonomy_db.db and searches with the taxid PRIMARY KEY to bring the taxid append it
# into a taxonomy list and then continue appending the parent id
# then this parent id will be the taxid.
# querying a database to select rows and retrieving row entries
def create_taxonomy_list3(tax_id):    #have as input the taxid list coming from the blast parsing function then for each taxid (BUT MUST KNOW THE CONTIG it was assigned to). 
	conn = sqlite3.connect('taxonomy_db.db')
	c = conn.cursor()
	taxonomy = []
	taxonomy.append(tax_id)
	#print "taxonomy from create taxonomy list function", str(taxonomy)
	while tax_id !=1:
		c.execute("SELECT (parentID) FROM Taxonomy WHERE taxid=:some_id",{"some_id": tax_id})
		all_rows = c.fetchone()
		if all_rows:
			parentID = all_rows[0]
			taxonomy.append(parentID)
			tax_id = parentID
		else:
			#print "this taxid" + str(all_rows) + "doesn't exist"   #DON'T include this taxid in the LCA tree
			break
	#print "new taxonomy"
	#print taxonomy
	#c.close()
	return taxonomy
	
if __name__ == "__main__":
	main()





