#! usr/bin/python

import os, re

import pandas as pd

donor_typing = ['A2', 'A26', 'C05', 'C01', 'B72', 'B27', 'DR12', 'DR1', 'DQ7', 'DQ5', "DR51", "DR52"]


recepient_ags = ['A2', 'B15', 'B27', 'Bw4', 'DR2']

conflicting_ags = ['A2', 'A0203', 'B27', 'Bw4']

optn_equis = ['A2', 'A0203', 'B15', 'B62', 'B63', 'B75', 'B76', 'B77', 'B27', 'Bw4', 'DR2', 'DR15', 'DR16']

donor_alleles = []

output_dict = {}




###### get list of locus from donor typing

locus_list = []
for ag in donor_typing:
	if ag[0:2].isalpha():
		locus = ag[0:2]
		locus_list.append(locus)
	else:
		locus = ag[0]
		locus_list.append(locus)	


final_locus_list = sorted(list(set(locus_list)))		

#print(final_locus_list)



###### Group all the data together by locus 
###### Group ags in donor typing as per locus


########################################### donor typing ##########################################
donor_aags = []
donor_bags = []
donor_cags = []
donor_drags = []
donor_dqags = []

for ag in donor_typing:
	if re.findall("A", ag):
		donor_aags.append(ag)

	if re.findall("B", ag):
		donor_bags.append(ag)

	if re.findall("C", ag):
		donor_cags.append(ag)

	if re.findall("DQ", ag):
		donor_dqags.append(ag)		


	if re.findall("DR", ag):
		donor_drags.append(ag)		


#print(donor_aags) 
#print(donor_bags)
#print(donor_cags) 
#print(donor_drags) 
#print(donor_dqags)


############################################ candidate's Unacceptable Antigens #############################################

cand_auags = []
cand_buags = []
cand_cuags = []
cand_druags = []
cand_dquags = []

for uag in recepient_ags:
	if re.findall("A", uag):
		cand_auags.append(uag)

	if re.findall("B", uag):
		cand_buags.append(uag)

	if re.findall("C", uag):
		cand_cuags.append(uag)

	if re.findall("DQ", uag):
		cand_dquags.append(uag)		


	if re.findall("DR", uag):
		cand_druags.append(uag)		


#print("AUA: " + str(cand_auags))
#print("BUA: " + str(cand_buags))
#print("CUA: " + str(cand_cuags)) 
#print("DRUA: " + str(cand_druags))
#print("DQUA: " + str(cand_dquags))

################################################################# Conflicting Antigens ###############################################################

conflicts_aags = []
conflicts_bags = []
conflicts_cags = []
conflicts_drags = []
conflicts_dqags = []

for cag in conflicting_ags:
	if re.findall("A", cag):
		conflicts_aags.append(cag)

	if re.findall("B", cag):
		conflicts_bags.append(cag)

	if re.findall("C", cag):
		conflicts_cags.append(cag)

	if re.findall("DQ", cag):
		conflicts_dqags.append(cag)		


	if re.findall("DR", cag):
		conflicts_drags.append(cag)		


##############################################################OPTN equivalents ###########################################################################

optne_aags = []
optne_bags = []
optne_cags = []
optne_drags = []
optne_dqags = []

for optne in optn_equis:
	if re.findall("A", optne):
		optne_aags.append(optne)

	if re.findall("B", optne):
		optne_bags.append(optne)

	if re.findall("C", optne):
		optne_cags.append(optne)

	if re.findall("DQ", optne):
		optne_dqags.append(optne)		


	if re.findall("DR", optne):
		optne_drags.append(optne)






for loc in final_locus_list:
	if loc == "A":
		output_dict[loc] = (", ".join(donor_aags)),  (", ".join(cand_auags)), (", ".join(optne_aags)), (", ".join(conflicts_aags))

	if loc == "B":
		output_dict[loc] = (", ".join(donor_bags)),  (", ".join(cand_buags)), (", ".join(optne_bags)), (", ".join(conflicts_bags)), 

	if loc == "C":
		output_dict[loc] = (", ".join(donor_cags)), (", ".join(cand_cuags)),  (", ".join(optne_cags)), (", ".join(conflicts_cags)),

	if loc == "DR":
		output_dict[loc] = (", ".join(donor_drags)), (", ".join(cand_druags)),  (", ".join(optne_drags)), (", ".join(conflicts_drags)),

	if loc == "DQ":
		output_dict[loc] = (", ".join(donor_dqags)), (", ".join(cand_dquags)),	(", ".join(optne_dqags)), (", ".join(conflicts_dqags)), 			


#print(output_dict)


df = pd.DataFrame.from_dict(output_dict, orient='index')
df.columns = ["Donor Typing", "Candidate's Unacceptable Antigens", "OPTN Equivalents", "Conflicting Antigens"]
df.index.name = "Locus"
	
#print(df)

def conflict_cell_red(X):
	


df_html = pd.DataFrame.to_html(df)

print(df_html)

###############################################################################################################################################################
#ag_group = dict(zip(locus_list, donor_typing))


#print(ag_group)

##### make a dataframe

