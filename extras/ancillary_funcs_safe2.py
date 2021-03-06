de 
#! usr/bin/python

import os, re
import vxm_hla


allele_to_ag_dict = {}
UNOS_conversion_table_filename = "UNOS_conversion_table_with_rules.csv"
UNOS_conversion_table_file = open(UNOS_conversion_table_filename, 'r')
for row in UNOS_conversion_table_file:
	expression_character = ""
	if row.startswith("Allele"):
		continue 
	else:
		allele = row.split(',')[0]
		allele_4d = vxm_hla.allele_truncate(allele)
		antigen = row.split(',')[1]
		bw4_6 = row.split(',')[3]
	
	allele_to_ag_dict[allele_4d] = [antigen, bw4_6]




def group_serotypes_per_locus(ag_list):
	a_locus_ag = []
	b_locus_ag = []
	c_locus_ag = []
	dr_locus_ag = []
	dq_locus_ag = []

	for ag in ag_list:
		if re.findall("A", ag):
			a_locus_ag.append(ag)

		if re.findall("B", ag):
			b_locus_ag.append(ag)

		if re.findall("C", ag):
			c_locus_ag.append(ag)

		if re.findall("DR", ag):
			dr_locus_ag.append(ag)		
		
		if re.findall("DQ", ag):
			dq_locus_ag.append(ag)		

	locus_sorted_ag_list = [", ".join(sorted(a_locus_ag))] + [", ".join(sorted(b_locus_ag))] + [", ".join(sorted(c_locus_ag))] + [", ".join(sorted(dr_locus_ag))] + [", ".join(sorted(dq_locus_ag))]


	return locus_sorted_ag_list	




def split_gl_string_per_locus(gl_string):
	a_string = ""
	b_string = ""
	c_string = ""
	dr_string = ""
	dq_string = ""

	gl_string_split = gl_string.split("^")	

	for string in gl_string_split:
		if string.startswith("A"):
			a_string = string.replace("+", " + ") 
			a_string = a_string.replace("/", "/ ")

		if string.startswith("B"):
			b_string = string.replace("+", " + ") 
			b_string = b_string.replace("/", "/ ")

		if string.startswith("C"):
			c_string = string.replace("+", " + ")
			c_string = c_string.replace("/", "/ ") 

		if string.startswith("DR"):
			dr_string = string.replace("+", " + ") 
			dr_string = dr_string.replace("/", "/ ")

		if string.startswith("DQ"):
			dq_string = string.replace("+", " + ") 
			dq_string = dq_string.replace("/", "/ ") 



	string_list = [a_string] + [b_string] + [c_string] + [dr_string] + [dq_string]	
	
	return string_list				


##### for alleles #################

def prob_dict_list_of_strings(allele_probs):

	a_alleles = [] 
	b_alleles = []
	c_alleles = []
	dr_alleles = []
	dq_alleles = []


	for i,j in allele_probs.items():
		if i.startswith("A"):
			ji = str(j)
			fi = i + ": " + ji 
			a_alleles.append(fi)

		if i.startswith("B"):
			ji = str(j)
			fi = i + ": " + ji 
			b_alleles.append(fi)


		if i.startswith("C"):
			ji = str(j)
			fi = i + ": " + ji 
			c_alleles.append(fi)

		if i.startswith("DQ"):
			ji = str(j)
			fi = i + ": " + ji 
			dq_alleles.append(fi)       

		if i.startswith("DR"):
			ji = str(j)
			fi = i + ": " + ji 
			dr_alleles.append(fi)



	list_of_allele_probs = [a_alleles] + [b_alleles] + [c_alleles] + [dr_alleles] + [dq_alleles] 

	return list_of_allele_probs


###### for antigens, conflicting antigens, optne equivalents

def prob_dict_list_of_strings_for_antigens(alleles_dict, ag_probs):

	a_ags = [] 
	b_ags = []
	c_ags = []
	dr_ags = []
	dq_ags = []


	for i, j in alleles_dict.items():
		bw_prob = ""
		ag_prob = ""

		ag = allele_to_ag_dict[i][0]
		bw = allele_to_ag_dict[i][1]
		prob1 = ag_probs[ag]
		ag_prob = ag + ": " + str(prob1)
		
		if bw != "NA":
			prob2 = ag_probs[bw]
			bw_prob = bw + ":" + str(prob2)
		
		if ag_prob.startswith("A"):
			a_ags.append(ag_prob)

		if ag_prob.startswith("B"):
			if bw_prob:
				b_ags.append(ag_prob + ", " + bw_prob)
			else:
				b_ags.append(ag_prob)	


		if ag_prob.startswith("C"):
			c_ags.append(ag_prob)

		if ag_prob.startswith("DR"):
			dr_ags.append(ag_prob)

		if ag_prob.startswith("DQ"):
			dq_ags.append(ag_prob)				


	list_of_ag_probs = [a_ags] + [b_ags] + [c_ags] + [dr_ags] + [dq_ags] 

	return list_of_ag_probs

def group_allele_codes_or_list_of_alleles_per_locus(donor_typing_list):
	donor_a_alleles = []
	donor_b_alleles = []
	donor_c_alleles = []
	donor_dr_alleles = []
	donor_dq_alleles = []
	

	for i in donor_typing_list:
		locus = i.split("*")[0]
		if locus == "A":
			donor_a_alleles.append(i)

		if locus == "B":
			donor_b_alleles.append(i)

		if locus == "C":
			donor_c_alleles.append(i)

		if (locus == "DRB1") or (locus == "DRB3") or (locus == "DRB4") or (locus == "DRB5"):
			donor_dr_alleles.append(i)

		if (locus == "DQB1") or (locus == "DQA1"):
			donor_dq_alleles.append(i)


	final_typing_list = [", ".join(sorted(donor_a_alleles))] + [", ".join(sorted(donor_b_alleles))] + [", ".join(sorted(donor_c_alleles))] + [", ".join(sorted(donor_dr_alleles))] + [", ".join(sorted(donor_dq_alleles))] 

	return final_typing_list



def conflicts_ags(allele_dict, conflicts):

	cag_prob_dict = {}

	for i,j in allele_dict.items():
		ag = allele_to_ag_dict[i][0]
		bw = allele_to_ag_dict[i][1]
		if i in conflicts:
			prob = conflicts[i]
			if bw in conflicts:
				prob2 = conflicts[bw]
				cag_prob_dict[i] = i + ": " + str(prob) + " " + bw + ": " + str(prob2)
			
			else:	
				cag_prob_dict[i] = i + ": " + str(prob)

		if ag in conflicts:
			prob = conflicts[ag]
			if i in cag_prob_dict.keys():
				cag_prob_dict[i] = cag_prob_dict[i] + " " + (ag + ": " + str(prob))

			else:
				cag_prob_dict[i] = 	ag + ": " + str(prob)

		else:
			cag_prob_dict[i] = ""


#print(cag_prob_dict)

	a_cags = []
	b_cags = []
	c_cags = []
	dr_cags = []
	dq_cags = []

	for i,j in cag_prob_dict.items():

		if i.startswith("A"):
			a_cags.append(j)

		if i.startswith("B"):
			b_cags.append(j)

		if i.startswith("C"):
			c_cags.append(j)

		if i.startswith("DR"):
			dr_cags.append(j)

		if i.startswith("DQ"):
			dq_cags.append(j)			

	list_of_cag_probs = [a_cags] + [b_cags] + [c_cags] + [dr_cags] + [dq_cags]
#print(list_of_cag_probs) 

	return list_of_cag_probs
