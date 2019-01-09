#! usr/bin/python


import os, re
import requests
import vxm_hla 
import itertools
from vxm_hla import allele_truncate, locus_string_geno_list, expand_ac, single_locus_allele_codes_genotype, gl_string_alleles_list, allele_code_to_allele_list, genotype_allele_ag_freq

import conversion_functions_for_VXM
from conversion_functions_for_VXM import  gl_string_ags, genotype_ags, allele_code_ags, unosagslist, convert_allele_list_to_ags, allele_freq, map_ag_to_bw46, agbw46

import reverse_conversion

from reverse_conversion import map_single_ag_to_alleles

UA_eq_dict = {}



UNOS_UA_eq_filename = "UNOS_UA_ag_equivalencies.csv"
UNOS_UA_eq_file = open(UNOS_UA_eq_filename, 'r')

for row in UNOS_UA_eq_file:
	if row.startswith("Antigen"):
		continue
	else:
		row = row.strip("\n")
		row_split = row.split(",")
		ua_ag = row_split[0]
		ua_ag_eqs = row_split[1:]
		ua_ag_eqs = list(filter(None, ua_ag_eqs))
		UA_eq_dict[ua_ag] = ua_ag_eqs
#print(UA_eq_dict)


def vxm_uags(donorags, candidateags):
	conflicts = []
	donor_ags_alleles = []
	UA_list = []
	donor_bws_list = []
	##### maps candidate's UA to OPTN equivalents ########

	for ag in candidateags:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	

	
	recepient_ags = [item for sublist in UA_list for item in sublist]

	for ag in donorags:
		if ag in agbw46.keys():
			donorags.append(agbw46[ag])

	donorags = list(set(donorags))

	for ag in donorags:
		bw = map_ag_to_bw46(ag)
		donor_bws_list.append(bw)
		
	for ag in donorags:
		alleles = map_single_ag_to_alleles(ag)
		donor_ags_alleles.append([ag])
		if alleles:
			donor_ags_alleles.append(alleles)
	
	#print(donor_ags_alleles)
	merged_dags_alleles = list(itertools.chain(*donor_ags_alleles))

	for ag in recepient_ags:
		if ag in merged_dags_alleles:
			conflicts.append(ag)

	return (donorags, recepient_ags, conflicts)


def vxm_hIresalleles(donorsAlleleList, candidateags):
	conflicts = []
	donorags = []

	UA_list = []
	

	for ag in candidateags:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = [item for sublist in UA_list for item in sublist]

	donorags = convert_allele_list_to_ags(donorsAlleleList)

	donorags_alleles = donorsAlleleList + donorags
	
	for ag in recepient_ags:
		if ag in donorags_alleles:
			conflicts.append(ag)

	return(donorags, recepient_ags, conflicts)


def vxm_gls(donor_gl_string, donor_ethnicity, recipient_UA_list):
	conflicts = []
	
	new_ag_probs = {}
	bw_prob = {}
	donor_ags = []
	output = gl_string_ags(donor_gl_string, donor_ethnicity)
	ag_output = output[0]
	ag_probs = genotype_allele_ag_freq(ag_output)
	
	for i,j in ag_probs.items():
		je = j
		if je > 1:
			je = 1
		ag_probs[i] = je	
		if i == "Bw4" or i == "Bw6":
			bw_prob[i] = je
	#print(bw_prob)

	#donor_alleles = vxm_hla.gl_string_alleles_list(donor_gl_string)
	allele_output = output[1]
	
	donor_allele_freqs = genotype_allele_ag_freq(allele_output)
	
	donor_alleles = vxm_hla.gl_string_alleles_list(donor_gl_string)
	#donor_allele_freqs = conversion_functions_for_VXM.allele_freq(donor_alleles, donor_ethnicity)
	

	for k in ag_probs.keys():
		donor_ags.append(k)

	UA_list = []
	for ag in recipient_UA_list:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = list(filter(None, [item for sublist in UA_list for item in sublist]))


	donor_alleles_ags = donor_ags + donor_alleles
	
	for ag in recepient_ags:
		if ag in donor_alleles_ags:
			conflicts.append(ag)
	
	conflicts = list(filter(None, conflicts))	
	conflict_ag_probs = {}

	for i in conflicts:
		if i in ag_probs.keys():
			conflict_ag_probs[i] = ag_probs[i]

		elif i in donor_allele_freqs.keys():
			conflict_ag_probs[i] = donor_allele_freqs[i]

		else:
			conflict_ag_probs[i] = 0	
	
	for i,j in conflict_ag_probs.items():
		if j > 1.00:
			j = 1.00
			conflict_ag_probs[i] = j


	return(donor_ags, recepient_ags, conflicts, conflict_ag_probs, donor_allele_freqs, ag_probs, bw_prob)








def vxm_allele_codes(allele_codes_list, donor_ethnicity, recepient_UA_list):
	conflicts = []
	ag_probs = {}
	new_ag_probs = {}
	bw_prob = {}
	donor_ags = []
	output = allele_code_ags(allele_codes_list, donor_ethnicity)

	
	ag_output = output[0]
	ag_probs = genotype_allele_ag_freq(ag_output)
	#print(ag_probs)
	

	'''for ag in ag_probs.keys():
		if ag in UA_eq_dict.keys():
			p=0
			for eqs in UA_eq_dict[ag]:
				if eqs in ag_probs.keys():
					p += ag_probs[eqs]
					ag_probs[ag] = p
		else:
			continue'''

	
	for i,j in ag_probs.items():
		je = j
		if je > 1:
			je = 1
		ag_probs[i] = je
		if i == "Bw4" or i == "Bw6":
			bw_prob[i] = je	
	allele_output = output[1]
	donor_alleles = vxm_hla.allele_code_to_allele_list(allele_codes_list)
	#donor_allele_freqs = conversion_functions_for_VXM.allele_freq(donor_alleles, donor_ethnicity)
	donor_allele_freqs = genotype_allele_ag_freq(allele_output)
	#print(donor_allele_freqs)


	for k in ag_probs.keys():
		donor_ags.append(k)

	UA_list = []
	for ag in recepient_UA_list:
		if ag in UA_eq_dict.keys():
			UA_list.append(UA_eq_dict[ag])
		else:
			UA_list.append([ag])	


	recepient_ags = list(filter(None, [item for sublist in UA_list for item in sublist]))

	donor_alleles_ags = donor_ags + donor_alleles

	for ag in recepient_ags:
		if ag in donor_alleles_ags:
			conflicts.append(ag)
	
	conflict_ag_probs = {}

	
	for i in conflicts:
		if i in ag_probs.keys():
			conflict_ag_probs[i] = ag_probs[i]

		elif i in donor_allele_freqs.keys():
			conflict_ag_probs[i] = donor_allele_freqs[i]

		else:
			conflict_ag_probs[i] = 0



	for i,j in conflict_ag_probs.items():
		if j > 1.00:
			j = 1.00
			conflict_ag_probs[i] = j

	#print(conflict_ag_probs)


	return(donor_ags, recepient_ags, conflicts, conflict_ag_probs, donor_allele_freqs, ag_probs, bw_prob)





