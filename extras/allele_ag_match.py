   #! usr/bin/python
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



#print(allele_to_ag_dict)

ags = {'A2': 1.0, 'A26': 1.0, 'C05': 1.0, 'C01': 1.0, 'B75': 0.5723, 'B27': 1.0, 
		'Bw6': 1.0, 'Bw4': 1.0, 'B62': 0.423, 'B72': 0.0047, 'DR12': 1.0, 'DR1': 1.0, 'DQ7': 1.0, 'DQ5': 1.0}

alleles = {'A*02:01': 0.0844, 'A*26:01': 0.0351, 'A*02:55': 0, 
			'A*26:07': 0.0001, 'C*05:01': 0.0079, 'C*01:02': 0.11, 'B*15:01': 0.0312, 
			'B*15:02': 0.0422, 'B*15:03': 0.0003, 'B*15:04': 0.0, 'B*27:05': 0.0082, 
			'DRB1*12:01': 0.0204, 'DRB1*01:01': 0.0264, 'DQB1*03:01': 0.1849, 'DQB1*05:01': 0.0809}


conflicts = {'A2': 1.0, 'Bw4': 1.0, 'A*02:01': 0.2716, 'B15': 1, 'DR1':0.1}



list_of_antigens_per_locus = [['A2: 1.0', 'A26: 1.0'], ['B72: 0.8535', 'B27: 0.9999', 'Bw6: 0.9999', 'Bw4: 0.9999', 
                              'B62: 0.1414', 'B75: 0.005'], ['C05: 1.0', 'C01: 1.0'], ['DR12: 1.0', 'DR1: 1.0'], ['DQ7: 1.0', 'DQ5: 1.0']]


list_of_alleles_per_locus = [['A*02:01: 0.1222', 'A*26:01: 0.0143', 'A*02:55: 0', 'A*26:07: 0.0'], ['B*15:01: 0.0104', 'B*15:02: 0.0004', 'B*15:03: 0.0631', 'B*15:04: 0.0', 'B*27:05: 0.0081'], 
							['C*05:01: 0.0326', 'C*01:02: 0.0078'], ['DRB1*12:01: 0.0373', 'DRB1*01:01: 0.026'], ['DQB1*03:01: 0.1851', 'DQB1*05:01: 0.1486']]




ag_list = []

allele_ag_dict = {}

for i, j in alleles.items():
	ag = allele_to_ag_dict[i][0]
	prob = ags[ag]

	ag_prob = ag + ": " + str(prob)
	ag_list.append(ag_prob)

#print(ag_list)	

a_ags = [] 
b_ags = []
c_ags = []
dr_ags = []
dq_ags = []


for i, j in alleles.items():
	ag = allele_to_ag_dict[i][0]
	prob = ags[ag]
	ag_prob = ag + ": " + str(prob)
	if ag_prob.startswith("A"):
		a_ags.append(ag_prob)

	if ag_prob.startswith("B"):
		b_ags.append(ag_prob)

	if ag_prob.startswith("C"):
		c_ags.append(ag_prob)

	if ag_prob.startswith("DR"):
		dr_ags.append(ag_prob)

	if ag_prob.startswith("DQ"):
		dq_ags.append(ag_prob)				


list_of_ag_probs = [a_ags] + [b_ags] + [c_ags] + [dr_ags] + [dq_ags] 

#print(list_of_ag_probs)






cag_prob_dict = {}

for i,j in alleles.items():
	ag = allele_to_ag_dict[i][0]
	if i in conflicts:
		prob = conflicts[i]
		cag_prob_dict[i] = i + ": " + str(prob)

	if ag in conflicts:
		prob = conflicts[ag]
		if i in cag_prob_dict.keys():
			cag_prob_dict[i] = cag_prob_dict[i] + " " + (ag + ": " + str(prob))

		else:
			cag_prob_dict[i] = 	ag + ": " + str(prob)

	else:
		cag_prob_dict[i] = ""


print(cag_prob_dict)

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
print(list_of_cag_probs) 