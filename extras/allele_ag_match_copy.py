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




allele_dict = {'A*02:01': 0.0844, 'A*26:01': 0.0351, 'A*02:55': 0, 
			'A*26:07': 0.0001, 'C*05:01': 0.0079, 'C*01:02': 0.11, 'B*15:01': 0.0312, 
			'B*15:02': 0.0422, 'B*15:03': 0.0003, 'B*15:04': 0.0, 'B*27:05': 0.0082, 
			'DRB1*12:01': 0.0204, 'DRB1*01:01': 0.0264, 'DQB1*03:01': 0.1849, 'DQB1*05:01': 0.0809}

ag_list = ['A2', 'A26', 'Bw4', 'A*02:01', 'C05', 'B15', 'B*15:01', 'DRB1*12:01']


cag_prob_dict = {}

for i,j in allele_dict.items():
	ag = allele_to_ag_dict[i][0]
	if i in ag_list:
		cag_prob_dict[i] = i 

	if ag in ag_list:
		if i in cag_prob_dict.keys():
			cag_prob_dict[i] = cag_prob_dict[i] + " " + "(" + ag + ")"

		else:
			cag_prob_dict[i] = 	ag 

	else:
		cag_prob_dict[i] = ""
	
print(cag_prob_dict)
	



