import re
upstream_file = "/project/functional-genomics/2019/group3/peak_comparison/upstream2000.fa"
output_file = "/project/functional-genomics/2019/group3/peak_comparison/upstream2000.bed"
with open(upstream_file,"r") as in_file, open(output_file, "w") as out_file:
	for line in in_file:
		if line.startswith('>'):
			temp = line.split()[1]
			temp_list = re.split('[:-]', temp)
			out_file.write('\t'.join(temp_list) + '\n')
