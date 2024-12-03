import os,sys

input_info = open(sys.argv[1])
out_info = open(sys.argv[2],"w")

gene_dic = {}

for each in input_info:
	each_ll = each.strip().split("\t")
	if each_ll[0]=="chr":
		out_info.write("\t".join(each_ll[7:])+"\n")
	else:
		if not each_ll[7] == "-":
			if not gene_dic.__contains__(each_ll[7]):
				each_nll = [x if x !="" else "NA" for x in each_ll]
				each_nll = [x if x !="" else "0" for x in each_ll]
				out_info.write("\t".join(each_nll[7:])+"\n")
				gene_dic[each_ll[7]] = 1

input_info.close()
out_info.close()