import os,sys

a = open(sys.argv[1])
b = open(sys.argv[2],"w")
i = int(sys.argv[3])

tmp_dic = {}
for each in a:
	each_list=each.strip().split()
	gene=each_list[i]
	# score=each_list[6]
	if not tmp_dic.__contains__(gene):
		# b.write(gene+"\t"+score+"\n")
		b.write(each)
	# tmp_dic[gene]=score
	tmp_dic[gene]=1

a.close()
b.close()