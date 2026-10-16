import os,sys

genome_file = sys.argv[1]
gene_file = sys.argv[2]
out_file = sys.argv[3]

genome_dic = {}
genome_info = open(genome_file)
for each in genome_info:
	each_list = each.strip().split()
	genome_dic[each_list[4]] = each
genome_info.close()

gene_info = open(gene_file)
out_info = open(out_file,"w")
for each_gene in gene_info:
	each_gene = each_gene.strip().split()[0]
	if genome_dic.has_key(each_gene):
		out_info.write(genome_dic[each_gene])
gene_info.close()
out_info.close()
