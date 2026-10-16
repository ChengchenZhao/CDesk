import os,sys
a = sys.argv[1]
b = open(sys.argv[2],"a")
b.write("sample_name\ttotal_tags\tCDS_Exons\t5'UTR_Exons\tt3'UTR_Exons\tIntrons\n")
for each_sample in a.split(",") :
	each_info = open(each_sample)
	sample_name = each_sample.split(".")[0]
	for each in each_info :
		if each.startswith("Total Tags"):
			total_tags = each.strip().split(" ")[-1]
		if each.startswith("CDS_Exons"):
			CDS_Exons = each[40:].split(" ")[0]
		if each.startswith("5'UTR_Exons"):
			FUTR_Exons = each[40:].split(" ")[0]
		if each.startswith("3'UTR_Exons"):
			TUTR_Exons = each[40:].split(" ")[0]
		if each.startswith("Introns"):
			Introns = each[40:].split(" ")[0]
	each_info.close()
	b.write("%s\t%s\t%s\t%s\t%s\t%s\n"%(sample_name,total_tags,CDS_Exons,FUTR_Exons,TUTR_Exons,Introns))

b.close()