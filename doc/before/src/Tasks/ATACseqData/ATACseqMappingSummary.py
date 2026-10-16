import os 
import sys

def mappingSummary(summary_info,output_info,output_name,ChromosomeDistribution_info):
	info = summary_info.readlines()
	TotalPairs = info[1].split("(")[0].strip()
	UnmappedPairs = info[2].split(" (")[0].strip()
	UniqMappedPairs = info[3].split(" (")[0].strip()
	MultiMappedPairs = info[4].split(" (")[0].strip()
	DiscordantlyAlignedPairs = info[7].split(" (")[0].strip()
	UnpairedReads = info[10].split("mates")[0].strip()
	UnmappedReads = info[11].split(" (")[0].strip()
	UniqMappedReads = info[12].split(" (")[0].strip()
	MultiMappedReads = info[13].split(" (")[0].strip()
	UniqAlignmentRate = float(UniqMappedPairs)/float(TotalPairs)*100
	OverallAlignmentRate = info[14].split("overall")[0].strip()
	MreadsRatio = (chrMReads(ChromosomeDistribution_info)/(float(TotalPairs)*2))*100
	
	output_info.write("{}\tpair\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}%\t{:.2f}%\n".format(output_name,TotalPairs,UnmappedPairs,UniqMappedPairs,MultiMappedPairs,DiscordantlyAlignedPairs,UnpairedReads,UnmappedReads,UniqMappedReads,MultiMappedReads,OverallAlignmentRate,UniqAlignmentRate,MreadsRatio))

def chrMReads(ChromosomeDistribution_info):
	M_reads = 0 
	for each in ChromosomeDistribution_info:
		each = each.strip().split("\t")
		if each[0] == "chrM":
			M_reads = int(each[1])
	return(M_reads)


def main():
	head_need = not os.path.exists(sys.argv[3])
	summary_info = open(sys.argv[1])
	ChromosomeDistribution_info = open(sys.argv[2])
	output_info = open(sys.argv[3],"a")
	output_name = sys.argv[4]
	if head_need:
		output_info.write("SampleName\tLayout\tTotalPairs\tUnmappedPairs\tUniqMappedPairs\tMultiMappedPairs\tDiscordantlyAlignedPairs\tUnpairedReads\tUnmappedReads\tUniqMappedReads\tMultiMappedReads\tOverallAlignmentRate\tUniqAlignmentRate\tMreadsRatio\n")
	mappingSummary(summary_info,output_info,output_name,ChromosomeDistribution_info)

	summary_info.close()
	ChromosomeDistribution_info.close()
	output_info.close()
main()