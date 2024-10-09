import os 
import sys

def mappingSummary(summary_info,output_info,output_name,):

	info = summary_info.readlines()
	TotalPairs = info[1].split("(")[0].strip()
	UnmappedPairs = info[2].split(" (")[0].strip()
	UniqMappedPairs = info[3].split(" (")[0].strip()
	MultiMappedPairs = info[4].split(" (")[0].strip()
	try:
		DiscordantlyAlignedPairs = info[7].split(" (")[0].strip()
		UnpairedReads = info[10].split("mates")[0].strip()
		UnmappedReads = info[11].split(" (")[0].strip()
		UniqMappedReads = info[12].split(" (")[0].strip()
		MultiMappedReads = info[13].split(" (")[0].strip()
		UniqAlignmentRate = float(UniqMappedPairs)/float(TotalPairs)*100
		OverallAlignmentRate = info[14].split("overall")[0].strip()
		t = "pair"
	except:
		DiscordantlyAlignedPairs = ""
		UnpairedReads = ""
		UnmappedReads = ""
		UniqMappedReads = ""
		MultiMappedReads = ""
		UniqAlignmentRate = 0.0
		OverallAlignmentRate = info[5].split("overall")[0].strip()
		t = "single"
			
	output_info.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.2f}%\n".format(output_name, t,TotalPairs,UnmappedPairs,UniqMappedPairs,MultiMappedPairs,DiscordantlyAlignedPairs,UnpairedReads,UnmappedReads,UniqMappedReads,MultiMappedReads,OverallAlignmentRate,UniqAlignmentRate))



def main():
	head_need = not os.path.exists(sys.argv[2])
	summary_info = open(sys.argv[1])
	output_info = open(sys.argv[2],"a")
	output_name = sys.argv[3]
	if head_need:
		output_info.write("SampleName\tLayout\tTotalPairs\tUnmappedPairs\tUniqMappedPairs\tMultiMappedPairs\tDiscordantlyAlignedPairs\tUnpairedReads\tUnmappedReads\tUniqMappedReads\tMultiMappedReads\tOverallAlignmentRate\tUniqAlignmentRate\n")
	mappingSummary(summary_info,output_info,output_name)

	summary_info.close()
	output_info.close()
main()