import os,sys
import math
import subprocess

def main():
	bed_file = sys.argv[1]
	bw_file = sys.argv[2]
	out_file = sys.argv[3]
	steps = int(sys.argv[4])

	bed_info = open(bed_file,"r")
	out_info = open(out_file,"w")

	for each in bed_info:
		each = each.strip().split()
		chr = each[0]
		start = int(each[1])
		end = int(each[2])
		cmd = "bigWigSummary {0} {1} {2} {3} {4}".format(bw_file,chr,start,end,steps)
		pipe = subprocess.Popen(cmd,stdout=subprocess.PIPE,shell='TRUE',universal_newlines='TRUE')
		(result,err) = pipe.communicate()

		if not result:
			out_info.write("\t".join(["NA"]*steps)+"\n")
		else:
			result = result.strip().split("\t")
			out_info.write("\t".join(result)+"\n")
	bed_info.close()
	out_info.close()

if __name__ == "__main__":
	main()