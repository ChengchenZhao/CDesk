#!usr/bin/python
#coding:utf-8
'''
stringtie output ctab data process
by xiuwenchao
modified by Chengchen Zhao on 20210712
'''

import warnings
warnings.filterwarnings("ignore")

import subprocess
import pandas as pd
import os,sys

def get_mapping_file_path(path):
	'''
	input:log path,
	output:log path list
	'''
	sample_names = subprocess.getoutput('ls %s'%path).split('\n')
	log_file_global_path_list = []
	for name in sample_names:
		log_file_global_path = '/'.join([path,name,'%s.Mapping.log'%name])
		if os.path.exists(log_file_global_path):
			log_file_global_path_list.append(log_file_global_path)
	return log_file_global_path_list

def get_mapping_info(mapping_log_path_list,filename):
	'''
	input: mapping log path
	output: mapping information file
	'''
	fout = open(filename,'w')
	print('\t'.join(['sample_name','total_read_pairs','mapped_read_pairs','unique_read_pairs','unique_read_pairs_percentage']),file=fout)
	for mapping_log_path in mapping_log_path_list:
		line_num = 1
		sample_name = '_'.join(mapping_log_path.strip().split('/')[-2:-1])
		for line in open(mapping_log_path,'r').readlines():
			if line_num==1:
				total_read_pairs = line.strip().split(' ')[0]
			if line_num==4:
				unique_read_pairs = line.strip().split(' ')[0]
				unique_read_pairs_percentage = line.strip().split('%')[0].split('(')[1]+'%'
			# if layout == "pair" and line_num==15: # pair end
			# 	mapped_read_pairs = int(round(float(line.strip().split('%')[0])*0.01*float(total_read_pairs)))
			# if layout == "single" and line_num==6: # single end
			# 	mapped_read_pairs = int(round(float(line.strip().split('%')[0])*0.01*float(total_read_pairs)))
			if line_num==6:
				if line.strip() == "----": # pair end
					pass
				else: # single end
					mapped_read_pairs = int(round(float(line.strip().split('%')[0])*0.01*float(total_read_pairs)))
			if  line_num==15:
				mapped_read_pairs = int(round(float(line.strip().split('%')[0])*0.01*float(total_read_pairs)))
			line_num += 1
		print('\t'.join([sample_name,total_read_pairs,str(mapped_read_pairs),unique_read_pairs,unique_read_pairs_percentage]),file=fout)
	fout.close()

def get_readdis_file_path(path):
	'''
	input:read distribution path,
	output:log path list
	'''
	sample_names = subprocess.getoutput('ls %s'%path).split('\n')
	readdis_file_global_path_list = []
	for name in sample_names:
		readdis_file_global_path = '/'.join([path,name,'%s.readsDistribution'%name])
		if os.path.exists(readdis_file_global_path):
			readdis_file_global_path_list.append(readdis_file_global_path)
	return readdis_file_global_path_list

def get_readdis_info(readdis_file_path_list,filename):
	'''
	input:read distribution path list,
	output:merge read distribution table
	'''
	fout = open(filename,'w')
	print('\t'.join(['sample_name','total_tags','CDS_Exons','fiveUTR_Exons','threeUTR_Exons','total_exons','introns','exon_percentage','intron_percentage']),file=fout)
	for readdis_file_path in readdis_file_path_list:
		sample_name = '_'.join(readdis_file_path.strip().split('/')[-3:-1])
		for line in open(readdis_file_path,'r').readlines():
			if line.startswith('Total Tags'):
				total_tags = line.strip().split(' ')[-1]
			if line.startswith('CDS_Exons'):
				CDS_Exons = [i for i in line.strip().split(' ') if i!=''][2]
			if line.startswith("5'UTR_Exons"):
				fiveUTR_Exons = [i for i in line.strip().split(' ') if i!=''][2]
			if line.startswith("3'UTR_Exons"):
				threeUTR_Exons = [i for i in line.strip().split(' ') if i!=''][2]
			if line.startswith('Introns'):
				introns = [i for i in line.strip().split(' ') if i!=''][2]
				total_exons = str(int(CDS_Exons)+int(fiveUTR_Exons)+int(threeUTR_Exons))
				exon_percentage = str(round((float(CDS_Exons)+float(fiveUTR_Exons)+float(threeUTR_Exons))/float(total_tags)*100,2))+'%'
				intron_percentage = str(round(float(introns)/float(total_tags)*100,2))+'%'
		print('\t'.join([sample_name,total_tags,CDS_Exons,fiveUTR_Exons,threeUTR_Exons,total_exons,introns,exon_percentage,intron_percentage]),file=fout)
	fout.close()

	
def get_ctab_file_path(path):
	'''
	input:exp path,
	output:ctab path list
	'''
	sample_names = subprocess.getoutput('ls %s'%path).split('\n')
	exp_file_global_path_list = []
	for name in sample_names:
		exp_file_global_path = '/'.join([path,name,'t_data.ctab'])
		if os.path.exists(exp_file_global_path):
			exp_file_global_path_list.append(exp_file_global_path)
	return exp_file_global_path_list

def merge_ctab_files(ctab_path_list,merge_file_name):
	'''
	input:ctab file absolute path list, output merge file names
	output: merge FPKM file
	'''
	first_file = True
	for ctab_path in ctab_path_list:
		sample_name = '_'.join(ctab_path.split('/')[-2:-1])
		if first_file:
			ctab = pd.read_csv(ctab_path,header=0,sep='\t')
			first_file = False
			ctab_merge = ctab[['chr','strand','start','end','t_name','num_exons','length','gene_id','FPKM']].rename(
				columns={'t_name':'refseq','gene_id':'gene_symbol','FPKM':sample_name})
		else:
			ctab = pd.read_csv(ctab_path,header=0,sep='\t')
			ctab_merge[sample_name] = ctab['FPKM']
	ctab_merge.to_csv(merge_file_name,index=False,sep='\t')


def create_gtf_path_file(paths_list,file_name):
	'''
	input:exp path,
	output:gtf path txt file
	'''
	names_list = []
	gtf_path_list = []
	for path in paths_list:
		sample_names = subprocess.getoutput('ls %s'%path).split('\n')
		exp_file_global_path_list = []
		for name in sample_names:
			exp_file_global_path = '/'.join([path,name,'%s.gtf'%name])
			if os.path.exists(exp_file_global_path):
				exp_file_global_path_list.append(exp_file_global_path)
				# names_list += [name+'_count']
				names_list.append(name)
		gtf_path_list += exp_file_global_path_list
	gtf_path_file = pd.DataFrame({'names_list':names_list,'gtf_path_list':gtf_path_list})
	gtf_path_file[['names_list','gtf_path_list']].to_csv(file_name,index=False,header=False,sep='\t')


path = sys.argv[1]

# 1.sample mapping information
mapping_file_path = '{}/log'.format(path)
mapping_file_paths = [mapping_file_path]

mapping_path_list = []
for mapping_file_path in mapping_file_paths:
	mapping_path_list += get_mapping_file_path(mapping_file_path)

get_mapping_info(mapping_path_list,'{}/mapping_info.txt'.format(path))

# 2.read distribution extract
readdis_file_path = '{}/log'.format(path)
readdis_file_paths = [readdis_file_path]

readdis_path_list = []
for readdis_file_path in readdis_file_paths:
	readdis_path_list += get_readdis_file_path(readdis_file_path)

get_readdis_info(readdis_path_list,'{}/read_distribution_info.txt'.format(path))

# 3.merge ctab files
ctab_file_path = '{}/'.format(path)
ctab_file_paths = [ctab_file_path]

ctab_path_list = []
for ctab_file_path in ctab_file_paths:
	ctab_path_list += get_ctab_file_path(ctab_file_path)

merge_ctab_files(ctab_path_list,'{}/merge_fpkm.txt'.format(path))


# 4.create gtf gtf_path file
gtf_file_path1 = '{}/'.format(path)
# gtf_file_paths = [gtf_file_path1,gtf_file_path2,gtf_file_path3]
gtf_file_paths = [gtf_file_path1]

create_gtf_path_file(gtf_file_paths,'{}/gtf_file_path.txt'.format(path))




