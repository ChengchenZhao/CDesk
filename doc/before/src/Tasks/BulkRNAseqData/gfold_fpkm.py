# -*- coding: UTF-8 -*-
# cython: language_level=3

'''
@File    :   gfold_fpkm.py
@Modif   :   2024/01/21
@Author  :   winter
@Version :   0.1
@Contact :   winter_lonely@126.com
@License :   (C)Copyright winter
@Desc    :   
'''

# %%
import concurrent.futures as fu
import pandas as pd
import sys
import os
# %%
dir_bam = sys.argv[1]
dir_result = sys.argv[2]
file_gtf = sys.argv[3]
# %%
assert os.path.exists(file_gtf), "GTF file not exists !"
assert os.path.exists(dir_bam), "Bam dir not exists !"
os.makedirs(dir_result, exist_ok=True)
# %%
# bam转sam
file_sam = list()
file_self = list()
for file_single in os.listdir(dir_bam):
    if file_single.lower().endswith(".sam"):
        file_sam.append(os.path.join(dir_bam, file_single))
        continue
    if not file_single.lower().endswith(".bam"):
        continue
    s_t = os.path.join(dir_bam, file_single)
    s_r = os.path.join(dir_result, "{}.sam".format(file_single[:-4]))
    cmd = [
        "samtools",
        "view",
        "-@",
        # str(os.environ.get("cdesk_task_num", 8)),
        "8",
        "-o",
        s_r,
        s_t
    ]
    os.system(" ".join(cmd))
    file_self.append(s_r)
    file_sam.append(s_r)

        
def task_done(f: fu.Future):
    # 线程异常处理
    e = f.exception()
    if e:
        print(e, flush=True)

# sam转gfold count
task_all = list()
file_counts = list()
with fu.ThreadPoolExecutor(max_workers=8) as executor:
    for file_single in file_sam:
        sample_name = os.path.basename(file_single)[:-4]
        s_t = os.path.join(dir_result, sample_name, "{}.gfold_count".format(sample_name))
        file_counts.append(s_t)
        cmd = [
            "cd",
            os.path.dirname(file_single),
            "&&"
            "gfold",
            "count",
            "-ann",
            file_gtf,
            "-tag",
            os.path.basename(file_single),
            "-o",
            s_t,
            "> {}/{}.gfold_count.log 2>&1".format(os.path.join(dir_result, "log", sample_name), sample_name)
        ]
        # 提交任务
        task_temp = executor.submit(
            os.system,
            " ".join(cmd)
        )
        # 任务回调
        task_temp.add_done_callback(task_done)
        task_all.append(task_temp)
    fu.wait(task_all)
# %%
data_counts = dict()
data_fpkm = dict()
for file_single in file_counts:
    if not file_single.endswith(".gfold_count"):
        continue
    sample_name = os.path.basename(file_single)[:-12]
    data_fpkm[sample_name] = dict()
    data_counts[sample_name] = dict()
    with open(file_single) as f:
        for line in f:
            line = line.strip().split("\t", 4)
            if len(line) != 5:
                continue
            gene_temp, _, count_temp, _, fpkm_temp = line
            data_fpkm[sample_name][gene_temp] = float(fpkm_temp)
            data_counts[sample_name][gene_temp] = int(count_temp)
# %%
def single(data: dict, file_result: str):
    df = pd.DataFrame(data)
    df.index.name = "gene_symbol"
    df.to_csv(
        file_result,
        sep="\t"
    )
# %%
single(data_counts, os.path.join(dir_result, "GeneSymbol_counts.txt"))
single(data_fpkm, os.path.join(dir_result, "GeneSymbol_fpkm.txt"))
# %%
for file_single in file_self:
    os.system("rm -rf \"{}\"".format(file_single))