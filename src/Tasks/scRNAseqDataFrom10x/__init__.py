from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
import pandas as pd
import numpy as np
import logging
import gzip
import time
import json
import re
import os

class Task(taskTemplate):
    """ scRNAseqDataFrom10x

    10x scRNAseq 数据预处理

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明
    

    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "02.PreProcess"

    def get_data(self, dir_data: str):
        data = dict()
        
        for root, _, files in os.walk(dir_data):
            for file_single in files:
                if not file_single.endswith(".fastq.gz"):
                    continue
                sample_temp = file_single.rsplit("_S", 1)
                if len(sample_temp) == 1:
                    continue
                sample_name_temp, ext_temp = sample_temp
                if "_L" not in ext_temp:
                    continue
                if "_R" not in ext_temp:
                    continue

                file_single_temp = os.path.join(dir_data, root, file_single)                
                # 文件夹名称为一整个样本数据
                sample_id_temp = os.path.basename(os.path.dirname(file_single_temp)).replace(
                    "-",
                    "_"
                ).replace(
                    ".",
                    "_"
                )
                dir_sample_temp = os.path.join(dir_data, root)
                if sample_id_temp not in data:
                    data[sample_id_temp] = dict()
                    data[sample_id_temp]["dir"] = dir_sample_temp
                    data[sample_id_temp]["sample"] = list()
                data[sample_id_temp]["sample"].append(sample_name_temp)
                assert data[sample_id_temp]["dir"] == dir_sample_temp, "dir error: sample has multi dirs: {}, {}".format(data[sample_id_temp]["dir"], dir_sample_temp)

        for k in data.keys():
            sample_name_temp = list(set(data[k]["sample"]))
            sample_name_temp.sort()
            data[k]["sample"] = sample_name_temp

        return data

    
    def counts_from_10x(self, dir_data: str, file_result: str):
        if os.path.exists(file_result):
            return


        file_barcode = os.path.join(dir_data, "barcodes.tsv.gz")
        file_features = os.path.join(dir_data, "features.tsv.gz")
        file_matrix = os.path.join(dir_data, "matrix.mtx.gz")
        
        for i in (file_barcode, file_features, file_matrix):
            assert os.path.exists(i), "file not exists: {}".format(os.path.basename(i))
        
        barcode = dict()
        with gzip.open(file_barcode) as f:
            n = 1
            for i in f:
                barcode[str(n)] = i.decode().strip()
                n += 1
        
        feature = dict()
        with gzip.open(file_features) as f:
            n = 1
            for i in f:
                feature[str(n)] = i.decode().strip().split("\t", 1)[0]
                n += 1
        
        matrix = dict()
        with gzip.open(file_matrix) as f:
            _ = f.readline()
            _ = f.readline()
            for i in f:
                nline, nrow, exp = i.decode().strip().split(" ", 2)
                feature_temp = feature[nline]
                barcode_temp = barcode[nrow]
                if barcode_temp not in matrix:
                    matrix[barcode_temp] = dict()
                matrix[barcode_temp][feature_temp] = np.uint16(int(exp))
        
        os.makedirs(os.path.dirname(file_result), exist_ok=True)

        df = pd.DataFrame(matrix)
        df.fillna(0, inplace=True)
        for i in df.columns:
            df[i] = df[i].astype(np.uint16)
        df.to_csv(file_result, sep="\t")


    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            include_introns: argSelectTemplate(
                arg_name="include_introns",
                arg_description="include_introns",
                arg_default="false",
                options={
                    "false": "False",
                    "true": "True"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="数据文件夹"
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        data_species = get_pipline_data(species)
        samples_dict = self.get_data(kwargs['return'][0])

        cdesk_task_num = int(os.environ.get("cdesk_task_num", 8))

        for k in samples_dict:
            logger.info("Running: {}, {}".format(k, samples_dict[k]))
            dir_result_temp = os.path.join(os.getcwd(), k)
            cmd = [
                "cd",
                os.getcwd(),
                "&&",
                "cellranger",
                "count",
                "--id={}".format(k),
                "--include-introns={}".format(include_introns),
                "--fastqs={}".format(samples_dict[k]["dir"]),
                "--sample={}".format(",".join(samples_dict[k]["sample"])),
                "--transcriptome={}".format(data_species["scRef10x"]),
                "--localcores={}".format(cdesk_task_num * 2),
                "--localmem={}".format(cdesk_task_num * 8),
                "--localvmem={}".format(cdesk_task_num * 8),
                "--disable-ui > {}.log 2>&1".format(dir_result_temp)
            ]
            self.exec_cmd(" ".join(cmd), logger)
            
            self.exec_cmd("cp -r {0}/outs {0}.temp".format(dir_result_temp), logger)
            self.exec_cmd("rm -rf {}".format(dir_result_temp), logger)
            self.exec_cmd("mv {0}.temp {0}".format(dir_result_temp), logger)

            logging.info("get counts ...")
            self.counts_from_10x(
                dir_data=os.path.join(dir_result_temp, "filtered_feature_bc_matrix"),
                file_result=os.path.join(dir_result_temp, "GeneSymbol_counts.filter.txt")
            )

            logger.info("Running: {} done !".format(k))

