from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
import logging
import time
import json
import re
import os

class Task(taskTemplate):
    """ BulkRNAseqData

    BulkRNAseqData的pipeline

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明
    
    #### 1. `mapping_info.txt`
    
    所有样本mapping的结果，mapping率较低则可能需要重新实验
    
    #### 2. `GeneSymbol_counts.txt`和`GeneSymbol_fpkm.txt`
    
    所有样本所有基因的`counts`或`fpkm`，每一列是一个样本，每一行是一个基因

    #### 3. `bam`文件夹等其余文件备用

    ### 三、任务作者

    Chengchen Zhao、Wenchao XIU、Dong LIU
    """
    def init(self) -> None:
        self.label = "02.PreProcess"

    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="数据文件夹"
                ),
            ), # type: ignore
            data_out=(
                dataTemplate(
                    data_name="GeneCounts",
                    data_description="基因counts",
                    data_file="GeneSymbol_counts.txt"
                ),
                dataTemplate(
                    data_name="GeneFpkm",
                    data_description="基因fpkm",
                    data_file="GeneSymbol_fpkm.txt"
                )
            )
        ):
        logger: logging.Logger = kwargs['logger']

        data_species = get_pipline_data(species)
        samples_dict = get_pipline_samples(kwargs['return'][0])

        for k in samples_dict:
            cmd = [
                os.path.join(os.path.dirname(__file__), "RunRNAseq.sh"),
                samples_dict[k]['data_path'],
                os.getcwd(),
                k,
                samples_dict[k]['type'],
                ",".join(map(
                    lambda d: os.path.basename(d),
                    samples_dict[k]["data"]
                )),
                data_species["mapping_index"],
                data_species["refseq_gtf"],
                data_species["refseq_bed"],
                data_species["chromInfo"]
            ]
            self.exec_cmd(" ".join(cmd), logger)
        cmd = [
            os.path.join(os.path.dirname(__file__), "summary.sh"),
            os.getcwd(),
            data_species["refseq_gtf"],
        ]
        self.exec_cmd(" ".join(cmd), logger)

