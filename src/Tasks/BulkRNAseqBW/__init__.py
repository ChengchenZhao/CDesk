from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
import logging
import time
import json
import re
import os

class Task(taskTemplate):
    """ BulkRNAseqBW

    BulkRNAseqData生成bigwig

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明
    
    ### 三、任务作者

    Chengchen Zhao、Dong LIU
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
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        data_species = get_pipline_data(species)
        samples_dict = get_pipline_samples(kwargs['return'][0])

        for k in samples_dict:
            cmd = [
                os.path.join(os.path.dirname(__file__), "RunRNAseqBW.sh"),
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
        logger.info("done")

