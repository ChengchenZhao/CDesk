from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
import logging
import time
import json
import re
import os

class Task(taskTemplate):
    """ ChIPseqDataWOControl

    ChIPseqDataWOControl的pipeline

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明
    
    #### 1. `MappingSummary_ChIPseq.txt`
    
    所有样本mapping的结果，mapping率较低则可能需要重新实验

    #### 2. `bw`和`bam`文件夹

    每个样本的`bigwig`或`bam`文件
    
    #### 3. `peak`文件夹

    包含所有样本的`peak`及`summits`

    其中每个样本根据`Fragments.bed`使用不同`p-value`(`1e-3`、`1e-5`、`1e-7`)跑一份分别为`e3`、`e5`、`e7`，`p-value`越小卡的越严，出来的peak越少
    
    并用`UniqSEReadsSorted.bed`利用`p-value=1e-7`跑一份为`HC`

    #### 4. 各个样本文件夹等其余文件备用

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
            shiftsize: argIntOrFloatTemplate(
                arg_name="PeakShiftSize",
                arg_description="小于多少bp的peak将被合并",
                arg_default=73,
                arg_min=1,
                arg_max=1000000,
                arg_step=1
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
                os.path.join(os.path.dirname(__file__), "RunChIPseqWOCtrl.sh"),
                samples_dict[k]['data_path'],
                os.getcwd(),
                k,
                samples_dict[k]['type'],
                ",".join(map(
                    lambda d: os.path.basename(d),
                    samples_dict[k]["data"]
                )),
                data_species["mapping_index"],
                data_species["effective_genome_size"],
                data_species["promoter_file"],
                data_species["chromInfo"],
                str(shiftsize)
            ]
            self.exec_cmd(" ".join(cmd), logger)

        logger.info("done")