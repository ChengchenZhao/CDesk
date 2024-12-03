from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
import pandas as pd
import logging
import time
import json
import re
import os

class Task(taskTemplate):
    """ Hi-CData

    Hi-CData的pipeline

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`数据文件夹
    
    ### 二、输出文件说明

    ### 三、任务作者

    Dong LIU、Chengchen ZHAO
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
            Fragment: argSelectTemplate(
                arg_name="限制性内切酶",
                arg_description="限制性内切酶种类",
                arg_default="mboi",
                options={
                    "/mnt/liudong/data/20230414_hic/genome/mboi_frag_mm10.bed": "mboi"
                }
            ), # type: ignore
            FragmentSeq: argSelectTemplate(
                arg_name="限制性内切酶序列",
                arg_description="限制性内切酶序列",
                arg_default="true",
                options={
                    "GATCGATC": "mboi(GATCGATC)"
                }
            ), # type: ignore
            ConfOther: argStrTemplate(
                arg_type="TextArea",
                arg_name="其余自定义配置",
                arg_description="每行一个配置，比如MAX_ITER=100，#开头的为注释，不算配置",
                arg_default="#"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="包括样本的数据文件夹，多个文件夹使用;间隔"
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        data_species = get_pipline_data(species)
        samples_dict = get_pipline_samples(kwargs['return'][0])

        # sort data
        logger.info("sort data ...")
        data_path = os.path.join(os.getcwd(), "sort_data")
        def ln_file(file_temp: str, data_path_now: str, replace_str: tuple):
            nonlocal self, logger
            cmd = [
                "ln",
                "-s",
                file_temp,
                os.path.join(data_path_now, os.path.basename(file_temp).replace(*replace_str)),
            ]
            self.exec_cmd(" ".join(cmd), logger)
        for k in samples_dict:
            data_path_now = os.path.join(data_path, k)
            if not os.path.exists(data_path_now):
                os.makedirs(data_path_now)
            if samples_dict[k]['type'] == "pair":
                ln_file(
                    file_temp=samples_dict[k]['data'][0],
                    data_path_now=data_path_now,
                    replace_str=("_1.f", "_R1.f")
                )
                ln_file(
                    file_temp=samples_dict[k]['data'][1],
                    data_path_now=data_path_now,
                    replace_str=("_2.f", "_R2.f")
                )
            else:
                ln_file(
                    file_temp=samples_dict[k]['data'][0],
                    data_path_now=data_path_now,
                    replace_str=("_", "_")
                )

        # config Hi-C pro
        logger.info("config Hi-C pro ...")
        file_config = os.path.join(os.getcwd(), "config-hicpro.txt")
        with open(os.path.join(os.path.dirname(__file__), "config-hicpro.txt")) as f:
            config_template = f.read()
        data_config = {
            "N_CPU": int(os.environ.get("cdesk_task_num", 8)),
            "JOB_NAME": "CDesk_HiC_{}".format(time.strftime("%Y%m%d%H%M%S", time.localtime())),
            "BOWTIE2_IDX_PATH": os.path.dirname(data_species['mapping_index']),
            "REFERENCE_GENOME": os.path.basename(data_species['mapping_index']),
            "GENOME_SIZE": data_species['chromInfo'],
            "GENOME_FRAGMENT": Fragment,
            "LIGATION_SITE": FragmentSeq
        }
        data_config.update({k.strip():v.strip() for k, v in map(
            lambda d: d.split("=", 1),
            filter(
                lambda dd: (not dd.strip().startswith("#")) or ("=" in dd),
                ConfOther.split("\n")
            )
        )})
        for k in data_config:
            config_template = re.sub("{}.*".format(k), "{} = {}".format(k, data_config[k]), config_template)
        with open(file_config, "w", encoding="UTF-8") as f:
            f.write(config_template)

        # run Hi-C pro
        logger.info("run Hi-C pro ...")
        dir_results = os.path.join(os.getcwd(), "result")
        cmd = [
            "HiC-Pro",
            "-c",
            file_config,
            "-i",
            data_path,
            "-o",
            dir_results
        ]
        self.exec_cmd(" ".join(cmd), logger)

        # run juicer
        dir_hic_results = os.path.join(dir_results, "hic_results", "hic")
        dir_hic_data = os.path.join(dir_results, "hic_results", "data")
        if not os.path.exists(dir_hic_results):
            os.makedirs(dir_hic_results)
        for sample_id in os.listdir(dir_hic_data):
            logger.info("run juicer {} ...".format(sample_id))
            cmd = [
                "bash",
                "hicpro2juicebox.sh",
                "-i",
                os.path.join(dir_hic_data, sample_id, "{}.allValidPairs".format(sample_id)),
                "-g",
                data_config["REFERENCE_GENOME"],
                "-j",
                "/mnt/liudong/software/juicer_tools.2.20.00.jar",
                "-o",
                dir_hic_results
            ]
            self.exec_cmd(" ".join(cmd), logger)
        
        logger.info("done")
        