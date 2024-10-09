from Tasks import taskTemplate, dataTemplate, returnTemplate, argSelectTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
import logging


class Task(taskTemplate):
    """ MNaseData

    MNaseData的pipeline

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明


    ### 三、任务作者

    Chengchen ZHAO、Dong LIU
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
        samples_dict = get_pipline_samples(kwargs['return'][0], need_gzip=False)
        for k in samples_dict:
            file_data = samples_dict[k]['data']
            cmd = [
                "CAM.py",
                "simple",
                "-a",
                file_data[0]
            ]
            if samples_dict[k]["type"] == "pair":
                cmd += [
                    "-b",
                    file_data[1],
                    "-t",
                    "PE"
                ]
            else:
                cmd += [
                    "-t",
                    "SE"
                ]
            cmd += [
                "-n",
                k,
                "-s",
                data_species["version"],
                "--fa",
                data_species["fa"],
                "--mapindex",
                data_species["mapping_index"]
            ]
            self.exec_cmd(" ".join(cmd), logger)

