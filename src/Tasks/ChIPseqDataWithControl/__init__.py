from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
import pandas as pd
import logging
import os

class Task(taskTemplate):
    """ ChIPseqDataWithControl

    ChIPseqDataWithControl

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    #### 2. 数据：`meta`
    
    样本的配对关系，excel文件，两列：treat、ctrl，每一行为一对配对的样本

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
    
    def parse_meta(self, file_meta: str):
        assert os.path.exists(file_meta), "file meta not exists: {}".format(file_meta)
        df = pd.read_excel(file_meta)
        cols_need = ["treat", "ctrl"]
        assert len(set(cols_need) - set(df.columns)) == 0, "file meta columns missing: need {}".format(cols_need)
        df = df[cols_need]
        sample_paired = list()
        sample_all = list()
        for _, se in df.iterrows():
            sample_paired.append((se["treat"], se["ctrl"]))
            sample_all.append(se["treat"])
            sample_all.append(se["ctrl"])
        return sample_paired, set(sample_all)
        

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
                dataTemplate(
                    data_name="meta",
                    data_description="样本的配对关系，excel文件，两列：treat和ctrl，每一行为一对配对的样本"
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        data_species = get_pipline_data(species)
        samples_dict = get_pipline_samples(kwargs['return'][0])

        sample_paired, sample_all = self.parse_meta(kwargs['return'][1].strip())
        sample_gap = sample_all - set(samples_dict.keys())
        assert len(sample_gap) == 0, "sample missing: {}".format(sample_gap)

        # 分开算每个样本
        for k in sample_all:
            cmd = [
                os.path.join(os.path.dirname(__file__), "RunChIPseqWithCtrl.split.sh"),
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
                data_species["chromInfo"]
            ]
            self.exec_cmd(" ".join(cmd), logger)

        # 样本合并
        for k in sample_paired:
            cmd = [
                os.path.join(os.path.dirname(__file__), "RunChIPseqWithCtrl.Combine.sh"),
                os.getcwd(),
                k[0],
                k[1],
                data_species["effective_genome_size"],
                data_species["promoter_file"],
                data_species["chromInfo"],
                str(shiftsize)
            ]
            self.exec_cmd(" ".join(cmd), logger)

        logger.info("done")