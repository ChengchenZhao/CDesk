from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, get_pipline_species

from .Task_DEGDetectionFC import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            is_enrich: argSelectTemplate(
                arg_name="富集分析",
                arg_description="是否进行富集分析",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            log2fc: argIntOrFloatTemplate(
                arg_name="log2fc",
                arg_description="log2FoldChange阈值",
                arg_default=1.0,
                arg_min=0.1,
                arg_max=20.0,
                arg_step=0.1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_exp",
                    data_description="基因表达矩阵txt，每一行是一个基因，每一列是一个样本。多个表达矩阵文件可以用英文逗号间隔开，不同表达矩阵文件利用基因合并",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_meta",
                    data_description="根据meta文件将样本分组计算差异表达基因"
                ),
            ), # type: ignore
            data_out=tuple()
        ):
            self.main(
                  species=species,
                  is_enrich=is_enrich,
                  log2fc=float(log2fc),
                  **kwargs
            )
Task.__doc__ = TaskSrc.__doc__
