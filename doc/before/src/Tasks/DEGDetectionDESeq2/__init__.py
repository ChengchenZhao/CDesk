from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, get_pipline_species

from .Task_DEGDetectionDESeq2 import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="species",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            cutoff: argIntOrFloatTemplate(
                arg_name="log2fc",
                arg_description="log2fc阈值",
                arg_default=1.0,
                arg_min=0.1,
                arg_max=10.0,
                arg_step=0.1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_counts",
                    data_description="需要计算样本的counts文件，多个counts文件利用英文逗号间隔，将利用每个文件的第一列合并",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_meta",
                    data_description="根据meta文件将样本分组计算差异表达基因"
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        
        self.main(
            species=species,
            cutoff=float(cutoff),
            **kwargs
        )
    
Task.__doc__ = TaskSrc.__doc__

