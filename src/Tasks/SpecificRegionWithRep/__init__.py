from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argStrTemplate

from .Task_SpecificRegionWithRep import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            cutoff: argIntOrFloatTemplate(
                arg_name="合并阈值",
                arg_description="合并重复时，重复region的百分比阈值。即在所有样本内的比例超过一定阈值才算Specific",
                arg_default=20,
                arg_min=0,
                arg_max=100,
                arg_step=1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_meta",
                    data_description="每类样本的bed文件路径。excel文件，存在tag、file两列。tag为样本标签，file为对应标签的bed文件。tag标签只能有两类。",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        self.main(
            cutoff=int(cutoff),
            **kwargs
        )
Task.__doc__ = TaskSrc.__doc__