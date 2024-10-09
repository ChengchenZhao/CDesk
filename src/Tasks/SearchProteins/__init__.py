from Tasks import returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, get_pipline_species

from .Task_SearchProteins import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            start: argIntOrFloatTemplate(
                arg_name="开始位置",
                arg_description="查找区间的开始位置",
                arg_default=1,
                arg_min=1,
                arg_max=1000000,
                arg_step=1
            ), # type: ignore
            end: argIntOrFloatTemplate(
                arg_name="结束位置",
                arg_description="查找区间的结束位置",
                arg_default=1,
                arg_min=1,
                arg_max=1000000,
                arg_step=1
            ), # type: ignore
            pat: argStrTemplate(
                arg_type="Input",
                arg_name="匹配的模式，多个模式直接使用英文逗号间隔",
                arg_description="使用正则书写匹配的模式",
                arg_default=""
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=tuple(), # type: ignore
            data_out=tuple()
        ):
        self.main(
            species=species,
            start=int(start),
            end=int(end),
            pat=str(pat),
            **kwargs
        )
Task.__doc__ = TaskSrc.__doc__