from Tasks import dataTemplate, returnTemplate, argSelectTemplate, argStrTemplate, get_pipline_species
from .Task_ChangeInRegions import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            bin_background: argStrTemplate(
                arg_type="Input",
                arg_name="BinBackground",
                arg_description="计算背景的bin size，单位kb，不同bin size间使用英文逗号间隔",
                arg_default="1,2,5,10"
            ), # type: ignore
            peak_file_ext: argStrTemplate(
                arg_type="Input",
                arg_name="Peak文件后缀",
                arg_description="Peak文件后缀，比如.e5_peaks",
                arg_default=".e5_peaks"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="metadata",
                    data_description="样本metadata，里面有三列，tag，sample1_ctrl，sample1_treat，每行的tag要求唯一，将两两tag计算change"
                ),
                dataTemplate(
                    data_name="bw_dir",
                    data_description="BigWig文件所在的文件夹，多个文件夹可以用英文分号间隔",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="peak_dir",
                    data_description="peak文件所在的文件夹，多个文件夹可以用英文分号间隔",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="regions",
                    data_description="计算给定regions区域的变化情况",
                    data_default="no",
                    is_multi=True
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        self.main(
            species=species,
            bin_background=bin_background,
            peak_file_ext=peak_file_ext,
            **kwargs
        )


Task.__doc__ = TaskSrc.__doc__