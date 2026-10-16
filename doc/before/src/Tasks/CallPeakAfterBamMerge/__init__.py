# -*- coding: UTF-8 -*-
# cython: language_level=3

'''
@File    :   setup.py
@Modif   :   2024/01/21
@Author  :   winter
@Version :   0.1
@Contact :   winter_lonely@126.com
@License :   (C)Copyright winter
@Desc    :   
'''
from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, get_pipline_species
from .Task_CallPeakAfterBamMerge import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            name: argStrTemplate(
                arg_type="Input",
                arg_name="结果名字",
                arg_description="合并后的样本名字",
                arg_default="sample"
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
                    data_name="Bams",
                    data_description="需要合并的bam文件",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        self.main(
            species=species,
            name=name,
            shiftsize=shiftsize,
            **kwargs
        )

Task.__doc__ = TaskSrc.__doc__