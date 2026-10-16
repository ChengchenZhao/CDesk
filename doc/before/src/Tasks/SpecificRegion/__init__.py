from Tasks import dataTemplate, returnTemplate, argSelectTemplate, argStrTemplate

from .Task_SpecificRegion import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            region_names: argStrTemplate(
                arg_type="TextArea",
                arg_name="Region名字",
                arg_description="每类regin的名字，每行为一个name，对应每个bed文件，auto为自动根据bed文件命名，利用.-_等分隔，去除前后相同的，即为名字",
                arg_default="auto"
            ), # type: ignore
            name_align: argStrTemplate(
                arg_type="Input",
                arg_name="标准Region类别",
                arg_description="以某一类region为标准，其余的做对比。no标准region则需要只有两类region，相互对比",
                arg_default="no"
            ), # type: ignore
            title: argStrTemplate(
                arg_type="Input",
                arg_name="结果名字",
                arg_description="结果名字，如果为无则自动根据类别名字计算",
                arg_default="no"
            ), # type: ignore
            palette: argStrTemplate(
                arg_type="TextArea",
                arg_name="颜色列表",
                arg_description="每行一个颜色",
                arg_default="#CE0013\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75\n#1E7CAF\n#46489A\n#0F231F\n#1187CD\n#b7704f\n#596032\n#2570a1\n#281f1d\n#de773f\n#525f42\n#2585a6\n#2f271d\n#c99979\n#5f5d46\n#1b315e\n#1d1626\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75"
            ), # type: ignore
            calc_cor: argSelectTemplate(
                arg_name="计算相关性",
                arg_description="是否计算相关性，根据peak计算bw的相关性",
                arg_default="true",
                options={
                    "true": "True",
                    "false": "False"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_beds",
                    data_description="各个bed文件，使用英文逗号间隔",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_summits",
                    data_description="各个bed文件对应的summits，两个region直接对比时，是否要用summit对齐，使用英文逗号间隔",
                    data_default="no",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_bws",
                    data_description="各个bed文件对应的bigWig文件，使用英文逗号间隔。如果多个bw，则画热图时是否需要其余的bigwig一起画",
                    data_default="no",
                    is_multi=True
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        self.main(
            region_names=region_names,
            name_align=name_align,
            title=title,
            palette_ori=palette,
            calc_cor=calc_cor,
            **kwargs
        )
Task.__doc__ = TaskSrc.__doc__