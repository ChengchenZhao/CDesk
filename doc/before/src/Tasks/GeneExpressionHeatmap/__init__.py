from Tasks import dataTemplate, returnTemplate, argStrTemplate

from .Task_GeneExpressionHeatmap import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            palette: argStrTemplate(
                arg_type="TextArea",
                arg_name="颜色列表",
                arg_description="每行一个颜色",
                arg_default="lightblue\n#c3e2a7\nred"
            ), # type: ignore
            sample_order: argStrTemplate(
                arg_type="TextArea",
                arg_name="样本顺序",
                arg_description="需要画热图的样本顺序，每行一个样本，auto为使用所有样本并自动排序",
                arg_default="auto"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_exp",
                    data_description="基因表达矩阵，每一行是一个基因，每一列是一个样本。多个文件自动根据基因合并",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_gene",
                    data_description="需要画热图的基因，文件内每一行为一个基因",
                ),
            ), # type: ignore
            data_out=tuple()
        ):
            self.main(
                palette_ori=palette,
                sample_order_ori=sample_order,
                **kwargs
            )
Task.__doc__ = TaskSrc.__doc__
