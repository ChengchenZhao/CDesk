from Tasks import returnTemplate, argStrTemplate

from .Task_GeneExpressionInJiazi import TaskSrc # type: ignore
import os

class Task(TaskSrc):
    def run(
            self,
            gene_symbols: argStrTemplate(
                arg_type="TextArea",
                arg_name="gene_symbols",
                arg_description="需要查看的基因，每一行为一个基因，要求为基因symbol",
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=tuple(), # type: ignore
            data_out=tuple()
        ):
            self.main(
                gene_list=list(filter(
                    lambda dd: len(dd) > 0,
                    map(
                        lambda d: d.strip(),
                        gene_symbols.split("\n")
                    )
                )),
                dir_base=os.getcwd(),
                file_jiazi="/mnt/liudong/data/Genome/jiazi.h5ad",
                dir_2d="/mnt/liudong/data/Genome/jiazi_2d",
                **kwargs
            )

Task.__doc__ = TaskSrc.__doc__
