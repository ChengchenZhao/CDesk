from Tasks import dataTemplate, returnTemplate, argSelectTemplate, get_pipline_species

from .Task_FunctionalEnrichment import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="GeneSymbols",
                    data_description="需要富集的基因，文件内每一行为一个基因，要求为基因symbol。多类基因放在不同文件内，使用英文逗号间隔",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
            self.main(
                species=species,
                **kwargs
            )
Task.__doc__ = TaskSrc.__doc__
