from Tasks import dataTemplate, returnTemplate, argStrTemplate

from .Task_CorrelationForSamples import TaskSrc # type: ignore


class Task(TaskSrc):
    def run(
            self,
            palette: argStrTemplate(
                arg_type="TextArea",
                arg_name="颜色列表",
                arg_description="每行一个颜色",
                arg_default="#CE0013\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75\n#1E7CAF\n#46489A\n#0F231F\n#1187CD\n#b7704f\n#596032\n#2570a1\n#281f1d\n#de773f\n#525f42\n#2585a6\n#2f271d\n#c99979\n#5f5d46\n#1b315e\n#1d1626\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75"
            ), # type: ignore
            sample_order: argStrTemplate(
                arg_type="TextArea",
                arg_name="样本顺序",
                arg_description="每行一个样本id，相关性热图内的样本顺序，auto为自动排序",
                arg_default="auto"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_fpkm",
                    data_description="fpkm文件。transcript fpkm时，前几列为chr、strand、start、end、refseq、num_exons、length、gene_symbol。genesymbol_fpkm时第一列为gene_symbol",
                    is_multi=True
                ), # type: ignore
                dataTemplate(
                    data_name="file_meta",
                    data_description="metadata的excel文件，存在experiment、sample和tag三列"
                ),
                dataTemplate(
                    data_name="file_genes",
                    data_description="利用给定的基因计算相关性，每个txt文件为一个基因列表，文件每行一个基因。默认为all，使用所有基因",
                    data_default="all",
                    is_multi=True
                )
            ),
            data_out=(
                dataTemplate(
                    data_name="ResultDir",
                    data_description="结果文件夹",
                    data_file="."
                ),
            )
        ):

        self.main(
            palette=self.handle_str_list(palette),
            sample_order=sample_order,
            **kwargs
        )

Task.__doc__ = TaskSrc.__doc__