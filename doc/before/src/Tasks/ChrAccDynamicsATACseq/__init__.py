from Tasks import dataTemplate, returnTemplate, argSelectTemplate, argStrTemplate, get_pipline_species

from .Task_ChrAccDynamicsATACseq import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            peak_file_ext: argStrTemplate(
                arg_type="Input",
                arg_name="Peak文件后缀",
                arg_description="Peak文件后缀，比如.e5_peaks",
                arg_default=".e5_peaks"
            ), # type: ignore
            summit_file_ext: argStrTemplate(
                arg_type="Input",
                arg_name="Summit文件后缀",
                arg_description="Summit文件后缀，比如.e5_summits",
                arg_default=".e5_summits"
            ), # type: ignore
            remove_cluster: argStrTemplate(
                arg_type="Input",
                arg_name="隐藏类别",
                arg_description="需要在结果内去除的类别，多个类别用英文分号间隔，比如CO9;OC1",
                arg_default="CO9;OC1"
            ), # type: ignore
            remove_motifs: argStrTemplate(
                arg_type="TextArea",
                arg_name="隐藏Motif",
                arg_description="需要在结果内去除的Motif，每行一个motif，如NF1:FOXA1(CTF,Forkhead)",
                arg_default="无"
            ), # type: ignore
            important_index: argSelectTemplate(
                arg_name="第一个是处理组",
                arg_description="是否metadata内第一个tag_group是处理组，否则的话第一个为对照组",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_meta",
                    data_description="样本metadata，里面有三列，sample/tag_time/tag_group，tag_group限制两类，每一类的tag_time的数量及顺序需要相同"
                ),
                dataTemplate(
                    data_name="dir_bw",
                    data_description="BigWig文件所在的文件夹，多个文件夹可以用英文分号间隔"
                ),
                dataTemplate(
                    data_name="dir_bed",
                    data_description="Bed文件所在的文件夹，多个文件夹可以用英文分号间隔"
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        self.main(
            species=species,
            peak_file_ext=peak_file_ext,
            summit_file_ext=summit_file_ext,
            remove_cluster=remove_cluster,
            remove_motifs=remove_motifs,
            important_index=important_index,
            **kwargs
        )


Task.__doc__ = TaskSrc.__doc__