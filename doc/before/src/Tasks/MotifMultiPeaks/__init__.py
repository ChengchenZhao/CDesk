from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, get_pipline_species

from .Task_MotifMultiPeaks import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            samples: argStrTemplate(
                arg_type="TextArea",
                arg_name="样本列表",
                arg_description="需要计算的样本列表，每行一个样本名字，样本名字为文件名字内“Peak文件后缀”前的部分",
                arg_default=""
            ), # type: ignore
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            summit_extend_bp: argIntOrFloatTemplate(
                arg_name="Summit扩大BP",
                arg_description="summit左右分别扩大多少BP计算Motif",
                arg_default=25,
                arg_min=0,
                arg_max=100000,
                arg_step=1
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
                arg_description="Summit文件后缀.e5_peaks",
                arg_default=".e5_summits"
            ), # type: ignore
            remove_motifs: argStrTemplate(
                arg_type="TextArea",
                arg_name="隐藏Motif",
                arg_description="需要在结果内去除的Motif，每行一个motif，如NF1:FOXA1(CTF,Forkhead)",
                arg_default="无"
            ), # type: ignore
            motif_length: argStrTemplate(
                arg_type="Input",
                arg_name="Motif长度",
                arg_description="需要查找的Motif长度，不同长度间使用英文逗号间隔",
                arg_default="8,10,12"
            ), # type: ignore
            cutoff_log10_pval: argIntOrFloatTemplate(
                arg_name="数据阈值-log10(pval)",
                arg_description="数据阈值-log10(pval)，决定数据是否画图展示，pval表示找到的该motif是否可信",
                arg_default=2.0,
                arg_min=0.1,
                arg_max=20.0,
                arg_step=0.1
            ), # type: ignore
            cutoff_log10_pval_star: argIntOrFloatTemplate(
                arg_name="星号阈值-log10(pval)",
                arg_description="星号阈值-log10(pval)，决定数据是否画图时标注星号，pval表示找到的该motif是否可信",
                arg_default=5.0,
                arg_min=0.1,
                arg_max=20.0,
                arg_step=0.1
            ), # type: ignore
            cutoff_ratio: argIntOrFloatTemplate(
                arg_name="数据阈值ratio",
                arg_description="数据阈值ratio，决定数据是否画图展示，ratio表示motif是否特异，即 % of Targets / % of Background",
                arg_default=2.0,
                arg_min=0.1,
                arg_max=20.0,
                arg_step=0.1
            ), # type: ignore
            cutoff_ratio_star: argIntOrFloatTemplate(
                arg_name="星号阈值ratio",
                arg_description="星号阈值ratio，决定数据是否画图时标注星号，ratio表示motif是否特异，即 % of Targets / % of Background",
                arg_default=3.0,
                arg_min=0.1,
                arg_max=20.0,
                arg_step=0.1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="peak_dir",
                    data_description="peak文件所在的文件夹，多个文件夹可以用英文逗号间隔。里面包含成对的peak和summit文件",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        self.main(
            samples=samples.strip().split("\n"),
            species=species,
            extend_bp=int(summit_extend_bp),
            peak_file_ext=peak_file_ext,
            summit_file_ext=summit_file_ext,
            remove_motifs=set(remove_motifs.strip().split("\n")),
            motif_length=motif_length,
            cutoff_log10_pval=float(cutoff_log10_pval),
            cutoff_log10_pval_star=float(cutoff_log10_pval_star),
            cutoff_ratio=float(cutoff_ratio),
            cutoff_ratio_star=float(cutoff_ratio_star),
            **kwargs
        )
    
Task.__doc__ = TaskSrc.__doc__


