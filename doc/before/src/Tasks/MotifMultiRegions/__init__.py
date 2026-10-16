from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, get_pipline_species

from .Task_MotifMultiRegions import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
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
                    data_name="bed_files",
                    data_description="查找motif的bed文件，使用英文逗号间隔开",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
            self.main(
                species=species,
                remove_motifs=set(remove_motifs.strip().split("\n")),
                motif_length=motif_length.strip().split(","),
                cutoff_log10_pval=float(cutoff_log10_pval),
                cutoff_log10_pval_star=float(cutoff_log10_pval_star),
                cutoff_ratio=float(cutoff_ratio),
                cutoff_ratio_star=float(cutoff_ratio_star),
                **kwargs
            )
Task.__doc__ = TaskSrc.__doc__

