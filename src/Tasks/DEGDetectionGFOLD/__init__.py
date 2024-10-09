from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, get_pipline_species

from .Task_DEGDetectionGFOLD import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            is_enrich: argSelectTemplate(
                arg_name="富集分析",
                arg_description="是否进行富集分析",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            cutoff: argIntOrFloatTemplate(
                arg_name="GFOLD阈值",
                arg_description="gfold阈值",
                arg_default=1.0,
                arg_min=0.1,
                arg_max=10.0,
                arg_step=0.1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="bam_or_sam_dir",
                    data_description="Bam或Sam文件夹，样本的bam或sam文件所在的文件夹，多个文件夹可以用英文逗号间隔开",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_meta",
                    data_description="根据meta文件将样本分组计算差异表达基因"
                ), # type: ignore
            ),
            data_out=tuple()
        ):
            self.main(
                species=species,
                is_enrich=is_enrich,
                cutoff=float(cutoff),
                **kwargs
            )
Task.__doc__ = TaskSrc.__doc__
