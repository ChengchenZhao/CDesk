from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples

from .Task_TADChange import TaskSrc # type: ignore


class Task(TaskSrc):
    def run(
            self,
            samples: argStrTemplate(
                arg_type="TextArea",
                arg_name="样本名字",
                arg_description="需要计算的样本的名字，每行一个样本名字，多个样本之间两两比较",
            ), # type: ignore
            resolution: argIntOrFloatTemplate(
                arg_name="分辨率",
                arg_description="domain文件的resolution，单位k",
                arg_default=1000,
                arg_min=1,
                arg_max=10000,
                arg_step=1
            ), # type: ignore
            domains_ext: argStrTemplate(
                arg_type="Input",
                arg_name="domains文件后缀",
                arg_description="domains文件后缀，如_domains.bed。整个文件的命名方式为sample_resolution_ext，如sample1_1m_domains.bed/sample2_10k_domains.bed",
                arg_default="_domains.bed"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="dir_data",
                    data_description="HiCExplorer结果，内包含domain.bed",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        self.main(
            samples=samples,
            resolution=int(resolution),
            domains_ext=domains_ext,
            **kwargs
        )
Task.__doc__ = TaskSrc.__doc__