from Tasks import dataTemplate, returnTemplate, get_pipline_samples, argStrTemplate

from .Task_CRISPERScreeningWithMAGeCk import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            sample_tests: argStrTemplate(
                arg_type="TextArea",
                arg_name="SampleTests",
                arg_description="处理的样本",
            ), # type: ignore
            sample_ctrls: argStrTemplate(
                arg_type="TextArea",
                arg_name="SampleCtrls",
                arg_description="对照的样本",
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="数据文件夹"
                ),
                dataTemplate(
                    data_name="file_lib",
                    data_description="sgRNA library file"
                ),
                dataTemplate(
                    data_name="file_none",
                    data_description="negative control sgRNA list"
                ),
            ), # type: ignore
            data_out=(
                dataTemplate(
                    data_name="ResultDir",
                    data_description="结果文件夹",
                    data_file="."
                ),
            )
        ):
    
        samples_dict = get_pipline_samples(kwargs['return'][0])
        name_tests = sample_tests.strip().split("\n")
        name_ctrls = sample_ctrls.strip().split("\n")
        
        sample_gap = set(name_tests + name_ctrls) - set(samples_dict.keys())
        assert len(sample_gap) == 0, "Sample {} not exists, all samples: {}".format(
            sample_gap,
            set(samples_dict.keys())
        )

        file_dict = dict()
        for k in samples_dict:
            file_dict[k] = samples_dict[k]["data"]
            file_dict[k].sort()

        self.main(
            file_lib=kwargs['return'][1].strip(),
            file_none=kwargs['return'][2].strip(),
            file_dict=file_dict,
            name_tests=name_tests,
            name_ctrls=name_ctrls,
            logger=kwargs['logger']
        )


Task.__doc__ = TaskSrc.__doc__