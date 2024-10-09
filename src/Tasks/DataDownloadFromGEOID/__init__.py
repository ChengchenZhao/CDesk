from Tasks import dataTemplate, returnTemplate, argSelectTemplate

from .Task_DataDownloadFromGEOID import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            download_type: argSelectTemplate(
                arg_name="下载类型",
                arg_description="需要下载的数据类型",
                arg_default="fastq",
                options={
                    "fastq": "Fastq",
                    "file": "SupplementaryFiles"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="IDs",
                    data_description="需要下载的GSE_ID/GSM_ID，文件每行一个ID"
                ),
            ), # type: ignore
            data_out=(
                dataTemplate(
                    data_name="DownloadDir",
                    data_description="下载的文件夹",
                    data_file="."
                ),
            )
        ):
        self.main(
            download_type=download_type,
            **kwargs
        )

Task.__doc__ = TaskSrc.__doc__

                


