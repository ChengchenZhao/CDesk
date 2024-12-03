from Tasks import dataTemplate, returnTemplate

from .Task_DataDownloadFromLinks import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="urls",
                    data_description="需要下载的链接，文件内每一行为一个文件链接"
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
        
        self.main(**kwargs)

Task.__doc__ = TaskSrc.__doc__
