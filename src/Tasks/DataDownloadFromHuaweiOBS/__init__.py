from Tasks import dataTemplate, returnTemplate, argStrTemplate

from .Task_DataDownloadFromHuaweiOBS import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            endpoint: argStrTemplate(
                arg_type="TextArea",
                arg_name="Endpoint",
                arg_description="填写Bucket所在地域的Endpoint，如https://obs.cn-north-4.myhuaweicloud.com",
            ), # type: ignore
            access_key_id: argStrTemplate(
                arg_type="TextArea",
                arg_name="AccessKeyId",
                arg_description="账号的AccessKey",
            ), # type: ignore
            secret_access_key : argStrTemplate(
                arg_type="TextArea",
                arg_name="SecretAccessKey",
                arg_description="账号的SecretAccessKey",
            ), # type: ignore
            path_obs : argStrTemplate(
                arg_type="TextArea",
                arg_name="OBS",
                arg_description="文件的路径，如obs://examplebucket",
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=tuple(),
            data_out=(
                dataTemplate(
                    data_name="DownloadDir",
                    data_description="下载的文件夹",
                    data_file="."
                ),
            ) # type: ignore
        ):
        
        self.main(
            endpoint=endpoint,
            access_key_id=access_key_id,
            secret_access_key=secret_access_key,
            path_obs=path_obs,
            **kwargs
        )


Task.__doc__ = TaskSrc.__doc__