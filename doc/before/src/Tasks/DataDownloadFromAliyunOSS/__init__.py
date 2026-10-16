from Tasks import dataTemplate, returnTemplate, argSelectTemplate, argStrTemplate

from .Task_DataDownloadFromAliyunOSS import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            access_mode: argSelectTemplate(
                arg_name="鉴权模式",
                arg_description="登录账号模式",
                arg_default="AK",
                options={
                    "AK": "AK"
                }
            ), # type: ignore
            endpoint: argStrTemplate(
                arg_type="TextArea",
                arg_name="Endpoint",
                arg_description="填写Bucket所在地域的Endpoint，如https://oss-cn-hangzhou.aliyuncs.com",
            ), # type: ignore
            access_keyid: argStrTemplate(
                arg_type="TextArea",
                arg_name="AccessKeyId",
                arg_description="账号的AccessKey",
            ), # type: ignore
            access_keysecret : argStrTemplate(
                arg_type="TextArea",
                arg_name="AccessKeySecret",
                arg_description="账号的AccessKeySecret",
            ), # type: ignore
            path_oss : argStrTemplate(
                arg_type="TextArea",
                arg_name="OSS",
                arg_description="文件的路径，如oss://examplebucket",
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
            access_mode=access_mode,
            endpoint=endpoint,
            access_keyid=access_keyid,
            access_keysecret=access_keysecret,
            path_oss=path_oss,
            **kwargs
        )


Task.__doc__ = TaskSrc.__doc__