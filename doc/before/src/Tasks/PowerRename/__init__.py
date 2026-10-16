from Tasks import dataTemplate, returnTemplate, argSelectTemplate, argStrTemplate

from .Task_PowerRename import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            rule_replace: argStrTemplate(
                arg_type="TextArea",
                arg_name="重命名规则",
                arg_description="重命名规则。一共有偶数行，每两行为一组。单数行是需要去除的内容，偶数行是被替换成为的内容",
                arg_default=""
            ), # type: ignore
            use_regex: argSelectTemplate(
                arg_name="使用正则",
                arg_description="使用正则进行替换，否则只把替换的数据当作普通字符",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            include_others: argSelectTemplate(
                arg_name="包含未命名文件",
                arg_description="不参与命名（不符合重命名规则）的文件是否也创建软链接（为了保持文件夹结构，防止部分文件丢失）",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            flatten_dir: argSelectTemplate(
                arg_name="去除各级目录",
                arg_description="去除各级目录，将所有结果文件放在同一个文件夹",
                arg_default="false",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="dir_file",
                    data_description="需要命名的文件所在的文件夹",
                    is_multi=False
                ),
            ), # type: ignore
            data_out=tuple()
        ):
            self.main(
                dir_file=kwargs['return'][0].strip(),
                rule_replace=rule_replace.strip(),
                use_regex=use_regex.lower() == "true",
                include_others=include_others.lower() == "true",
                flatten_dir=flatten_dir.lower() == "true",
                logger=kwargs['logger']
            )
Task.__doc__ = TaskSrc.__doc__

