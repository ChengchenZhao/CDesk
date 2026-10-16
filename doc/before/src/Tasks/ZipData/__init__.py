from Tasks import taskTemplate, dataTemplate, returnTemplate, argSelectTemplate
import logging
import time
import os

class Task(taskTemplate):
    """ ZipData

    将文件夹内的所有数据压缩为zip，方便下载结果数据

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`数据文件夹
    
    需要下载的文件夹数据
    
    #### 2. 参数：`link`快捷方式
    
    是否需要将快捷方式转换为文件
    
    ### 二、输出文件说明

    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "04.Tool"
    

    def get_same_dir(self, dirs: list[str]):
        assert len(dirs) > 0, "No dirs"

        if len(dirs) == 1:
            return os.path.dirname(dirs[0]), [os.path.basename(dirs[0])]

        dirs_split = list(map(
            lambda d: d.split(os.sep),
            dirs
        ))
        dir_base = ""
        for i in zip(*dirs_split):
            is_same = True
            for ii in range(1, len(i)):
                if i[ii] == i[0]:
                    continue
                is_same = False
                break
            if is_same:
                dir_base = "{}{}{}".format(dir_base, i[0], os.sep)
            else:
                break
        
        dirs = list(map(
            lambda d: d[len(dir_base):],
            dirs
        ))
        return dir_base, dirs


    def run(
            self,
            link: argSelectTemplate(
                arg_name="快捷方式",
                arg_description="是否需要将快捷方式转换为文件",
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
                    data_name="DataPath",
                    data_description="需要下载的文件夹数据",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']
        dir_result = os.getcwd()

        logger.info("get dir: {}".format(kwargs['return'][0]))
        
        dir_datas = tuple(map(
            lambda d: d.strip(),
            kwargs['return'][0].split(",")
        ))

        for i in dir_datas:
            assert os.path.exists(i), "dir not exists: {}".format(i)
        
        dir_base, dir_datas = self.get_same_dir(dir_datas)


        logger.info("dir: {}, {}".format(dir_base, dir_datas))

        dir_str = " ".join(dir_datas)
        file_result = os.path.join(dir_result, "data.{}.zip".format(time.strftime("%Y.%m.%d.%H.%M.%S", time.localtime())))

        logger.info("start ...")
        if str(link) == "true":
            cmd = "cd {} && zip -r {} {}".format(dir_base, file_result, dir_str)
        else:
            cmd = "cd {} && zip -r -y {} {}".format(dir_base, file_result, dir_str)
        
        self.exec_cmd(cmd, logger)
        logger.info("result data: {}".format(file_result))
        