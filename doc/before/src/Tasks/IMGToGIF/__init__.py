from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, argStrTemplate, argTimeTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
from copy import deepcopy
from urllib.parse import unquote
import numpy as np
import imageio
import logging
import time
import json
import re
import os

class Task(taskTemplate):
    """ IMGToGIF

    将批量PNG/JPG/JPEG图片转为GIF

    ### 一、输入参数说明
    
    #### 1. 数据：`dir_imgs`
    
    png图片所在文件夹
    
    ### 二、输出文件说明
    
    各个下载的文件

    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "04.Tool"
        
    def run(
            self,
            duration: argIntOrFloatTemplate(
                arg_name="duration",
                arg_description="相邻图片之前的时间间隔，单位毫秒",
                arg_default=50,
                arg_min=1,
                arg_max=100000,
                arg_step=1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="dir_imgs",
                    data_description="需要合并的PNG图片所在的文件夹"
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        dir_imgs = kwargs['return'][0]
        assert os.path.exists(dir_imgs), "dir_imgs文件夹不存在"
        
        dir_name = os.path.basename(dir_imgs)

        file_all = list(filter(
            lambda d: d.lower().endswith((".png", ".jpg", ".jpeg")),
            os.listdir(dir_imgs)
        ))
        try:
            file_all.sort(key=lambda d: int(d[len(dir_name):-4]))
        except:
            file_all.sort()

        im_list = list()
        for file_single in file_all:
            im_temp = imageio.imread(os.path.join(dir_imgs, file_single))
            if file_single.lower().endswith(".png"):
                im_temp = im_temp * np.array([1, 1, 1, 0], dtype=np.uint8) + np.array([0, 0, 0, 255], dtype=np.uint8)
            im_list.append(im_temp)
        imageio.mimsave(
            "{}.gif".format(dir_name),
            im_list,
            "GIF",
            loop=0,
            duration=float(duration),
        )

                


