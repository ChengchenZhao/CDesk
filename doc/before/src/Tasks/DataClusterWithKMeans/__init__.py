from Tasks import dataTemplate, returnTemplate, argStrTemplate, argIntOrFloatTemplate, argSelectTemplate

from .Task_DataClusterWithKMeans import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            n_clusters: argIntOrFloatTemplate(
                arg_name="n_clusters",
                arg_description="需要聚为多少个类别",
                arg_default=10,
                arg_min=1,
                arg_max=100,
                arg_step=1
            ), # type: ignore
            zscore: argSelectTemplate(
                arg_name="zscore",
                arg_description="数据是否需要进行Z-Score",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            log2: argSelectTemplate(
                arg_name="log2",
                arg_description="数据是否需要进行log2(data+1)",
                arg_default="true",
                options={
                    "true": "是",
                    "false": "否"
                }
            ), # type: ignore
            method: argSelectTemplate(
                arg_name="method",
                arg_description="聚类阈值的三种方法: cos(cos) eu(euclidean metriccor) cor(pearson correlation)",
                arg_default="eu",
                options={
                    "eu": "欧式距离",
                    "cor": "皮尔逊相关系数",
                    "cos": "余弦相关性",
                }
            ), # type: ignore
            palette: argStrTemplate(
                arg_type="TextArea",
                arg_name="颜色列表",
                arg_description="每行一个颜色",
                arg_default="#CE0013\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75\n#1E7CAF\n#46489A\n#0F231F\n#1187CD\n#b7704f\n#596032\n#2570a1\n#281f1d\n#de773f\n#525f42\n#2585a6\n#2f271d\n#c99979\n#5f5d46\n#1b315e\n#1d1626\n#16557A\n#C7A609\n#87C232\n#008792\n#A14C94\n#15A08C\n#8B7E75"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_data",
                    data_description="需要聚类的数据文件，除了第一列，其余列为样本。如果输入多个则根据第一列自动进行合并",
                    is_multi=True
                ),
                dataTemplate(
                    data_name="file_meta",
                    data_description="样本metadata，里面有三列，sample/tag_time/tag_group，tag_group为每一类处理，每一类的tag_time的数量及顺序需要相同"
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

        self.main(
            file_data=kwargs["return"][0].strip().split(","),
            file_meta=kwargs["return"][1].strip(),
            n_clusters=n_clusters,
            zscore=zscore,
            log2=log2,
            method=method,
            palette=self.handle_str_list(palette),
            logger=kwargs["logger"]
        )


Task.__doc__ = TaskSrc.__doc__