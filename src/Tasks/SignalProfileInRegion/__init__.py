from Tasks import dataTemplate, returnTemplate, argIntOrFloatTemplate, argStrTemplate

from .Task_SignalProfileInRegion import TaskSrc # type: ignore

class Task(TaskSrc):
    def run(
            self,
            extend_bp: argIntOrFloatTemplate(
                arg_name="扩大BP(k)",
                arg_description="region两边各扩大多少k bp画图",
                arg_default=5,
                arg_min=1,
                arg_max=100,
                arg_step=1
            ), # type: ignore
            extend_split: argIntOrFloatTemplate(
                arg_name="扩大点数",
                arg_description="region两边各扩大多少k bp画图，这扩大部分取多少个bin计算信号值",
                arg_default=100,
                arg_min=10,
                arg_max=1000,
                arg_step=10
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
                    data_name="file_bed",
                    data_description="需要计算的region，bed文件至少为四列，为chrom/start/end/strand",
                ),
                dataTemplate(
                    data_name="file_bws",
                    data_description="各个bw文件，使用英文逗号间隔",
                    is_multi=True
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        self.main(
            extend_bp=int(extend_bp),
            extend_split=int(extend_split),
            palette_ori=palette,
            **kwargs
        )
Task.__doc__ = TaskSrc.__doc__