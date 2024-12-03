from Tasks import taskTemplate, dataTemplate, returnTemplate, argStrTemplate
import pandas as pd
import logging
import os

class Task(taskTemplate):
    """ TableComb

    table文件合并，注意所有文件内不能存在相同的样本

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`数据文件夹
    
    包括所有文件的数据文件夹，多个文件夹使用英文逗号间隔
    
    #### 2. 参数：`data_ext`文件后缀
    
    所有需要的数据文件相同的后缀，根据此后缀识别为需要合并的数据文件，多个后缀之间利用,间隔
    
    #### 3. 参数：`data_no_ext`排除的文件后缀
    
    所有排除的数据文件相同的后缀，根据此后缀识别为不需要合并的数据文件，多个后缀之间利用,间隔
    
    #### 4. 参数：`columns_on`依据列名
    
    依据哪些列合并，每行为一个列名。注意这些列必须在所有数据文件内存在。其余的列将会被合并
    
    ### 二、输出文件说明
    
    #### 1. `Table.combine.txt`
    
    合并之后的数据文件

    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "04.Tool"

    def run(
            self,
            data_ext: argStrTemplate(
                arg_type="Input",
                arg_name="文件后缀",
                arg_description="所有需要的数据文件相同的后缀，根据此后缀识别为需要合并的数据文件，多个后缀之间利用,间隔",
                arg_default="GeneSymbol_fpkm.txt"
            ), # type: ignore
            data_no_ext: argStrTemplate(
                arg_type="Input",
                arg_name="排除的文件后缀",
                arg_description="所有排除的数据文件相同的后缀，根据此后缀识别为不需要合并的数据文件，多个后缀之间利用,间隔",
                arg_default="."
            ), # type: ignore
            columns_on: argStrTemplate(
                arg_type="TextArea",
                arg_name="依据列名",
                arg_description="依据哪些列合并，每行为一个列名。注意这些列必须在所有数据文件内存在。其余的列将会被合并",
                arg_default="gene_symbol"
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="包括所有文件的数据文件夹，多个文件夹使用英文逗号间隔",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        data_ext = tuple(data_ext.split(","))
        data_no_ext = tuple(data_no_ext.split(","))
        file_all = list()
        for data_path_single in kwargs['return'][0].split(","):
            assert os.path.exists(data_path_single), f"Dir not exists: {data_path_single}"
            for root, _, files in os.walk(data_path_single):
                for file_single in files:
                    if not file_single.endswith(data_ext):
                        continue
                    if file_single.endswith(data_no_ext):
                        continue
                    file_all.append(os.path.join(root, file_single))
        
        assert len(file_all) > 0, f"No file endswith: {data_ext}"
        logger.info(f"Got files: {file_all}")

        logger.info(f"Merge file: {file_all[0]}")
        df = pd.read_table(file_all[0])
        cols_on = set(filter(
            lambda d: len(d) > 0,
            map(
                lambda dd: dd.strip(),
                columns_on.split("\n")
            )
        ))
        cols_0 = set(df.columns)
        for i in range(1, len(file_all)):
            logger.info(f"Merge file: {file_all[i]}")
            df_temp = pd.read_table(file_all[i])
            cols_1 = set(df_temp.columns)
            cols_on_temp = cols_1 & cols_on
            cols_samples_gap = (cols_1 & cols_0) - cols_on
            assert len(cols_samples_gap) == 0, f"Duplicate samples: {cols_samples_gap}, {file_all[i]}"
            assert len(cols_on) == len(cols_on_temp), f"file columns error: {file_all[0]}{cols_on} != {file_all[i]}{cols_on_temp}"
            df = pd.merge(
                left=df,
                right=df_temp,
                on=list(cols_on),
                how="outer"
            )
            del df_temp
        df.fillna(0, inplace=True)
        df.to_csv("Table.combine.txt", sep="\t", index=False)
        logger.info("done")
        