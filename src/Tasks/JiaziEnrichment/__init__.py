from Tasks import taskTemplate, dataTemplate, returnTemplate
from matplotlib.font_manager import fontManager
import matplotlib as mpl
import logging
import os

mpl.use("Agg")
fontManager.addfont(os.path.join(os.path.dirname(__file__), "..", "arial.ttf"))
mpl.rc('pdf', fonttype=42)
mpl.rcParams['font.sans-serif'] = ["Arial"]


class Task(taskTemplate):
    """ JiaziEnrichment

    根据基因列表文件和甲子数据进行富集分析，与GO富集分析类似，只是结果term换成了甲子数据标注的cell type

    ### 一、输入参数说明
    
    #### 1. 数据：`GeneSymbols`
    
    需要富集的基因，文件内每一行为一个基因，要求为基因symbol，可以根据输入框右边的辅助按钮选择文件

    多类基因放在不同文件内，使用英文逗号间隔
    
    ### 二、输出文件说明
    
    1、`txt`为富集结果，`pdf`为对应的画图

    2、画图结果为前12个符合的term，p-value闽值为0.05，qvalue阔值为0.2，最少基因数量为10(Count = 10)，因此可能存在空的画图结果
    
    3、富集结果字段说明
    
      - pvalue：富集的p值
      
      - p.adjust：校正之后的p值
      
      - qvalue：q值
      
      - genes：输入的做富集分析的基因中富集到这个term上面的基因
      
      - Count：输入的做富集分析的基因中富集到这个term上面的基因的数目（genes的数量，dotplot内点的大小）
      
      - GeneRatio：这里是一个分数，分子是富集到这个term上的基因的数目（genes的数量），分母是所有输入的做富集分析的gene的数目（dotplot的横轴）
    
    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "03.Analysis"

    def enrich_single(self, file_single: str, result_file: str, logger: logging.Logger):
        logger.info("Enrich {} running ...".format(os.path.basename(file_single)[:-4]))
        
        cmd = [
            "/usr/bin/Rscript-4.2",
            os.path.join(os.path.dirname(__file__), "jiazi_enrich.r"),
            "-o",
            result_file,
            file_single
        ]
        self.exec_cmd(" ".join(cmd), logger)

        logger.info("Enrich {} complete !".format(os.path.basename(file_single)[:-4]))

    def run(
            self,
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="GeneSymbols",
                    data_description="需要富集的基因，文件内每一行为一个基因，要求为基因symbol。多类基因放在不同文件内，使用英文逗号间隔",
                    is_multi=True
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']
        file_genes = tuple(filter(
            lambda d: os.path.exists(d),
            kwargs['return'][0].split(",")
        ))
        dir_result = os.getcwd()

        func_list = list()
        args_list = list()
        os.makedirs(dir_result, exist_ok=True)
        for file_single in file_genes:
            func_list.append(self.enrich_single)
            args_list.append(dict(
                file_single=file_single,
                result_file=os.path.join(dir_result, os.path.basename(file_single)[:-4]),
                logger=logger
            ))
        self.multi_task(func_list, args_list, logger)
        
        