from Tasks import taskTemplate, dataTemplate, returnTemplate, argIntOrFloatTemplate, argSelectTemplate, get_pipline_data, get_pipline_species
from itertools import product
import logging
import os

class Task(taskTemplate):
    """ GetRegionGenes

    根据region获取gene

    ### 一、输入参数说明
    
    
    ### 二、输出文件说明


    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "04.Tool"

    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            ratio: argIntOrFloatTemplate(
                arg_name="百分比",
                arg_description="region覆盖超过promoter/genebody一定比例该基因才算，0则为只要覆盖即算（至少1bp）",
                arg_default=0,
                arg_min=0,
                arg_max=100,
                arg_step=1
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_beds",
                    data_description="计算的region文件",
                    is_multi=True
                ), # type: ignore
            ),
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']
        dir_base = os.getcwd()

        data_species = get_pipline_data(species)
        file_genome = {
            "Promoter": data_species["promoter_file"],
            "GeneBody": data_species["genebody_file"],
        }

        file_beds = kwargs['return'][0].strip(",").split(",")

        ratio_str = "-f {}".format(ratio) if ratio > 0 else ""

        for g, file_bed in product(file_genome.keys(), file_beds):
            gene_name = os.path.basename(file_bed)[:-4]
            file_res = os.path.join(dir_base, "{}.{}.txt".format(gene_name, g))

            logger.info("running {} {} ...".format(gene_name, g))

            assert os.path.exists(file_genome[g]), "Genome {} bed not exists".format(file_genome[g])
            assert os.path.exists(file_bed), "file bed not exist: {}".format(file_bed)
            
            self.exec_cmd("intersectBed -wa {} -a {} -b {} | cut -f 5 | sort | uniq > {}".format(
                ratio_str,
                file_genome[g],
                file_bed,
                file_res
            ), logger)