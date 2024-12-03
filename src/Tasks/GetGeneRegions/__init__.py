from Tasks import taskTemplate, dataTemplate, returnTemplate, argSelectTemplate, get_pipline_data, get_pipline_species
from itertools import product
import pandas as pd
import logging
import os

class Task(taskTemplate):
    """ GetGeneRegions

    获取基因对应的region：promoter/genebody

    ### 一、输入参数说明
    
    
    ### 二、输出文件说明


    ### 三、任务作者

    Dong LIU
    """
    def init(self) -> None:
        self.label = "04.Tool"


    def read_genes(self, file_gene: str):
        with open(file_gene) as f:
            genes_temp = set(filter(
                lambda d: d != '',
                map(
                    lambda dd: dd.strip(),
                    f.readlines()
                )
            ))
        return genes_temp



    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="file_genes",
                    data_description="计算的基因列表",
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
        for i in ("corepromoter", "DistalRegion", "TSS"):
            file_temp = data_species["promoter_file"].replace("promoter", i)
            if not os.path.exists(file_temp):
                continue
            file_genome[i] = file_temp

        file_res_all = set()
        for g, file_gene in product(file_genome.keys(), kwargs['return'][0].strip(",").split(",")):
            gene_name = os.path.basename(file_gene)[:-4]
            file_res = os.path.join(dir_base, "{}.{}.bed".format(gene_name, g))

            logger.info("running {} {} ...".format(gene_name, g))

            assert file_res not in file_res_all, "Gene {} duplicate".format(gene_name)
            assert os.path.exists(file_genome[g]), "Genome {} bed not exists".format(file_genome[g])
            assert os.path.exists(file_gene), "file gene not exist: {}".format(file_gene)
            
            file_res_all.add(file_res)
            gene_need = self.read_genes(file_gene)

            df = pd.read_table(file_genome[g], header=None)
            df.columns = ["chr", "start", "end", "trans", "gene_symbol", "strand"]
            df = df[df["gene_symbol"].isin(gene_need)]
            df[["chr", "start", "end", "strand", "gene_symbol"]].drop_duplicates().to_csv(
                file_res,
                index=False,
                header=False,
                sep="\t"
            )
            del df
            del gene_need