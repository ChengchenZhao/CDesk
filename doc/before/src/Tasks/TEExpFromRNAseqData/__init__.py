from Tasks import taskTemplate, dataTemplate, returnTemplate, argSelectTemplate, get_pipline_data, get_pipline_species, get_pipline_samples
import pandas as pd
import logging
import os

class Task(taskTemplate):
    """ TEExpFromRNAseqData

    TransposableElementFromRNAseq的pipeline

    ### 一、输入参数说明
    
    #### 1. 数据：`DataPath`
    
    测序文件(fq、fq.gz、fastq或fastq.gz)所在的文件夹，可以根据输入框右边的辅助按钮选择文件夹
    
    该文件夹内的所有测序文件都会被运行。

    #### 2. 参数：物种
    
    该测序文件对应的物种，选择即可。
    
    ### 二、输出文件说明

    ### 三、任务作者

    Dong LIU、Chengchen ZHAO
    """
    def init(self) -> None:
        self.label = "02.PreProcess"

    def calc_gene_length(self, file_bed: str) -> dict:
        df = pd.read_table(file_bed, header=None)
        cols = list(df.columns)
        cols_need = ["chr", "start", "end", "strand", "gene"]
        df.columns = cols_need + cols[5:]
        df = df[cols_need]
        df['len'] = df["end"].astype(int) - df["start"].astype(int) + 1
        return df.groupby("gene").agg({"len": sum}).to_dict()["len"]
        

    def calc_fpkm_single(self, df: pd.DataFrame, length_gene: dict):
        sample_all = list(df.columns)
        len_se = df.index.map(lambda d: length_gene[d])

        for i in sample_all:
            sum_temp = df[i].sum()
            df[i] = df[i] / len_se * (10 ** 9) / sum_temp

        return df
        

    def calc_fpkm(self, file_exp: str, file_res: str, file_annote_gene: str, file_annote_te: str, is_split: bool = True):
        if os.path.exists(file_res):
            return

        df = pd.read_table(file_exp, index_col=0)

        length_gene = self.calc_gene_length(file_annote_gene)
        length_te = self.calc_gene_length(file_annote_te)

        if is_split:
            df_gene = df[df.index.isin(set(length_gene.keys()))]
            df_te = df[df.index.isin(set(length_te.keys()))]
            df_gene.to_csv(file_res.replace(".fpkm.split", ".counts").replace(".txt", ".Gene.txt"), sep="\t")
            df_te.to_csv(file_res.replace(".fpkm.split", ".counts").replace(".txt", ".TE.txt"), sep="\t")


            df_gene = self.calc_fpkm_single(df_gene, length_gene)
            df_te = self.calc_fpkm_single(df_te, length_te)
            
            df_gene.to_csv(file_res.replace(".split.txt", ".Gene.txt"), sep="\t")
            df_te.to_csv(file_res.replace(".split.txt", ".TE.txt"), sep="\t")
            df = pd.concat([df_gene, df_te], axis=0)
        else:
            length_gene.update(length_te)
            df = df[df.index.isin(set(length_gene.keys()))]
            df = self.calc_fpkm_single(df, length_gene)
        
        df.to_csv(file_res, sep="\t")

    def te_csv_to_txt(self, file_te: str, logger: logging.Logger, file_annote_gene: str, file_annote_te: str):
        file_res = file_te.replace(".csv", ".txt")
        df = pd.read_csv(file_te, index_col=0).T
        df.index.name = "gene_symbol"
        df.to_csv(
            file_res,
            sep="\t"
        )
        self.exec_cmd("rm -rf {}".format(file_te), logger)

        for i in (True, False):
            self.calc_fpkm(
                file_exp=file_res,
                file_res=file_res.replace("counts", "fpkm{}".format(".split" if i else "")),
                file_annote_gene=file_annote_gene,
                file_annote_te=file_annote_te,
                is_split=i
            )

    def run(
            self,
            species: argSelectTemplate(
                arg_name="物种",
                arg_description="数据的物种",
                arg_default="mouse",
                options=get_pipline_species()
            ), # type: ignore
            data_type: argSelectTemplate(
                arg_name="数据类型",
                arg_description="输入数据类型",
                arg_default="fastq",
                options={
                    "fastq": "Fastq",
                    "bam": "Bam",
                    "scbam": "scBam"
                }
            ), # type: ignore
            **kwargs
        ) -> returnTemplate(
            data_in=(
                dataTemplate(
                    data_name="DataPath",
                    data_description="数据文件夹"
                ),
            ), # type: ignore
            data_out=tuple()
        ):
        logger: logging.Logger = kwargs['logger']

        assert os.path.exists(kwargs['return'][0]), "data path not exists: {}".format(kwargs['return'][0])

        data_species = get_pipline_data(species)

        if data_type == "fastq":
            samples_dict = get_pipline_samples(kwargs['return'][0])

            for k in samples_dict:
                cmd = [
                    os.path.join(os.path.dirname(__file__), "Mapping.sh"),
                    samples_dict[k]['data_path'],
                    os.getcwd(),
                    k,
                    samples_dict[k]['type'],
                    ",".join(map(
                        lambda d: os.path.basename(d),
                        samples_dict[k]["data"]
                    )),
                    data_species["mapping_index"],
                    data_species["species"],
                    data_species["promoter_file"],
                    data_species["chromInfo"]
                ]
                self.exec_cmd(" ".join(map(
                    lambda d: "\"{}\"".format(d),
                    cmd
                )), logger)
        elif data_type == "bam":
            self.exec_cmd("ln -s {} {}".format(
                kwargs['return'][0],
                os.path.join(os.getcwd(), "bam")
            ), logger)
        elif data_type == "scbam":
            dir_data = kwargs['return'][0]
            for root, _, files in os.walk(dir_data):
                for file_single in files:
                    if not file_single.endswith(".bam"):
                        continue
                    file_bam_temp = os.path.join(dir_data, root, file_single)
                    name_bam_temp = "{}.{}".format(
                        os.path.basename(os.path.dirname(file_bam_temp)),
                        file_single[:-4]
                    )
                    self.exec_cmd(" ".join([
                        "scTE",
                        "-i",
                        file_bam_temp,
                        "-o",
                        "TE.counts.{}".format(name_bam_temp),
                        "-x",
                        data_species["TE_idx"],
                        "-p",
                        str(os.environ.get("cdesk_task_num", 8)),
                        "--hdf5 False"
                    ]), logger)
                    self.te_csv_to_txt(
                        file_te=os.path.join(os.getcwd(), "TE.counts.{}.csv".format(name_bam_temp)),
                        file_annote_gene=data_species["TE_idx_bed_gene"],
                        file_annote_te=data_species["TE_idx_bed_te"],
                    )
        else:
            logger.error("data type error !")
            
            return


        dir_bam = os.path.join(os.getcwd(), "bam")
        file_bam_all = list(filter(
            lambda d: d.endswith(".bam") and (os.path.getsize(d) > 1000),
            map(
                lambda dd: os.path.join(dir_bam, dd),
                os.listdir(dir_bam)
            )
        ))
        if len(file_bam_all) == 0:
            logger.error("No valid bam file !")
            return
        self.exec_cmd(" ".join([
            "scTE",
            "-i",
            ",".join(file_bam_all),
            "-o",
            "TE.counts",
            "-x",
            data_species["TE_idx"],
            "-p",
            str(os.environ.get("cdesk_task_num", 8)),
            "--hdf5 False -CB False -UMI False"
        ]), logger)
        
        self.te_csv_to_txt(
            file_te=os.path.join(os.getcwd(), "TE.counts.csv"),
            file_annote_gene=data_species["TE_idx_bed_gene"],
            file_annote_te=data_species["TE_idx_bed_te"],
        )


