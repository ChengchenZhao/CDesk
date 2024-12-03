'''
@File    :   pdf_combine.r
@Modif   :   2024/01/21
@Author  :   winter
@Version :   0.1
@Contact :   winter_lonely@126.com
@License :   (C)Copyright winter
@Desc    :   
'''

library(qpdf)
sys_argv <- commandArgs(T)
file_path <- sys_argv[1]

file_list <- list.files(
    path=file_path,
    pattern=".*pdf",
    full.names=T
)

pdf_combine(file_list, output=paste(file_path,"combine.pdf",sep="/"))