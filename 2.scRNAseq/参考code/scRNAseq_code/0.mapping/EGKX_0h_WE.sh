###################################################################################################
################################           1.MERGE FILES           ################################
###################################################################################################
bam_path=/data/zhaochengchen/Work/2.ChickEmbryo/results/0.mapping/bam
log_path=/data/zhaochengchen/Work/2.ChickEmbryo/results/0.mapping/log
matrix_path=/data/zhaochengchen/Work/2.ChickEmbryo/results/0.mapping/matrix
plot_path=/data/zhaochengchen/Work/2.ChickEmbryo/results/0.mapping/plot

mkdir ${bam_path}
mkdir ${log_path}
mkdir ${matrix_path}
mkdir ${plot_path}

result_path=/data/huangxingnan/Data/scRNA/1.ChickenEmbryo/script/batch1
samtools merge -@ 10 -o EGKX_0h_WE_c1.sam ${result_path}/EGKX_0h_WE_r1/mapping/EGKX_0h_WE_r1.sam ${result_path}/EGKX_0h_WE_r2/mapping/EGKX_0h_WE_r2.sam ${result_path}/EGKX_0h_WE_r3/mapping/EGKX_0h_WE_r3.sam ${result_path}/EGKX_0h_WE_r4/mapping/EGKX_0h_WE_r4.sam
cat ${result_path}/EGKX_0h_WE_r1_1.fastq ${result_path}/EGKX_0h_WE_r2_1.fastq ${result_path}/EGKX_0h_WE_r3_1.fastq ${result_path}/EGKX_0h_WE_r4_1.fastq > EGKX_0h_WE_c1_1.fastq
cat ${result_path}/EGKX_0h_WE_r1_2.fastq ${result_path}/EGKX_0h_WE_r2_2.fastq ${result_path}/EGKX_0h_WE_r3_2.fastq ${result_path}/EGKX_0h_WE_r4_2.fastq > EGKX_0h_WE_c1_2.fastq
result_path2=/data/huangxingnan/Data/scRNA/1.ChickenEmbryo/script/batch2
samtools merge -@ 10 -o EGKX_0h_WE_c2.sam ${result_path}/EGKX_0h_WE_r5/mapping/EGKX_0h_WE_r5.sam ${result_path2}/EGKX_0h_WE_r6/mapping/EGKX_0h_WE_r6.sam ${result_path2}/EGKX_0h_WE_r7/mapping/EGKX_0h_WE_r7.sam ${result_path2}/EGKX_0h_WE_r8/mapping/EGKX_0h_WE_r8.sam
cat ${result_path}/EGKX_0h_WE_r5_1.fastq ${result_path2}/EGKX_0h_WE_r6_1.fastq ${result_path2}/EGKX_0h_WE_r7_1.fastq ${result_path2}/EGKX_0h_WE_r8_1.fastq > EGKX_0h_WE_c2_1.fastq
cat ${result_path}/EGKX_0h_WE_r5_2.fastq ${result_path2}/EGKX_0h_WE_r6_2.fastq ${result_path2}/EGKX_0h_WE_r7_2.fastq ${result_path2}/EGKX_0h_WE_r8_2.fastq > EGKX_0h_WE_c2_2.fastq
samtools merge -@ 10 -o EGKX_0h_WE_c3.sam ${result_path2}/EGKX_0h_WE_r9/mapping/EGKX_0h_WE_r9.sam ${result_path2}/EGKX_0h_WE_r10/mapping/EGKX_0h_WE_r10.sam ${result_path2}/EGKX_0h_WE_r11/mapping/EGKX_0h_WE_r11.sam ${result_path2}/EGKX_0h_WE_r12/mapping/EGKX_0h_WE_r12.sam
cat ${result_path2}/EGKX_0h_WE_r9_1.fastq ${result_path2}/EGKX_0h_WE_r10_1.fastq ${result_path2}/EGKX_0h_WE_r11_1.fastq ${result_path2}/EGKX_0h_WE_r12_1.fastq > EGKX_0h_WE_c3_1.fastq
cat ${result_path2}/EGKX_0h_WE_r9_2.fastq ${result_path2}/EGKX_0h_WE_r10_2.fastq ${result_path2}/EGKX_0h_WE_r11_2.fastq ${result_path2}/EGKX_0h_WE_r12_2.fastq > EGKX_0h_WE_c3_2.fastq
samtools merge -@ 10 -o EGKX_0h_WE_c4.sam ${result_path2}/EGKX_0h_WE_r13/mapping/EGKX_0h_WE_r13.sam ${result_path2}/EGKX_0h_WE_r14/mapping/EGKX_0h_WE_r14.sam ${result_path2}/EGKX_0h_WE_r15/mapping/EGKX_0h_WE_r15.sam ${result_path2}/EGKX_0h_WE_r16/mapping/EGKX_0h_WE_r16.sam
cat ${result_path2}/EGKX_0h_WE_r13_1.fastq ${result_path2}/EGKX_0h_WE_r14_1.fastq ${result_path2}/EGKX_0h_WE_r15_1.fastq ${result_path2}/EGKX_0h_WE_r16_1.fastq > EGKX_0h_WE_c4_1.fastq
cat ${result_path2}/EGKX_0h_WE_r13_2.fastq ${result_path2}/EGKX_0h_WE_r14_2.fastq ${result_path2}/EGKX_0h_WE_r15_2.fastq ${result_path2}/EGKX_0h_WE_r16_2.fastq > EGKX_0h_WE_c4_2.fastq
samtools merge -@ 10 -o EGKX_0h_WE_c5.sam ${result_path2}/EGKX_0h_WE_r17/mapping/EGKX_0h_WE_r17.sam ${result_path2}/EGKX_0h_WE_r18/mapping/EGKX_0h_WE_r18.sam ${result_path2}/EGKX_0h_WE_r19/mapping/EGKX_0h_WE_r19.sam ${result_path2}/EGKX_0h_WE_r20/mapping/EGKX_0h_WE_r20.sam
cat ${result_path2}/EGKX_0h_WE_r17_1.fastq ${result_path2}/EGKX_0h_WE_r18_1.fastq ${result_path2}/EGKX_0h_WE_r19_1.fastq ${result_path2}/EGKX_0h_WE_r20_1.fastq > EGKX_0h_WE_c5_1.fastq
cat ${result_path2}/EGKX_0h_WE_r17_2.fastq ${result_path2}/EGKX_0h_WE_r18_2.fastq ${result_path2}/EGKX_0h_WE_r19_2.fastq ${result_path2}/EGKX_0h_WE_r20_2.fastq > EGKX_0h_WE_c5_2.fastq
result_path3=/data/huangxingnan/Data/scRNA/1.ChickenEmbryo/script/batch3
samtools merge -@ 10 -o EGKX_0h_WE_c6.sam ${result_path3}/EGKX_0h_WE_r21/mapping/EGKX_0h_WE_r21.sam ${result_path3}/EGKX_0h_WE_r22/mapping/EGKX_0h_WE_r22.sam ${result_path3}/EGKX_0h_WE_r23/mapping/EGKX_0h_WE_r23.sam ${result_path3}/EGKX_0h_WE_r24/mapping/EGKX_0h_WE_r24.sam
cat ${result_path3}/EGKX_0h_WE_r21_1.fastq ${result_path3}/EGKX_0h_WE_r22_1.fastq ${result_path3}/EGKX_0h_WE_r23_1.fastq ${result_path3}/EGKX_0h_WE_r24_1.fastq > EGKX_0h_WE_c6_1.fastq
cat ${result_path3}/EGKX_0h_WE_r21_2.fastq ${result_path3}/EGKX_0h_WE_r22_2.fastq ${result_path3}/EGKX_0h_WE_r23_2.fastq ${result_path3}/EGKX_0h_WE_r24_2.fastq > EGKX_0h_WE_c6_2.fastq

DrSeq run -c EGKX_0h_WE_c1.conf -f --clean
DrSeq run -c EGKX_0h_WE_c2.conf -f --clean
DrSeq run -c EGKX_0h_WE_c3.conf -f --clean
DrSeq run -c EGKX_0h_WE_c4.conf -f --clean
DrSeq run -c EGKX_0h_WE_c5.conf -f --clean
DrSeq run -c EGKX_0h_WE_c6.conf -f --clean

for sample in EGKX_0h_WE_c1 EGKX_0h_WE_c2 EGKX_0h_WE_c3 EGKX_0h_WE_c4 EGKX_0h_WE_c5 EGKX_0h_WE_c6
do
    samtools view -@ 10 -Sb ${sample}.sam > ${bam_path}/${sample}.bam
    cp ${sample}/progress_log.txt ${log_path}/${sample}.progress_log.txt
    cp ${sample}/summary/plots/${sample}_map_summary.txt ${log_path}/${sample}.map_summary.txt
    
    cp ${sample}/expmatrix/${sample}_qcmat.txt ${matrix_path}/
    cp ${sample}/expmatrix/${sample}_expmat.txt ${matrix_path}/

    cp ${sample}/summary/plots/${sample}*.pdf ${plot_path}/
    cp ${sample}/summary/${sample}_summary.pdf ${plot_path}/

    gzip ${sample}_1.fastq
    gzip ${sample}_2.fastq
    # rm ${sample}.sam
done