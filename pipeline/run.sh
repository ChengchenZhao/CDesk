bash NewBashFile.sh RunRNAseq
bash NewBashFile.sh RunATACseq
bash NewBashFile.sh RunChIPseq
bash NewBashFile.sh RunCUTRUN
bash NewBashFile.sh RunRRBS

su
cp /mnt/zhaochengchen/bin/pipeline/RNA_seq_preprocessing.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/RemoveDuplicateGenes.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/ReadsDistribution.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/read_distribution.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/pipeline/RunRNAseq /usr/local/bin/
chmod +x /usr/local/bin/RunRNAseq
chmod +x /usr/local/bin/RNA_seq_preprocessing.py
chmod +x /usr/local/bin/RemoveDuplicateGenes.py
chmod +x /usr/local/bin/read_distribution.py

cp /mnt/zhaochengchen/bin//GetBigwigSignal.py /usr/local/bin/
cp /mnt/zhaochengchen/bin//GeneSymbolUniq.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/pipeline/RunATACseq /usr/local/bin/
chmod +x /usr/local/bin/RunATACseq
chmod +x /usr/local/bin/GetBigwigSignal.py
chmod +x /usr/local/bin/GeneSymbolUniq.py


cp /mnt/zhaochengchen/bin/pipeline/RunCUTRUN /usr/local/bin/
chmod +x /usr/local/bin/RunCUTRUN

cp /mnt/zhaochengchen/bin/pipeline/RunRRBS /usr/local/bin/
chmod +x /usr/local/bin/RunRRBS

cp /mnt/zhaochengchen/bin/DataDownloadFromGEO_wget.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/DataDownloadFromGEO_wget_new.py /usr/local/bin/
cp /mnt/zhaochengchen/bin/DataDownloadFromGEO_prefetch.py /usr/local/bin/
chmod +x /usr/local/bin/DataDownloadFromGEO_wget.py
chmod +x /usr/local/bin/DataDownloadFromGEO_wget_new.py
chmod +x /usr/local/bin/DataDownloadFromGEO_prefetch.py

cp /mnt/zhaochengchen/bin/pipeline/RunChIPseq /usr/local/bin/
chmod +x /usr/local/bin/RunChIPseq

cp /mnt/zhaochengchen/bin/pipeline/RunChIPseqWOCtrl /usr/local/bin/
chmod +x /usr/local/bin/RunChIPseqWOCtrl



