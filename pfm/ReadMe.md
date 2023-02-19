##### 🙈: yukaiquan
##### 📧: 1962568272@qq.com
# 1 软件安装
```bash
conda install -c bioconda hmmer
```
# 2 全基因组hmm搜索
Pfam-A.hmm为最新的pfm35数据库
```bash
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/atlantica.fasta.output Pfam-A.hmm 01_pep/atlantica.fasta >01_pep/atlantica.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/eriantha.fasta.output Pfam-A.hmm 01_pep/eriantha.fasta >01_pep/eriantha.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/OTv2.fasta.output Pfam-A.hmm 01_pep/OTv2.fasta >01_pep/OTv2.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/Sang.fasta.output Pfam-A.hmm 01_pep/Sang.fasta >01_pep/Sang.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/Sang_ins.fasta.output Pfam-A.hmm 01_pep/Sang_ins.fasta >01_pep/Sang_ins.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/Sang_lon.fasta.output Pfam-A.hmm 01_pep/Sang_lon.fasta >01_pep/Sang_lon.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/SAU_ins.fasta.output Pfam-A.hmm 01_pep/SAU_ins.fasta >01_pep/SAU_ins.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/SAU_lon.fasta.output Pfam-A.hmm 01_pep/SAU_lon.fasta >01_pep/SAU_lon.fasta.log 2>&1
hmmsearch --cut_tc --cpu 16 --domtblout 01_pep/strigosa.fasta.output Pfam-A.hmm 01_pep/strigosa.fasta >01_pep/strigosa.fasta.log 2>&1


nohup sh work.sh &
```
# 3 数据整理

