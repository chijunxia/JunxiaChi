# 论文

---
## 处理流程
```
目标简述：
    通过寻找Mix中差异表达的microRNA,从两个方面进行研究：
        1.植物抵抗(植物microRNA进入真菌体内影响真菌的致病功能的进行[需要选择能匹配到植物的])
        # 2.真菌靶向(真菌microRNA进入植物体内使植物致病) --->由于暂时没有发现真菌microRNA数据集搁置
    Start：植物抵抗：
        a.从c707_18to25.fastq中筛选能匹配到植物基因组的部分reads -> Mix_mapPlants.fastq
        b.保留MixMapPlants.fastq中相比较c697_18to25.fastq明显增强和数量比较多的reads —> Mix_mapPlants_upAmount.fastq
        c.这里认为Mix_mapPlants_upAmount.fastq为假设起作用[注一]的植物microRNA, 分别靶基因预测[大豆和大豆疫霉]
        d.进行功能富集分析
        e.模型构建
       
    注一：植物microRNA起两方面作用：
        1.真菌侵染植物后，真菌的起作用的microRNA导致植物应激反应，使植物某些microRNA大量增加来[被动]反抗真菌侵染
        2.真菌侵染植物后，植物某些microRNA通过侵入真菌[主动]反抗真菌侵染
```
---

## Done
+ 大豆和大豆疫霉数据
+ 找到接头序列adapter
+ 数据预处理工作
+ 完成序列操作(匹配-去冗余等)
+ 标准化
+ 靶基因预测(RnaHybrid)
+ 靶基因预测(Tapir)
+ 火山图
+ 计算选定集合的集合(ttest和pvalue)

## ToDo
+ 总体流程图及各章节流程图


## 资源索引
数据来源论文：https://www.researchgate.net/publication/272376050_Coordination_of_MicroRNAs_PhasiRNAs_and_NB-LRR_Genes_in_Response_to_a_Plant_Pathogen_Insights_from_Analyses_of_a_Set_of_Soybean_Rps_Gene_Near-Isogenic_Lines
Tapir命令:tapir_hybrid sszInUseq.fasta sequence.fasta | hybrid_parser > sszTapir1.txt
大豆基因组数据库SoyBean：https://www.soybase.org/dlpages/
大豆疫霉基因组：http://fungi.ensembl.org/
samtools使用方法：https://www.cnblogs.com/emanlee/p/4316581.html






成果
1．	分析了大豆被大豆侵染后差异表达的大豆microRNA的特点
2．	基于1的特点，构建了机器学习预测模型，用于预测未知microRNA与大豆抵抗疫霉侵染作用是否相关以及相关程度
3．	通过将大豆microRNA靶向大豆疫霉的mRNA，进行功能富集分析，分析大豆对于大豆疫霉反向侵染的可能性，为大豆和大豆疫霉的跨界调控研究提供数据依据
4.	为植物抵抗真菌致病性侵染提供了通用的数据分析方案