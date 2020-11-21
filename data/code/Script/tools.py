#coding=utf-8
import logging
import csv
from utils import ReadFile
logging.getLogger().setLevel(logging.INFO)

rd = ReadFile()


def cal_GC():
    """计算GC含量"""
    data = open('sszAll.txt', 'r').readlines()
    count = 1
    tcount = 0
    ls = []
    for key in data:
        tcount = 0
        if count%2 == 0:
            for i in key.strip():
                if i=='G' or i=='C':
                    tcount = tcount+1
            ls.append([key.strip(), tcount/len(key.strip())])
        count = count+1
    newls = sorted(ls, key=lambda k: k[1])
    print(newls)


def cal_18T25():
    """统计长度位于18-25之间的数据信息"""
    f = open(r'../../../DATA/FASTQ/c697.fastq', 'r')
    fo = open(r'../../../DATA/FASTQ/c697_1201.fastq', 'a+')
    lenc = 1
    while lenc < 35313924:
        logging.info(lenc)
        l1 = f.readline()
        l2 = f.readline()
        l3 = f.readline()
        l4 = f.readline()
        if len(l2.strip()) <= 25 and len(l2.strip()) >=18:
            fo.writelines(l1)
            fo.writelines(l2)
            fo.writelines(l3)
            fo.writelines(l4)
        lenc = lenc + 4
    f.close()
    fo.close()

    f1 = open(r'../../../DATA/FASTQ/c707.fastq', 'r')
    fo1 = open(r'../../../DATA/FASTQ/c707_1201.fastq', 'a+')
    lenc1 = 1
    while lenc1 < 89623084:
        logging.info(lenc1)
        l1 = f1.readline()
        l2 = f1.readline()
        l3 = f1.readline()
        l4 = f1.readline()
        if len(l2.strip()) <= 25 and len(l2.strip()) >= 18:
            print(len(l2.strip()))
            print(l2)
            fo1.writelines(l1)
            fo1.writelines(l2)
            fo1.writelines(l3)
            fo1.writelines(l4)
        lenc1 = lenc1 + 4
    f1.close()
    fo1.close()


def findOne():
    """确实是否有异常数据"""
    f = open(r'../../../DATA/FASTQ/c707_1117.fastq', 'r')
    lenc = 1
    while lenc <= 24295716:
        l1 = f.readline()
        l2 = f.readline()
        l3 = f.readline()
        l4 = f.readline()
        if len(l2.strip()) != len(l4.strip()):
            logging.info(lenc)
        lenc = lenc + 4
    logging.info("END")


def Aliases():
    """根据Aliases将XM改为Gly格式"""
    data = open('Aliases/Soybean.3847.protein.aliases.v11.0.txt', 'r').readlines()
    gene_result = []
    xm_result = []
    union_result = []
    for key in data:
        if key.strip().split("\t")[len(key.strip().split("\t"))-1].split(" ")[0] == "BLAST_UniProt_DR_GeneID" or key.strip().split("\t")[len(key.strip().split("\t"))-1].split(" ")[0] == "Ensembl_UniProt_DR_GeneID":
            gene_result.append([key.strip().split("\t")[0], key.strip().split("\t")[1]])
        if key.strip().split("\t")[len(key.strip().split("\t"))-2].split("_")[0] == "XM" and key.strip().split("\t")[len(key.strip().split("\t"))-1] == "Ensembl_RefSeq_short":
            xm_result.append([key.strip().split("\t")[0], key.strip().split("\t")[1]])
    print(len(gene_result))
    print(len(xm_result))

    for g in gene_result:
        for x in xm_result:
            if x[0] == g[0]:
                union_result.append([g[0], [g[1], x[1]]])
    print(len(union_result))

    xm = open('Target/TargetSoybeanXM.txt', 'r').readlines()
    xm_target = []
    for t in xm:
        xm_target.append(t.strip())

    for u in union_result:
        if u[1][1] in xm_target:
            print(u[1][0])


def T1to2():
    data = open('Target/G1T2.csv', 'r').readlines()
    for d in data:
        print(d.strip().split(" ")[len(d.strip().split(" "))-1].strip())

def AliasSojae():
    data = open('Target/TargetSojae.txt', 'r').readlines()
    raw = []
    for k in data:
        raw.append(k.strip())

    AliasData = open('Aliases/Sojae.67593.protein.aliases.v11.0.txt', 'r')
    for ad in AliasData:
        if ad.strip().split("\t")[len(ad.strip().split("\t"))-1] == "Ensembl_protein_id":
            if ad.strip().split("\t")[0].split(".")[1] in raw:
                print(ad.strip().split("\t")[0].split(".")[1])


def get_U():
    """统计"""
    A = []
    B = []
    C = []
    Adata = open('Mapping/NewC707_18T25.fastq', 'r').readlines()
    Acount = 0
    for ad in Adata:
        Acount = Acount + 1
        if Acount%4 == 2:
            A.append(ad.strip())
    logging.info(len(A))

    Bdata = open('Mapping/NewC697_18T25.fastq', 'r').readlines()
    Bcount = 0
    for bd in Bdata:
        Bcount = Bcount + 1
        if Bcount%4 == 2:
            B.append(bd.strip())
    logging.info(len(B))

    dA = rd.countDownD(A)
    logging.info(len(dA))
    dB = rd.countDownD(B)
    logging.info(len(dB))

    rootIn = "Mapping/New_c7_c6.csv"
    rootNo = "Mapping/New_C1C2FcPvalue_No_200121.csv"

    I = []
    Indata = open('Mapping/707_In_Soybean_Not_Sojae.fastq', 'r').readlines()
    for inD in Indata:
        I.append(inD.strip())
    logging.info(len(I))
    dI = rd.countDownD(I)
    logging.info(len(dI))

    rd.countDownAll(dB, dA, rootIn, rootNo, dI)


def cut_value():
    csvData = open('Mapping/New707In697_count_Rate_1201.csv', 'r')
    opt_sample = open('Mapping/NavSample_seq.fasta', 'a+')
    reader = csv.reader(csvData)
    pre_title = ">NewSSZNav-"
    count = 0
    Allcount = 0
    Acount = 0
    Qcount = 0
    Rcount = 0
    for item in reader:
        if reader.line_num == 1:
            continue
        else:
            if float(item[3]) == 0:
                Allcount += 1
                if Allcount%102 == 0:
                    Acount += 1
                    #opt_sample.write(pre_title+str(Acount)+"\n")
                    opt_sample.write(str(item[0]).replace('T', 'U')+"\n")
                    #print(item)
                    count += 1

    logging.info("U: %d", count)

    csvData2 = open('Mapping/New707No697_count_Rate_1201.csv', 'r')
    reader2 = csv.reader(csvData2)
    count2 = 0
    for item in reader2:
        if reader2.line_num == 1:
            continue
        else:
            if float(item[2]) > 100:
                Allcount += 1
                #opt_sample.write(pre_title + str(Allcount) + "\n")
                #opt_sample.write(str(item[0]).replace('T', 'U') + "\n")
                #print(item)
                count2 += 1
    print(count2)
    opt_sample.close()


def write_In():
    data = open("Mapping/New707In697_count_Rate.csv").readlines()
    AllSeq = []
    count = 0
    for d in data:
        count += 1
        if count == 1:
            continue
        if float(d.strip().split(",")[2]) > 200 and float(d.strip().split(",")[3]) >10:
            AllSeq.append(d.strip().split(",")[0])

    outfile = open('New_732.csv', 'w')
    Twodata = open("Mapping/New_c7_c6.csv").readlines()
    for t in Twodata:
        if t.strip().split(",")[0].strip() in AllSeq:
            print(t.strip())
            outfile.write(t)
    outfile.close()


def static_In():
    import numpy as np
    # data = open("Mapping/New_732.csv").readlines()
    # l_707, l_697 = [], []
    # for d in data:
    #     l_707.append(int(d.strip().split(",")[1].strip()))
    #     l_697.append(int(d.strip().split(",")[2].strip()))

    import pandas as pd
    import matplotlib.pyplot as plt
    data_all = pd.read_csv("Mapping/New_732.csv")
    print(data_all.shape)
    # data_all.plot(kind='bar', title='1')
    # data_all.plot(kind='box', title='2')
    data_all['dif']=data_all['c7']-data_all['c6']
    import seaborn as sns
    sns.distplot(data_all['dif'])
    plt.title('Difference distribution')
    plt.xlabel("Difference between two sets of data")
    plt.show()

    from scipy import stats
    t, p_2tail = stats.ttest_rel(data_all['c7'], data_all['c6'])
    p_1tail = p_2tail/2
    print(t)
    print(p_1tail)


def get_soybean_protein_id():
    count = 0
    path = "TargetSoybean/NewSoybeanOpr1_768.txt"
    data = open(path).readlines()
    SoybeanList = []
    SoybeanLists = []
    SoybeanDict = {}
    for key in data:
        count += 1
        if (count-5)%16 == 2:
            temp = key.strip().split(".")[-2]+"."+key.strip().split(".")[-1].strip()
            SoybeanDict[temp] = 1 if not SoybeanDict.__contains__(temp) else SoybeanDict[temp]+1
    SoybeanDict = dict(sorted(SoybeanDict.items(), key=lambda item:-item[1]))
    print(len(SoybeanDict))

    for i in SoybeanDict.keys():
        SoybeanList.append(i)
    #SoybeanList = list(set(SoybeanList))
    print(len(SoybeanList))
    count = 0
    for k in SoybeanList:
        # print(k.split(".")[0])
        count += 1
        if count >= 2000 :
            print(k.split(".")[0])
        else:
            continue

def get_sojae_protein_id():
    Tlines = open("TargetSojae/TSojae1_768.txt").readlines()
    tlines = []
    tcount = 0
    for t in Tlines:
        tcount += 1
        if (tcount-5)%16 == 2:
            tlines.append(t.strip().split("         ")[1].strip())
    print(len(tlines))
    tlines = list(set(tlines))
    print(len(tlines))

    Mlines = open("TargetSojae/Sojae1214.fasta").readlines()
    mlines = []
    for m in Mlines:
        if m[0] == ">":
            if m.strip().split(">")[1].split(" ")[0].strip() in tlines:
                mlines.append(m.strip().split("locus_tag=")[1].split("]")[0].strip() if len(m.strip().split("locus_tag="))>1 else None)
    print(len(mlines))
    Adata = open("Aliases/Sojae.67593.protein.aliases.v11.0.txt")
    alines = []
    for a in Adata:
        if len(a.strip().split("PHYSODRAFT_"))>1:
            if a.strip().split("psoj:")[1].split("	")[0] in mlines:
                alines.append(a.strip().split("	")[0].split(".")[1])
                # print(a.strip().split("psoj:")[1].split("	")[0])
                # print(a.strip().split("	")[0].split(".")[1])
    print(len(alines))

    count = 0
    for i in alines:
        if count < 100:
            count += 1
            print(i)


def get_protein_seq():
    count = 0
    iden = 0
    str = "ATGATAGGAAGAGCCGACATCGAAGGATCAAAAAGCAACGTCGCTATGAACGCTTGGCTGCCACAAGCCAGTTATCCCTGTGGTAACTTTTCTGACACCTCTAGCTTCAAATTCCGAAGGACTAAAGGATCGATAGGCCACGCTTTCACGGTTCGTATTCGTACTGGAAATCAGAATCAAACGAGCTTTTACCCTTTTGTTCCACACGAGATTTCTGTTCTCGTTGAGCTCATCTTAGGACACCTGCGTTATCTTTTAACAGATGTGCCGCCCCAGCCAAACTCCCCACCTGACAATGTCTTCCGCCCGGATCGACCGGCCGAAGCCGACCTTGGGTCCAAAAAGAGGGGCAGTGCCCCGCTTCCGATTCACGGAATAAGTAAAATAACGTTAAAAGTAGTGGTATTTCACTTTCGCTGTTTCCAGCTCCCACTTATCCTACACCTCTCAAGTCATTTCACAAAGTCGGACTAG"
    data = open("TargetSoybean/SoybeanCDS.fna").readlines()
    for key in data:
        if iden == 1:
            if key[0] == ">":
                iden = 2
            else:
                print(key.strip())
        if iden == 2:
            break
        if key[0] == ">" and iden == 0:
            count += 1
            if key.split(">")[1].split(" ")[0] == "glyma.Wm82.gnm1.ann1.Glyma13g11866.1":
                print("glyma.Wm82.gnm1.ann1.Glyma13g11866.1")
                iden = 1
    print(count)


def NewTapirV1():
    # 记录microRNA
    optseq = []
    OptSeq = open("Mapping/OptSample_seq.fasta").readlines()
    for opt in OptSeq:
        temp = list(opt.strip()[::-1])
        for i in range(len(temp)):
            if temp[i] == 'G':
                temp[i] = 'C'
            elif temp[i] == 'C':
                temp[i] = 'G'
            elif temp[i] == 'U':
                temp[i] = 'A'
            else:
                temp[i] = 'T'
        optseq.append("".join(temp))
    print(len(optseq))

    mRnaSeq = {}
    MRNASeq = open("TargetSoybean/SoybeanCDS.fna").readlines()
    seq = ""
    key = ""
    count = 0
    for mRna in MRNASeq:
        if mRna[0] == '>':
            count += 1
            if count == 2:
                print("Key:", key, ",Seq:", seq)
                break
            key = mRna.strip()
        else:
            seq += mRna
    print(len(mRnaSeq))

    ###demo
    str = "ATGATAGGAAGAGCCGACATCGAAGGATCAAAAAGCAACGTCGCTATGAACGCTTGGCTGCCACAAGCCAGTTATCCCTGTGGTAACTTTTCTGACACCTCTAGCTTCAAATTCCGAAGGACTAAAGGATCGATAGGCCACGCTTTCACGGTTCGTATTCGTACTGGAAATCAGAATCAAACGAGCTTTTACCCTTTTGTTCCACACGAGATTTCTGTTCTCGTTGAGCTCATCTTAGGACACCTGCGTTATCTTTTAACAGATGTGCCGCCCCAGCCAAACTCCCCACCTGACAATGTCTTCCGCCCGGATCGACCGGCCGAAGCCGACCTTGGGTCCAAAAAGAGGGGCAGTGCCCCGCTTCCGATTCACGGAATAAGTAAAATAACGTTAAAAGTAGTGGTATTTCACTTTCGCTGTTTCCAGCTCCCACTTATCCTACACCTCTCAAGTCATTTCACAAAGTCGGACTAG"
    s1 = "CCGACCTTGGGTCCAAAAAGAG"
    for i in range(0, len(str) - len(s1) - 1):
        if str[i:i + len(s1)] == s1:
            print("OOO")


def staticsSoybean():
    data = open("TargetSoybean/NewSoybeanOpr1_768.txt").readlines()
    count = 0
    temps = {}
    headers = ["mRNA", "count", "score", "mfe", "misMatch", "GU"]
    temp = []
    iden = ""
    for k in data:
        count += 1
        if count < 7:
            continue
        if (count - 5) % 16 == 2:
            if iden != "":
                temps[iden] = temp
            temp = []
            iden = k.strip().split(".")[-2] + "." + k.strip().split(".")[-1].strip()
            if temps.__contains__(iden):
                temp.append(iden)
                temp.append(temps[iden][1]+1)
            else:
                temp.append(iden)
                temp.append(1)
        if (count - 5) % 16 == 3:
            if temps.__contains__(iden):
                temp.append(temps[iden][2] + float(k.strip().split("          ")[1]))
            else:
                temp.append(float(k.strip().split("          ")[1]))
        if (count - 5) % 16 == 4:
            if temps.__contains__(iden):
                temp.append(round(temps[iden][3] + float(k.strip().split("-")[1]), 1))
            else:
                temp.append(float(k.strip().split("-")[1]))
        if (count - 5) % 16 == 8:
            if temps.__contains__(iden):
                temp.append(temps[iden][4] + int(k.strip().split("       ")[1]))
            else:
                temp.append(int(k.strip().split("       ")[1]))
        if (count - 5) % 16 == 9:
            if temps.__contains__(iden):
                temp.append(temps[iden][5] + int(k.strip().split("             ")[1]))
            else:
                temp.append(int(k.strip().split("             ")[1]))

    contents = []
    result = list(temps.values())
    for res in result:
        coun = res[1]
        for i in range(2, 6):
            res[i] = round(res[i]/coun, 1)

    sortLists = []
    mean = 0
    for i in result:
        mean = i[1] + i[2] + i[3]/10 + i[4] + i[5]
        sortLists.append([i[0], mean])

    sortLists = sorted(sortLists, key=lambda k: -k[1])
    seqCount = 0
    for i in range(len(sortLists)):
        if sortLists[i][1] > 10:
            seqCount += 1
            print(sortLists[i][0])
    print(seqCount)

    count = [i[1] for i in result]
    score = [i[2] for i in result]
    mfe = [i[3] for i in result]
    misMatch = [i[4] for i in result]
    GU = [i[5] for i in result]


    # 写文件
    # with open("TargetSoybean/NewSoybeanStaticRumdupMean.csv", "w") as f:
    #     f_csv = csv.writer(f)
    #     f_csv.writerow(headers)
    #     f_csv.writerows(list(temps.values()))

    # 作图
    # import matplotlib.pyplot as plt
    #
    # plt.subplot(221)
    # plt.title('Score distribution')
    # plt.ylabel("score_value")
    # plt.xlabel("mRNA_sequence")
    # plt.scatter([i + 1 for i in range(len(score))], score, marker='.')
    # plt.subplot(222)
    # plt.title('Mfe distribution')
    # plt.ylabel("mfe_value")
    # plt.xlabel("mRNA_sequence")
    # plt.scatter([i + 1 for i in range(len(mfe))], mfe, marker='+')
    # plt.subplot(223)
    # plt.title('misMatch distribution')
    # plt.ylabel("misMatch_value")
    # plt.xlabel("mRNA_sequence")
    # plt.scatter([i + 1 for i in range(len(misMatch))], misMatch, marker='.')
    # plt.subplot(224)
    # plt.title('GU distribution')
    # plt.ylabel("GU_value")
    # plt.xlabel("mRNA_sequence")
    # plt.scatter([i + 1 for i in range(len(GU))], GU, marker='+')
    #
    # count.sort(reverse=True)
    # f, (ax, ax2) = plt.subplots(2, 1, sharex=True)
    # ax.set_ylim(12, 4000)  # outliers only
    # ax2.set_ylim(0, 11)  # most of the data
    # ax.spines['bottom'].set_visible(False)
    # ax2.spines['top'].set_visible(False)
    # ax.xaxis.tick_top()
    # ax.tick_params(labeltop='off')  # don't put tick labels at the top
    # ax2.xaxis.tick_bottom()
    # d = .015  # how big to make the diagonal lines in axes coordinates
    # # arguments to pass to plot, just so we don't keep repeating them
    # kwargs = dict(transform=ax.transAxes, color='k', clip_on=False)
    # ax.plot((-d, +d), (-d, +d), **kwargs)  # top-left diagonal
    # ax.plot((1 - d, 1 + d), (-d, +d), **kwargs)  # top-right diagonal
    # kwargs.update(transform=ax2.transAxes)  # switch to the bottom axes
    # ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    # ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal
    # ax.scatter([i + 1 for i in range(10)], count[:10], marker='.', c='r')
    # ax2.scatter([i + 1 for i in range(len(count))], count, marker='.', c='r')
    # plt.show()


if __name__ == '__main__':
    # static_In()
    get_sojae_protein_id()
    # get_soybean_protein_id()
    ### staticsSoybean()

    # staticSojae()

    # data = open("TargetSoybean/SoybeanCDS.fna").readlines()
    # res = ["Glyma05g01470"]
    # iden = 0
    # count = 0
    # for key in data:
    #     if iden == 1:
    #         if key[0] == ">":
    #             iden = 2
    #         else:
    #             print(key.strip())
    #     if iden == 2:
    #         break
    #     if iden == 0 and key[0]==">":
    #         count += 1
    #         if key.strip().split(" ")[0].split(".")[-2] in res:
    #             print(key.strip().split(" ")[0].split(".")[-2].strip())
    #             iden = 1
    # NewTapirV1()
    # cut_value()
    # get_U()
    # AliasSojae()a
    # T1to2()
    # demo()
    # Aliases()
    # cal_18T25()


