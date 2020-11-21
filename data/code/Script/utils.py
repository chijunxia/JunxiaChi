#coding=utf-8
from __future__ import division
from Decor import Log
import csv
import logging
from config import *
logging.getLogger().setLevel(logging.INFO)


class ReadFile:
    def __init__(self):
        self.count = 0
        self.Seqs = []
        self.Dicts = {}

    # 本函数用于将文件内容存储到list中
    @Log.loginfo
    def countDownL(self, path):
        l = []
        logging.info("Start L")
        for line in open(path):
            l.append(line.strip())
            #logging.info(line.strip().split(',')[0])
        logging.info("List: %s", len(l))
        return l

    @Log.loginfo
    def countDownD(self, l):
        d = {}
        for k in l:
            d[k] = d.get(k, 0) + 1
        logging.info("D ok!")
        return d

    @Log.loginfo
    def countDownLs(self, l, path):
        logging.info("Start Ls")
        d = {}
        ls = []

        for k in l:
            d[k] = d.get(k, 0) + 1
        logging.info("Ls d ok!")

        for key in d:
            ls.append((key, len(key), d[key]))
        logging.info("Ls ok!")

        newls = sorted(ls, key=lambda k: -k[2])

        headers = ['Reads','Length','Count']
        logging.info("Start Write to csv")
        with open(path, "w") as f:
            f_csv = csv.writer(f)
            f_csv.writerow(headers)
            f_csv.writerows(ls)
        logging.info("Ls Write SUCCESS")
        return newls

    def hotMoun(self, fold, pvalue):
        import pandas as pd  # Data analysis
        import numpy as np  # Scientific computing
        import matplotlib.pyplot as plt  # Plotting
        import seaborn as sns  # Statistical visualization

        fold1 = pd.DataFrame(fold)
        fold1['a'] = range(len(fold1))
        fold1 = fold1.set_index('a')

        pvalue1 = pd.DataFrame(pvalue)

        result = pd.concat([fold1, pvalue1], axis=1, ignore_index=True)

        result.columns = ['fold', 'pvalue']
        result['log(pvalue)'] = -np.log10(result['pvalue'])
        result['sig'] = 'normal'
        result['size'] = np.abs(result['fold']) / 10
        result.loc[(result.fold > 1) & (result.pvalue < 0.05), 'sig'] = 'up'
        result.loc[(result.fold < -1) & (result.pvalue < 0.05), 'sig'] = 'down'

        ax = sns.scatterplot(x="fold", y="log(pvalue)",
                             hue='sig',
                             hue_order=('down', 'normal', 'up'),
                             palette=("#377EB8", "grey", "#E41A1C"),
                             data=result)
        ax.set_ylabel('-log(pvalue)', fontweight='bold')
        ax.set_xlabel('FoldChange', fontweight='bold')
        plt.show()

    def countDownAll(self, sinD, mixD, Ain, Ano, inD):
        import math
        logging.info("Start All")
        lsIn = []
        lsNo = []
        lsflag = 0
        lscount = 0
        all_count = 0

        fold = []
        pvalue = []

        for kmix in mixD:
            lscount = lscount+1
            if kmix in sinD.keys():
                if kmix in inD.keys():
                    # fold.append(math.log(mixD.get(kmix)/sinD.get(kmix), 2))
                    # pvalue.append(mixD.get(kmix)/4108080)
                    # lsIn.append((kmix, len(kmix), mixD.get(kmix), mixD.get(kmix)/sinD.get(kmix)-1))
                    # all_count += int(mixD.get(kmix))
                    lsIn.append((kmix, mixD.get(kmix), sinD.get(kmix)))
                lsflag = 1
            # if lsflag == 0:
            #     if kmix in inD.keys():
            #         lsNo.append((kmix, len(kmix), mixD.get(kmix)))
            # lsflag = 0
        logging.info(all_count)
        # self.hotMoun(fold, pvalue)

        newlsIn = sorted(lsIn, key=lambda k: -k[2])
        newlsNo = sorted(lsNo, key=lambda k: -k[2])

        headers = ['Reads', 'Count7', 'Count6']
        logging.info("Start Write to first csv")
        print(len(lsIn))
        with open(Ain, "w", newline='') as f:
            f_csv = csv.writer(f)
            f_csv.writerow(headers)
            f_csv.writerows(newlsIn)
        logging.info("Ls Write first SUCCESS")
        logging.info(len(newlsIn))
        return
        # headers = ['Reads', 'Length', 'Count']
        # logging.info("Start Write to two csv")
        # print(len(newlsNo))
        # with open(Ano, "w") as f:
        #     f_csv = csv.writer(f)
        #     f_csv.writerow(headers)
        #     f_csv.writerows(newlsNo)
        # logging.info("Ls Write two SUCCESS")

    @Log.loginfo
    def cutNoDump(self):
        lDup = []
        for line in open(r"data/c707_rmdup_result_only2.fastq"):
            lDup.append(line.strip())
        logging.info("List: %s", len(lDup))

        lNo = []
        for line in open(r"data/RateNoSort.csv"):
            lNo.append(line.strip().split(",")[0])

        lInter = list(set(lDup).intersection(set(lNo)))
        lNo_Inter = list(set(lNo).difference(set(lInter)))

        logging.info(len(lNo_Inter))

        NoDict = {}
        NoDictCutDup = {}
        for line in open(r"data/RateNoSort.csv"):
            NoDict[line.strip().split(",")[0]] = line.strip().split(",")[2]

        logging.info(len(NoDict))

        count = 0
        for key in NoDict:
            if key not in lNo_Inter:
                count = count+1
                NoDictCutDup[key] = NoDict[key]
                logging.info(count)
        logging.info(len(NoDictCutDup))

        newNoCutDup = []
        for key in NoDictCutDup:
            newNoCutDup.append((key, len(key), NoDictCutDup[key]))

        headers = ['Reads', 'Length', 'Count']
        logging.info("Start Write to two csv")
        with open(r"data/RateNoCutDup.csv", "w") as f:
            f_csv = csv.writer(f)
            f_csv.writerow(headers)
            f_csv.writerows(newNoCutDup)

    @Log.loginfo
    def ToFasta(self):
        TAG = ">ssz-in"
        l = []
        count = 0
        for line in open(r"data/RateInSort.csv"):
            count = count+1
            if count>=1000 and count <2000:
                l.append(TAG + str(count))
                l.append(line.strip().split(",")[0].replace('T','U'))
            elif count<1000:
                continue
            else:
                break
        fp = open("data/sszInUseq1000_1999.fasta", 'w+')
        for i in l:
            fp.write(i)
            fp.write("\n")
        fp.close()


    def get_geneId(self):
        tapirs = TapirFiles
        TaFile = sequence
        urls = ReadFile.get_term(tapirs, TaFile)
        logging.info("url-list: %s", len(urls))
        urls = list(set(urls))
        logging.info("url-set: %s", len(urls))


        from bs4 import BeautifulSoup
        import requests
        geneIds = []
        logging.info("Start Get geneIds")
        headers = {'User-Agent': 'User-Agent:Mozilla/5.0'}
        lenL = len(urls)
        count = 1;
        countN = 0;
        fp = open("data/geneIdsLists.txt", 'w+')
        for u in urls:
            with requests.get(url=r"https://www.ncbi.nlm.nih.gov/search/all/?term="+u, headers=headers) as response:
                data2 = response.content.decode()
                soup = BeautifulSoup(data2, 'html.parser')
                try:
                    count = count + 1
                    links = soup.find_all('div', class_="ncbi-doc-id")
                    link = links[0].find_all('ul', class_="ncbi-inline-list")[0].get_text()
                    geneIds.append(link.strip()[8:])
                    fp.write(link.strip()[8:])
                    fp.write("\n")
                    # print(link.strip()[8:])
                    print(count+countN+"("+ count +") in " +lenL)
                except:
                    countN = countN+1
                    continue
        fp.close()

    @staticmethod
    def get_term(tapirs, sequence):
        l = []
        for i in tapirs:
            count = 0
            for line in open(i):
                count = count+1
                if count%16 == 2:
                    l.append(line.strip()[15:])
                    # logging.info(line.strip()[15:])
        logging.info("len-list: %s", len(l))
        l = list(set(l))
        logging.info("len-set: %s", len(l))
        urlParms = []
        for line in open(sequence):
            if line[0:5] == ">lcl|":
                if len(line.split("[")[1].strip()[10:27]) == 17:
                    urlParms.append(line.split("[")[1].strip()[10:27])
                # print(line.split("[")[1].strip()[10:27])
        logging.info(">lcl| match success!")
        return urlParms

    @staticmethod
    def get_protainId():
        tapirs = TapirFiles
        seqs = sequence
        l = []
        for i in tapirs:
            count = 0
            for line in open(i):
                count = count + 1
                if count % 16 == 2:
                    l.append(line.strip()[15:])
        l = list(set(l))
        logging.info("len-set: %s", len(l))
        protainIds = []

        outfile = open('Script/ssz1Protauns.txt', 'w')
        for line in open(sequence):
            if line[0:5] == ">lcl|":
                if line.split(" ")[0][1:] in l:
                    #print(line[line.find("protein_id")+11:].split("]")[0])
                    #protainIds.append(line[line.find("protein_id")+11:].split("]")[0].split(".")[0])
                    outfile.write(line[line.find("protein_id")+11:].split("]")[0].split(".")[0] + '\n')
        outfile.close()


if __name__ == '__main__':
    print("A")

