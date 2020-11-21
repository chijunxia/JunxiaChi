import logging
import csv
from ML.feature import feature
logging.getLogger().setLevel(logging.INFO)


def preProcess():
    # 计算最小自由能
    mdata = open('../Script/Mapping/Mfe/OptMfe.txt', 'r').readlines()
    count = 0
    mfe = []
    logging.info("START CALCULATE MFE")
    for key in mdata:
        count = count + 1
        if count % 3 == 2:
            s = key.strip()
        if count % 3 == 0:
            s = s + "," + key.split(" ")[-1].split("-")[-1].split(")")[0]
            mfe.append(s)
    logging.info("END CALCULATE MFE")

    # 计算Motif
    logging.info("START Pre MOTIF")
    level1 = ["A", "G", "C", "U"]
    level2 = []
    level3 = []
    for i in level1:
        for j in level1:
            level2.append(i + j)
    for i in level1:
        for j in level2:
            level3.append(i + j)
    logging.info("END Pre MOTIF")

    # 转码
    encode = {
        "N": [0, 0, 0, 0],
        "A": [0, 0, 0, 1],
        "G": [0, 0, 1, 0],
        "C": [0, 1, 0, 0],
        "U": [1, 0, 0, 0]
    }
    headers = feature
    data = open('../Script/Mapping/OptSample_seq.fasta', 'r').readlines()
    ls = []

    logging.info("START CALCULATE ALLContent")
    for key in data:
        tcount = 0
        temp = []
        for i in key.strip():
            if i == 'G' or i == 'C':
                tcount = tcount + 1
        mt = ""
        for m in mfe:
            if key.strip() == m.split(",")[0]:
                mt = round(float(m.split(",")[1]))
        temp = [key.strip(), len(key.strip()), round(tcount / len(key.strip()) * 10), mt]

        addrls = []
        for ik in range(len(key.strip())):
            addrls.extend(encode[key[ik]])
        for ik in range(len(key.strip()), 25):
            addrls.extend(encode["N"])
        addrls.extend(encode[key[0]])
        addrls.extend(encode[key[1]])
        addrls.extend(encode[key[len(key.strip()) - 2]])
        addrls.extend(encode[key[len(key.strip()) - 1]])
        temp.extend(addrls)
        # level1
        for l1 in level1:
            acount = 0
            t1 = ""
            for i1 in key.strip():
                t1 = t1 + i1
                if i1 == l1:
                    acount = acount + 1
            temp.append(acount)
        # level2
        for l2 in level2:
            acount = 0
            for i2 in range(0, len(key.strip()) - 3):
                if l2 == key[i2:i2 + 2]:
                    acount = acount + 1
            temp.append(acount)
        # level3
        for l3 in level3:
            acount = 0
            for i3 in range(0, len(key.strip()) - 4):
                if l3 == key[i3:i3 + 3]:
                    acount = acount + 1
            temp.append(acount)
        temp.append(0)
        ls.append(temp)
    logging.info("END CALCULATE ALLContent")

    logging.info("START Writting")
    with open("../Script/Mapping/Mfe/AllSample_0314.csv", "w", newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
        f_csv.writerows(ls)
    logging.info("END Writting")


if __name__=="__main__":
    preProcess()

