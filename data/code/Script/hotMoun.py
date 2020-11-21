# -*- coding: utf-8 -*-
from __future__ import division

import numpy as np
from scipy import stats


def hotMoun():
    import pandas as pd  # Data analysis
    import numpy as np  # Scientific computing
    import matplotlib.pyplot as plt  # Plotting
    import matplotlib.colors as colors  # Coloring
    import seaborn as sns  # Statistical visualization

    fold = [-1, -2, 0, -3, 6, 5, 8, 1]
    pvalue = [0.03, 0.04, 0.02, 0.005, 0.001, 0.04, 0.05, 0.001]

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


def trandn_p():
    means = [0.0539, 4, 8, 3, 6, 9, 1]
    stds = [5, 4, 8, 3, 6, 7, 9]
    ms = [0, 4.1, 7, 2, 5, 8, 0]
    n = 7
    output = []
    for mean, std, m in zip(means, stds, ms):
        tt = (mean - m) / (std / np.sqrt(float(n)))  # t-statistic for mean
        pval = stats.t.sf(np.abs(tt), n - 1) * 2  # two-sided pvalue = Prob(abs(t)>tt)
        print('t-statistic = %6.3f pvalue = %6.4f' % (tt, pval))
        output.append(format(pval))
    print("\t".join(output))


def frac(n):
    r = 1
    if n==0 or n==1:
        return 1
    for i in range(1, n+1):
        r *= i
    return r


def get_pvalue():
    data = open('Mapping/New_c7_c6.csv').readlines()
    N = 6399922
    M = 4108080
    pvalues = []
    import scipy.special
    count = 0

    for d in data:
        k = int(d.strip().split(",")[1])
        n = k + int(d.strip().split(",")[2])
        # p = scipy.special.comb(M, k)*scipy.special.comb(N-M, n-k)/scipy.special.comb(N, n)
        print(scipy.special.comb(M, k),":",scipy.special.comb(N-M, n-k),":",scipy.special.comb(N, n))

        count += 1
        if count>1:
            break

    print(int(scipy.special.perm(3, 2)))
    print(int(scipy.special.comb(3, 2)))


if __name__ == '__main__':
    print()
