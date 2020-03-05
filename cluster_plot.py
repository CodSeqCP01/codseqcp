#!/usr/bin/env python

import sys, codseqcp
from collections import Counter
import matplotlib.pyplot as plt

seq_dict = codseqcp.parsefasta(open(sys.argv[1]))

Clusters = []
for key, seq in seq_dict.items():
    sub_lst = []
    flag = "_".join(key.split("_")[:2])
    Clusters.append(flag)

occur = Counter([val for clust, val in Counter([clust for clust in Clusters]).items()])
plt.bar(occur.keys(), occur.values(), color='b')
plt.ylabel('Number of homologous loci', fontsize='large', fontweight='bold')
plt.xlabel('Number of samples per locus', fontsize='large', fontweight='bold')
plt.show()






    
