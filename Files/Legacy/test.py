# import networkx as nx
# import matplotlib.pyplot as plt
# import Network_Repair_Functions as nf
#
# f = open(r'MapRules.txt','r')
# lines = f.readlines()
# f.close()
#
# # We construct a network based on these rules
# G,nodes = nf.form_network(lines)
#
# nx.draw(G)
# plt.savefig("path.png")

### Write dictionary to csv with keys as headers

import csv
import numpy as np

data1 = np.arange(10)
data2 = np.arange(10)*2
data3 = np.arange(10)*3

writefile = 'test.csv'

fieldnames = ['ALK','MDK', 'TrkA', 'NGF', 'TrkB', 'BDNF', 'Differentiation', 'Apoptosis', 'Proliferation', \
              'Angiogenesis', 'DNADamage', 'p53', 'MDM2', 'MAPK', 'p27', 'FoxO', 'AKT', 'Ras', 'MYCN', 'MTOR', 'IP3']
with open( writefile, 'w' ) as f:
    writer = csv.writer(f)
    writer.writerow(fieldnames)
    writer.writerows(ALK,MDK, TrkA, NGF, TrkB, BDNF, Differentiation, Apoptosis, Proliferation, Angiogenesis, \
                     DNADamage, p53, MDM2, MAPK, p27, FoxO, AKT, Ras, MYCN, MTOR, IP3)