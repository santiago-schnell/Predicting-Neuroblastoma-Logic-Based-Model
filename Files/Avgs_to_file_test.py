import pandas as pd

print
x = pd.read_csv('GalleyRulesRemoved_TRKBON_DL.csv')
df = x.loc[x['CycleLength'] == 1]
print df