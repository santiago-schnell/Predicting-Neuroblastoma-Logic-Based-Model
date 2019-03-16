# Export data to csv
import csv
import pandas as pd

import json
headers = ['numbers', 'letters']
d = [[1,2,3,4,5,6,7,8,9,10,11],['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']]



# for i in range(3):
#     with  open('testappend.csv','a') as out:
#         csv_out = csv.writer(out)
#         csv_out.writerow(headers)
#         for row in d:
#             csv_out.writerow(row)