import BooleanFunct as bf

df = bf.datahandler('Allison.txt', geneon=[], geneoff=[], csv_out_txt='Allison_Test_1.csv',debug=True, nodes_for_averages=['AMPK','mTORC1','Leu','Sestrin'], graph=True)

bf.findlc('Allison_Test_1.csv', savename='Allison_Test_1_lc.csv')

