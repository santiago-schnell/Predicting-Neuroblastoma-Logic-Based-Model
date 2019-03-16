from boolean2 import Model, util, state
import matplotlib.pyplot as plt
import networkx
import csv
import numpy as np
import pandas as pd



Bool1 = """
#----------Initial States----------#
#6-gene input signature
ALK = True
MDK = False
TrkA = Random
NGF = True
TrkB = False
BDNF = False

#Model outcome states
Differentiation = Random
Apoptosis = Random
Proliferation = Random
Angiogenesis = Random

#Nodes at Issue -did some fixing base upon quick lit review
DNADamage = Random
p53 = Random

#Other Nodes
MDM2 = Random
MAPK = Random
p27 = Random
FoxO = Random
AKT = Random
Ras = Random
MYCN = Random
MTOR = Random
IP3 = Random

#---------updating rules----------#

#6-gene input signature
ALK* = ALK
MDK* = MDK
TrkA* = not MYCN
NGF* = NGF
TrkB* = TrkB
BDNF* = BDNF

#Model outcome states
Differentiation* = MAPK
Apoptosis* = (p53 and not AKT) or (TrkA and not NGF)
Proliferation* = (not p27 and not p53) or IP3
Angiogenesis* = MTOR

#Nodes at Issue -did some fixing base upon quick lit review
DNADamage* = DNADamage
p53* = DNADamage and not MDM2

#Other Nodes
MDM2* = p53
MAPK* = (MDK and ALK) or Ras
p27* = FoxO or not MYCN
FoxO* = not AKT
AKT* = (MDK and ALK) or (BDNF and TrkB)
Ras* = NGF and TrkA
MYCN* = (AKT or Ras) and not TrkA
MTOR* = AKT
IP3* = BDNF and TrkB

"""


### An Asynchronous model


def BoolModPlot(timecourse, samples, booldata, str_mode):
    coll = util.Collector()
    for i in range(samples):
        model = Model(text= booldata, mode= str_mode)
        model.initialize()
        model.iterate(steps=timecourse)
        print model.states

        # takes all nodes
        nodes = model.nodes
        coll.collect(states=model.states, nodes=nodes)

    # --------------- Detect Cycles ------------------#
    print model.report_cycles()
    # Return results
    avgs = coll.get_averages(normalize=True)
    data = pd.DataFrame(avgs)
    pd.DataFrame.to_csv(data,"test.csv")
    # ------------- Retrieve the Data -------------#
    # 6-gene input signature
    ALK = avgs.get('ALK')
    MDK = avgs.get('MDK')
    TrkA = avgs.get('TrkA')
    NGF = avgs.get('NGF')
    TrkB = avgs.get('TrkB')
    BDNF = avgs.get('BDNF')

    # Model outcome states
    Differentiation = avgs.get('Differentiation')
    Apoptosis = avgs.get('Apoptosis')
    Proliferation = avgs.get('Proliferation')
    Angiogenesis = avgs.get('Angiogenesis')

    # Nodes at Issue
    DNADamage = avgs.get('DNADamage')
    p53 = avgs.get('p53')

    # Other Nodes
    MDM2 = avgs.get('MDM2')
    MAPK = avgs.get('MAPK')
    p27 = avgs.get('p27')
    FoxO = avgs.get('FoxO')
    AKT = avgs.get('AKT')
    Ras = avgs.get('Ras')
    MYCN = avgs.get('MYCN')
    MTOR = avgs.get('MTOR')
    IP3 = avgs.get('IP3')

      ### Time axis (x)
    t = range(0, timecourse + 1)

    # Create plots with pre-defined labels. Try to make this a for loop
    fig, ax = plt.subplots()
    ax.plot(t, Differentiation, label='Differentiation')
    ax.plot(t, Apoptosis, label='Apoptosis')
    ax.plot(t, Angiogenesis, label='Angiogenesis')
    ax.plot(t, Proliferation, label='Proliferation')
    # ax.plot(t, TrkA, label='TrkA')
    # ax.plot(t, TrkB, label = 'TrkB')
    # ax.plot(t, MYCN, label = 'MYCN')
    # ax.plot(t, NGF, label = 'NGF')
    # ax.plot(t, MDK, label = 'MDK')
    # ax.plot(t, ALK, label = 'ALK')
    # ax.plot(t, Ras, label = 'Ras')
    # ax.plot(t, AKT, label = 'AKT')
    # ax.plot(t, FoxO, label = 'FoxO')
    # ax.plot(t, p27, label = 'P27')
    # ax.plot(t, p53, label = 'P53')

    # Legends
    legend = ax.legend(loc=0, shadow=True, fontsize='medium')
    plt.xlabel('Iterations')
    plt.ylabel('On Proportion')
    plt.title("Asynchronous Updating Model")

    # beautify
    legend.get_frame().set_facecolor('#C0C0C0')

    # View Plot
    plt.show()

    return avgs

BoolModPlot(100,30,Bool1,'async')





