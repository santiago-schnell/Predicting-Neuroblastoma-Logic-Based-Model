from boolean2 import Model, util, state, network
import matplotlib.pyplot as plt

Bool1 = file('Test_Rules.txt').read()

### An Asynchronous model

def BoolMod(timecourse, booldata, str_mode):
    coll = util.Collector()
    model = Model(text= booldata, mode= str_mode)
    initializer = state.all_initial_states(model.nodes, limit=None)

    # the data is the inital data, the func is the initializer
    n=0
    for data, initfunc in initializer:
        # shows the initial values
        print ((n)/float(32768))*100, "% done."
        n = n + 1
        model.initialize(missing=initfunc)
        model.iterate(steps=timecourse)
        # takes all nodes
        nodes = model.nodes
        coll.collect(states=model.states, nodes=nodes)
    # --------------- Detect Cycles ------------------#
    print model.report_cycles()
    # Return results
    avgs = coll.get_averages(normalize=True)
    data = pd.DataFrame(avgs)
    pd.DataFrame.to_csv(data, "test.csv")
    print avgs
    print model.fp()
    # ------------- Retrieve the Data -------------#
    # 6-gene input signature
    #ALK = avgs.get('ALK')
    #MDK = avgs.get('MDK')
    #TrkA = avgs.get('TrkA')
    #NGF = avgs.get('NGF')
    #TrkB = avgs.get('TrkB')
    #BDNF = avgs.get('BDNF')

    # Model outcome states
    Differentiation = avgs.get('Differentiation')
    Apoptosis = avgs.get('Apoptosis')
    Proliferation = avgs.get('Proliferation')
    Angiogenesis = avgs.get('Angiogenesis')

    # Nodes at Issue
    #DNADamage = avgs.get('DNADamage')
    #p53 = avgs.get('p53')

    # Other Nodes
    #MDM2 = avgs.get('MDM2')
    #MAPK = avgs.get('MAPK')
    #p27 = avgs.get('p27')
    #FoxO = avgs.get('FoxO')
    #AKT = avgs.get('AKT')
    #Ras = avgs.get('Ras')
    #MYCN = avgs.get('MYCN')
    #MTOR = avgs.get('MTOR')
    #IP3 = avgs.get('IP3')

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
    # x.plot(t, AKT, label = 'AKT')
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

BoolMod(10,Bool1,'async')