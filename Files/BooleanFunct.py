from boolean2 import Model, state, tokenizer
import csv
import pandas as pd
import matplotlib.pyplot as plt


# TODO: Figure out why the heck the timecourse var is not stopping the model iteration when using the async setting
def datahandler(str_file, str_mode='sync', timecourse=10, csv_out_txt="default.csv", geneon=[], geneoff=[], dumpevery=1000,\
                nodes_for_averages=['Apoptosis', 'Proliferation', 'Angiogenesis', 'Differentiation'], debug=False, graph=False):
    '''
    A function to generate a list of all possible starting states of a model. Note that this function is a
    modification of the all_ss_model() function and allows for specific nodes to be fixed on or off.

    str_file: as a string, write the name of the file containing the rules for the boolean model
    str_mode: enter the mode of model that is going to be used, ie sync or async. Default is sync
    timecourse: The number of time steps that the model will progress through. Default is ten
    csv_out_text: the results are output to this CSV file. Default is 'default.csv.'
    geneon: Takes a list of nodes as strings. These nodes will be forced ON for all simulations
    geneoff: Takes a list of nodes as strings. These nodes will be forced OFF for all simulations
    dumpevery: utilizes the dumper subfunction to write the data to file frequently to avoid running out of RAM.
        Default is to dump every 1000 simulations
    nodes_for_averages: Input a list of nodes as strings. The program will keep track of the average ON state of these
        nodes. The default values are for use with the BPC model of neuroblastoma (2018)
    debug: Default as false. When True will print the progress of the simulation. Will slightly slow rate of simulation
        when on
    graph: Boolean literal, if true will display a
    :return:
    A list of dictionaries with all possible starting states of the network.
    '''

    def dumper():
        '''
        sub-function
        Will dump the data into a CSV file of name csv_out_txt, keeps memory burden low for extremely large networks.
        :return: none
        '''
        with open(csv_out_txt, 'a') as csvfile:
            datadumper = csv.writer(csvfile, lineterminator='\n')
            for row in d:
                datadumper.writerow(row)
        return None

    # Header names for pandas dataframe and CSV data dump file
    headers = ['Index', 'CycleLength', 'CycleFingerprint', 'FirstState', 'LastState', 'FirstFP', 'SteadyState', str(nodes_for_averages[0]),\
               str(nodes_for_averages[1]), str(nodes_for_averages[2]), str(nodes_for_averages[3])]

    # Init a CSV to contain the data generated over the course of the run
    with  open(csv_out_txt,'wb') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(headers)

    # Variable saving the imported list of rules for use by the model.
    Bool1 = file(str_file).read()

    # Generate a new version of the model that allows the genes to be fixed on or off by removing updating rules for
    # them. This will not override the initialized value however, and that must be updated to be assigned on or off
    # separately.
    on = geneon
    off = geneoff
    Bool2 = tokenizer.modify_states(Bool1, on, off)
    model = Model(text=Bool2, mode=str_mode)
    initializer = state.all_initial_states(model.nodes, limit=None)

    # Utilized in the for loops for the loading bar
    load_status_divisor = len(model.nodes) - len(geneoff) - len(geneon)

    # The BooleanNet Data Collector. Here it is implemented to gather data on the states of the nodes in the model.
    coll = util.Collector()

    # Wiley Stoeber. 5/8/18. Create a modified version of the initializer that will pass over initial states that
    # contradict the gene set mode
    initializer_new = []

    # An initializer for single KO mutants.
    # TODO I need to make it so this initializer can adapt to any number of KO/OE genes. As it stands it can only do up
    # to a double. Not helpful! Need to make recursive?
    if geneoff != [] or geneon != []:
        if len(geneoff) + len(geneon) == 1:
            d = []
            p = 0
            if debug:
                n = 0
            for data_init in initializer:
                data = data_init[0]
                for i in range(len(geneoff)):
                    if not data[geneoff[i]]:
                        initializer_new.append(data_init)
                        if debug:
                            print 'Single Gene Knockout, Geneoff'
                for i in range(len(geneon)):
                    if data[geneon[i]]:
                        initializer_new.append(data_init)
                        if debug:
                            print 'Single Gene Overexpression, Geneon'

        # An initializer for double KO when some on and some off
        elif len(geneoff) >= 1 and len(geneon) >= 1:
            d = []
            p = 0
            if debug:
                n = 0
            for data_init in initializer:
                data = data_init[0]
                if not data[geneoff[0]]:
                    if data[geneon[0]]:
                        initializer_new.append(data_init)
                        if debug:
                            print 'Gene off then gene on, KO then Overexpression'
                elif data[geneon[0]]:
                    if not data[geneoff[0]]:
                        initializer_new.append(data_init)
                        if debug:
                            print 'Gene on then gene off, Overexpression then KO'
                else:
                    pass


        # An initializer for double KO when both genes are on or off
        elif len(geneoff) == 2 or len(geneon) == 2:
            d = []
            p = 0
            if debug:
                n = 0
            for data_init in initializer:
                data = data_init[0]
                i = 0
                while i <= 1:
                    if not data[geneoff[i]]:
                        i += 1
                        if not data[geneoff[i]]:
                            initializer_new.append(data_init)
                            if debug:
                                print 'Double Gene Knockout'
                            i += 1
                            break
                        break
                    break

            for data_init in initializer:
                data = data_init[0]
                i = 0
                while i <= 1:
                    if data[geneon[i]]:
                        i += 1
                        if data[geneon[i]]:
                            initializer_new.append(data_init)
                            if debug:
                                print 'Double Gene Overexpression'
                            i += 1
                            break
                        break
                    break
        for data, initfunc in initializer_new:
            # Fixes genes on or off (True or False) at their starting state.
            for i in range(len(geneoff)):
                data.update({str(geneoff[i]): False})
            for i in range(len(geneon)):
                data.update({str(geneon[i]): True})

            # Initialize the model with the given pre-computed initial conditions stored in the data vfariable.
            # for a given model with Z nodes, there are 2 to the power of Z starting states.
            # todo: switch to dynamic length if i is going to go over 10, becomes much less efficient quickly
            model.initialize(defaults=data)
            for i in range(timecourse):
                model.iterate(steps=i)
                e = model.detect_cycles()
                nodes = nodes_for_averages

                if e[1] == 1:
                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    detect_states = 0 - e[1]

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]
                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break

                elif e[1] > 1:

                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    detect_states = 0 - e[1]

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug == True:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]

                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break

                elif e[1] == 0 and i+1 == timecourse:
                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    detect_states = 0 - e[1]

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]
                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break

            # Iterate the model
            n += 1
            p += 1

            # Dump data to file if the iterator p has reached the dumpevery variable. Keeps RAM burden low. Resets the
            # p counter to zero and overwrites the data that had been written to file.
            if p == dumpevery:
                dumper()
                p = 0
                d = []

        # A catchall at the end of the set of simulations that dumps any remaining data that might not have been
        # enough to trigger the prior dumper.
        dumper()
        d = []

    else:
        d = []
        p = 0
        if debug:
            n = 0
        for data, initfunc in initializer:
            # Fixes genes on or off (True or False) at their starting state.
            # Initialize the model with the given pre-computed initial conditions stored in the data variable.
            # for a given model with Z nodes, there are 2 to the power of Z starting states.
            model.initialize(defaults=data)
            for i in range(timecourse):
                model.iterate(steps=i)
                e = model.detect_cycles()
                nodes = nodes_for_averages
                print e

                if e[1] == 1:
                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    detect_states = 0 - e[1]

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    # Affix the limit cycle or FPA that defines the simulation to a new column
                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]
                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break

                elif e[1] > 1:

                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    detect_states = 0 - e[1]

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]
                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break


                # This function makes sure that the entire user specified timecourse has passed before writing
                # that the model did not reach a steady state.
                elif e[1] == 0 and i+1 == timecourse:
                    # append the number of index at which cycles began and the length of the cycle to data
                    d.append(list(model.detect_cycles()))
                    # append the model fingerprint for the entirety of the simulation
                    d[p].append(model.fp())
                    # append the first state of the model, coincides with assigned state
                    d[p].append(model.first)
                    # append the last state of the model, coincides with one part of LC or FPA
                    d[p].append(model.last)
                    # detect states is used to figure out how many states to write to the collector object.
                    detect_states = 0 - e[1]
                    # append the an int that represents the models starting state, useful for later sort functions
                    # especially in excel.
                    d[p].append(model.fp()[0])

                    # Affix the limit cycle or FPA that defines the simulation to a new column.
                    a = model.fp()
                    b = a[detect_states:]
                    c = sorted(b)
                    d[p].append(c)

                    # Collect the averages over the course of the run.
                    coll.collect(states=model.states[detect_states:], nodes=nodes)
                    # Console output for debugging and progress tracking
                    if debug:
                        print 'The fingerprint is', d[p][2]
                        print 'The cycle length is', d[p][1]
                        prc = ((n) / float(pow(2, load_status_divisor))) * 100
                        print '%.2f' % prc + "% done."
                        print '\n'
                    break

            # Iterate the model
            if debug:
                n += 1
            p += 1

            # Dump data to file if the iterator p has reached the dumpevery variable. Keeps RAM burden low. Resets the
            # p counter to zero and overwrites the data that had been written to file.
            if p == dumpevery:
                dumper()
                p = 0
                d = []
        # A catchall at the end of the set of simulations that dumps any remaining data that might not have been enough
        # to trigger the prior dumper.
        dumper()
        d = []

    # Generate a pandas DataFrame of the dump file.

    df = pd.read_csv(csv_out_txt)

    # Create count list of df
    fpa = df[df.CycleLength == 1]
    fpa["SteadyState"] = fpa['SteadyState'].str.strip('[]').astype(int)
    counts = fpa['SteadyState'].value_counts().sort_index()
    counts = pd.DataFrame([counts], columns=["SS","NumberInFPA"])
    print counts

    counts['Proportion'] = counts['NumberInFPA'] / float(pow(2, load_status_divisor))



    # If there are nodes of interest in the nodes_for_averages will output the average on state of the nodes to the
    # dataframe  then update the CSV to contain them
    if nodes_for_averages is not []:
        for i in range(len(nodes_for_averages)):
            avgs = coll.get_averages(normalize=True)
            a1 = avgs[nodes_for_averages[i]]
            df.set_value(0, nodes_for_averages[i], a1[0])
        df.to_csv(csv_out_txt)

    #TODO: I need to create a method of exporting an excel file that is processed so that the user does not need to process them manually to get system state information

    # Graphical output
    # Will export a multicolor graph that shows the average on proportion of
    if graph is True:
        nfa = nodes_for_averages
        x = []
        for i in range(len(nfa)):
            x.append(nfa[i])
        gr = df[x]
        gr = gr.dropna()
        plot = gr.plot(kind='bar')
        plot.set_xlabel("Groups")
        plot.set_ylabel("On Proportion")
        plt.show()

    return df


def findfpa(file=None, dataframe=None, savename='fpadefault.csv'):
    """
    A function that will generate results similar to those reported in Kulesa et al (2018), specifically by generated
    to_csv(csv_out_txt) histograms of the Probability of entering a given steady state. This function can take an import
    from previously generated data stored a .csv
    :param file: A .csv file generated using the datahandler() to be passed into this function for analysis
    :return: A sorted object  and an excel file of the same data
    """

    if file != None:
        df = pd.read_csv(file)
    elif dataframe != None:
        df = dataframe
    fpa = df[df.CycleLength == 1]
    fpa["SteadyState"] = fpa['SteadyState'].str.strip('[]').astype(int)
    counts = fpa['SteadyState'].value_counts().sort_index()
    print counts
    counts.to_csv(savename)

    return counts


def findlc(file=None, dataframe=None, savename='lcdefault.csv'):
    """
    A function that generates a list of the unique limit cycles in a model. Takes output of datahandler
    :param dataframe: A dataframe generated in the all_ss_w_model_fix() to be passed into this function for analysis.
    :param file: A .csv file generated using the all_ss_w_model_fix() to be passed into this function for analysis.
    :return: Counts variable
    """
    if file is not None:
        df = pd.read_csv(file)
    elif dataframe is not None:
        df = dataframe
    lc = df[df.CycleLength >= 2]
    lc['SteadyState'] = lc['SteadyState'].str.strip('[]').astype(int)
    counts = lc['SteadyState'].value_counts().sort_index()
    counts.to_csv(savename)

    return counts


# An incomplete
# weights = np.ones_like(fpa.CycleLength) / float(len(df.CycleLength))
# # the histogram of the data
# plt.hist(x=fpa.SteadyState, bins=10, align='mid')
# plt.xlabel('SteadyState Number')
# plt.ylabel('Proportion in State')
# plt.title('Histogram')
# plt.show()









