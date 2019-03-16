from boolean2 import Model, state, tokenizer, util
import pandas as pd
import csv

def dict_to_state_list(dict):
    """
    A function that converts a dictionary of nodes and there values to two lists.
    :param dict: Input a list of dictionaries to convert the data to 01010101 form
    :return: data: will return a list of gapless binary data
             label: will return a tuple containing the names of the nodes in the same
             order as the binary data.
    """
    data = []
    for i in range(len(dict)):
        a = dict[i].values()
        for i in range(len(a)):
            if a[i] == True:
                a[i] = 1
            elif a[i] == False:
                a[i] = 0
            else:
                print 'Error: Input must be as a list of dict'
        c = str(a)
        c = str(c.translate(None, '. ,[]'))
        data.append(c)
    labels = tuple(dict[1].keys())
    return data, labels

def all_ss_model_w_fix(str_file, str_mode, timecourse, csv_out_txt, geneon=[], geneoff=[], dumpevery=1000,\
                       nodes_for_averages=['Apoptosis', 'Proliferation', 'Angiogenesis', 'Differentiation']):
    '''
    A function to generate a list of all possible starting states of a model. Note that this function is a
    modification of the all_ss_model() function and allows for specific nodes to be fixed on or off.
    str_file: as a string, write the name of the file containing the rules for the boolean model
    str_mode: enter the mode of model that is going to be used, ie sync or async
    numnodes: the number of nodes in the model
    :return:
    A list of dictionaries with all possible starting states of the network.
    '''

    def dumper():
        '''
        will dump the data into a CSV file of name csv_out_txt
        :return: none
        '''
        with open(csv_out_txt, 'a') as csvfile:
            datadumper = csv.writer(csvfile, lineterminator='\n')
            for row in d:
                datadumper.writerow(row)
        return

    # header names for pandas dataframe and CSV data dump file
    headers = ['Index', 'CycleLength', 'CycleFingerprint', 'FirstState', 'LastState', str(nodes_for_averages[0]),\
               str(nodes_for_averages[1]), str(nodes_for_averages[2]), str(nodes_for_averages[3])]

    # Init a CSV to contain the data generated over the course of the run
    with  open(csv_out_txt,'wb') as out:
        csv_out = csv.writer(out)
        csv_out.writerow(headers)


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

    if geneoff != [] or geneon != []:
        d = []
        p = 0
        n = 0
        for data_init in initializer:
            data = data_init[0]
            for i in range(len(geneoff)):
                if data[str(geneoff[i])]:
                    initializer_new.append(data_init)
            for i in range(len(geneon)):
                if data[str(geneon[i])]:
                    initializer_new.append(data_init)

        for data, initfunc in initializer_new:
            # Fixes genes on or off (True or False) at their starting state.
            for i in range(len(geneoff)):
                data.update({str(geneoff[i]): False})
            for i in range(len(geneon)):
                data.update({str(geneon[i]): True})

            # Initialize the model with the given pre-computed initial conditions stored in the data variable.
            # for a given model with Z nodes, there are 2 to the power of Z starting states.
            model.initialize(defaults=data)
            model.iterate(steps=timecourse)
            e = model.detect_cycles()
            nodes = nodes_for_averages

            d.append(list(model.detect_cycles()))
            d[p].append(model.fp())
            d[p].append(model.first)
            d[p].append(model.last)
            detect_states = 0 - e[1]
            coll.collect(states=model.states[detect_states:], nodes=nodes)
            # Console output for debugging and progress tracking
            print 'The fingerprint is', d[p][2]
            print 'The cycle length is', d[p][1]

            prc = ((n) / float(pow(2, load_status_divisor))) * 100
            print '%.2f' % prc + "% done."
            print '\n'

            # Iterate the model
            n += 1
            p += 1

            if p == dumpevery:
                dumper()
                p = 0
                d = []

        dumper()
        d = []
    else:
        d = []
        p = 0
        n = 0
        for data, initfunc in initializer:
            # Fixes genes on or off (True or False) at their starting state.
            for i in range(len(geneoff)):
                data.update({str(geneoff[i]): False})
            for i in range(len(geneon)):
                data.update({str(geneon[i]): True})

            # Initialize the model with the given pre-computed initial conditions stored in the data variable.
            # for a given model with Z nodes, there are 2 to the power of Z starting states.
            model.initialize(defaults=data)
            model.iterate(steps=timecourse)
            e = model.detect_cycles()
            nodes = nodes_for_averages

            d.append(list(model.detect_cycles()))
            d[p].append(model.fp())
            d[p].append(model.first)
            d[p].append(model.last)
            detect_states = 0 - e[1]
            coll.collect(states=model.states[detect_states:], nodes=nodes)

            # Console output for debugging and progress tracking
            print 'The fingerprint is', d[p][2]
            print 'The cycle length is', d[p][1]

            prc = ((n) / float(pow(2, load_status_divisor))) * 100
            print '%.2f' % prc + "% done."
            print '\n'

            # Iterate the model
            n += 1
            p += 1

            if p == dumpevery:
                dumper()
                p = 0
                d = []

        dumper()
        d = []

    # Generate a pandas dataframe of the dump file

    df = pd.read_csv(csv_out_txt)

    # Add the avg on state of the nodes specified by the user. They will output to new columns in row 1 of the data
    if nodes_for_averages != []:
        avgs = coll.get_averages(normalize=True)
        a1 = avgs[nodes_for_averages[0]]
        b1 = avgs[nodes_for_averages[1]]
        c1 = avgs[nodes_for_averages[2]]
        d1 = avgs[nodes_for_averages[3]]

        # Set the averages values in the dataframe.

        df.set_value(0, 'Apoptosis', a1[0])
        df.set_value(0, 'Proliferation', b1[0])
        df.set_value(0, 'Differentiation', d1[0])
        df.set_value(0, 'Angiogenesis', c1[0])

        df.to_csv(csv_out_txt)

    return df



