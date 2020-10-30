# Predicting-Neuroblastoma-Logic-Based-Model
This repository contains a function to allow for the recreation of the model from the paper 
"Predicting neuroblastoma using developmental signals and a logic-based model" 
(_Biophys Chem_ **238** (2018), 30-38)  If you use our computer code, please cite our paper.  For 
more details about the paper, please visit https://doi.org/10.1016/j.bpc.2018.04.004

Notes about the computer code
-----------------------------

The datahandler() function locates in /Files/BooleanFunct.py can be called by working .py documents in the 
/Files/ folder to simulate a boolean network model whose rules are specified by the user in a .txt file.
Fixed point attractors and limit cycles can be found using the lc() and fpa() functions in Booleanfunct.py.
There is an option when using the datahandler() to print out a graphical representation of the on proportion
of four nodes of interest. NOTE. At this time, the user must specify four nodes to use this feature, even if 
there are fewer than four nodes that the researcher is interested in.

In the future, I hope to automate the processing of the datahandler data into graphs like those in the paper, 
specifically that show specific fixed point attractor by the likelyhood of entering that attractor starting 
from a random system state

To process the data gathered by the datahandler, manual data cleanup is conducted in excel. To identify all 
possible fixed point attractors in a simulated network, you should run the datahandler() with your specified 
parameters and open the resulting CSV file in excel. When the excel sheet is open, first remove any rows 
that have CycleLength > 1 using the filter data tool. Now, to 

The user can change the ending system state data to individual columns by using the text-to-columns excel 
tool and first delineating by : and not keeping the first column. And then doing the text-to-columns excel 
tool again this time separating by , and (space) and excluding the names of the nodes from the transformation.

To figure out how likely it was to enter a given steady fpa, use the fpa(). the first column is the SS number 
and the second is the number of occurances. Divide the number of occurances by the number of 2^(number of nodes) 
to get proportion on.

Combine this data and generate a histogram that refers to the specific steady states and their percentages, 
and make sure the steady state identities are tabulated and easily accessable for further analysis.

For more details, we strongly recommend reading Gribbin Will's Master in Physiology Capstone dissertation 
in root directory. 
