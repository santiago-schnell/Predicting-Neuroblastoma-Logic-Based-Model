import boolean2
from boolean2 import util, state, network

# updating rules

rules = """
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

# create the model
model = boolean2.Model( text=rules, mode='async')

# generates all states, set limit to a value to keep only the first that many states
# when limit is a number it will take the first that many initial states
initializer = state.all_initial_states( model.nodes, limit=None )

# the data is the inital data, the func is the initializer
for data, initfunc in initializer:
    # shows the initial values
    model.initialize(missing=initfunc)
    model.iterate(5)
