import pandas as pd
import numpy as np

table = pd.read_csv('master.csv')

flim = table[table['F_kw'] > 0.6] #Imposes a fidelity per kernel width limit of 0.6
flim.to_csv('master_fidlim.csv')