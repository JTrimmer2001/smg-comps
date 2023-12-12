from astropy.table import Table as tb
import numpy as np
import matplotlib.pyplot as plt
import os

table = tb.read('master.csv')
mask06 = (table['F_kw']>0.6)
mask08 = (table['F_kw']>0.8)
fidelity0_6 = table[mask06]
fidelity0_8 = table[mask08]

fig, ax = plt.subplots()
f00=ax.scatter(table['X'], table['Y'],s=0.5,c='0.8',label='F_kw > 0')
f06=ax.scatter(fidelity0_6['X'],fidelity0_6['Y'],marker='2',c='r',label='F_kw > 0.6')
f08=ax.scatter(fidelity0_8['X'],fidelity0_8['Y'],marker='1',c='b',label='F_kw > 0.8')

ax.set(xlabel='Image X',ylabel='Image Y')
ax.legend(handles=[f00,f06,f08])

if os.path.isdir('plots/special/combi/') == False:
    os.makedirs('plots/special/combi/')

fig.savefig('plots/special/combi/XYfid.pdf',bbox_inches='tight')