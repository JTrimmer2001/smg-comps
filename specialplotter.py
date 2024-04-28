from astropy.table import Table as tb
import numpy as np
import matplotlib.pyplot as plt
import os
import pandas as pd
import astrofuncs as af

'''table = tb.read('master.csv')
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

fig.savefig('plots/special/combi/XYfid.pdf',bbox_inches='tight')'''


#new fidelity histogram/line graph

'''table = pd.read_csv('master.csv')
mask = (table['source']=='aless49')
fid = table[mask]
#fid = table[table['window']=='spw2123']

fig, ax = plt.subplots()
ax.hist(fid['F'],bins=20,range=(0,1))
ax.set_xlabel('Fidelity')
ax.set_yscale('log')

ax.plot(fid['SNR'],fid['F'])
ax.set_xlim(left=3,right=10)
ax.set_ylim(bottom=-0.1,top=1.1)
ax.set_ylabel('Fidelity')
ax.set_xlabel('SNR')
'''

# Separation histogram post-F catalogue

table = pd.read_csv('newSepsMaster.csv')
fig,ax = plt.subplots()
ax.hist(table['Separation'],bins=40,range=(0,10),histtype='step',color='k')
ax.set_xticks(range(0,11))
ax.set_xlim(left=0,right=10)
ax.set_yscale('log')
ax.set_xlabel('Separation (arcsec)')
ax.set_ylabel('Matches Log(N)')


af.tableformer(ax)
fig.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.show()