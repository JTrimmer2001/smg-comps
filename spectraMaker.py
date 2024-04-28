from turtle import color
import comm
import matplotlib.pyplot as plt
import pandas as pd
import astrofuncs as af

df = pd.read_csv('FCatalogues/aless41_ID31.tsv',sep='\t',comment='#')

fig, ax = plt.subplots()

ax.plot(df['X'],df['Y']*1e3,c='k')
#ax.set_xlim(left=-44000,right=-38500)
ax.set_xlabel('Radio Velocity (km/s)')
ax.set_ylabel('Flux Density (mJy)')
ax.axhline(y=0,color='gray',linestyle='--')
ax.axvline(x=-4.14e4,color='b',linestyle='--')
ax.axvline(x=-4.02e4,color='b',linestyle='--')

af.tableformer(ax)
plt.show()