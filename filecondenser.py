from astropy.table import Table as tb
from astropy.table import vstack
import numpy as np
import matplotlib.pyplot as plt
import generictrawler as trawl
import pandas as pd

cats = trawl.trawler('mastercats')
tablemade = False

for f in cats:
    parts = f.split('_')
    source = parts[0]
    window = parts[1]
    table = tb.read('mastercats/'+f,format='ascii')

    table['source'] = source
    table['window'] = window

    if tablemade == False:
        master = table
        tablemade=True
    else:
        master = vstack([master,table])

df = master.to_pandas()

df.to_csv('master.csv')
