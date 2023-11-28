import astropy.table as Table
import numpy as np
import matplotlib.pyplot as plt
import generictrawler as trawl

cats = trawl.trawler('mastercats')
tablemade = False

for f in cats:
    parts = f.split('_')
    source = parts[0]
    table = Table.read('mastercats/'+f,format='ascii')

    table['source'] = source

    if tablemade == False:
        master = table
    else:
        master = Table.vstack([master,table])

Table.write('master.dat',format='ascii')


