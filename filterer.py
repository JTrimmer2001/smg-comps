import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.io import fits
import generictrawler as gt
import astropy.wcs as wcs

y=0

table = pd.read_csv('master_fidlim.csv')

'''flim = table[table['F_kw'] > 0.6] #Imposes a fidelity per kernel width limit of 0.6
flim.to_csv('master_fidlim.csv')'''

sauces = pd.unique(table['source'])

for i in sauces:
    img = table[table['source'] == i]

    winlist = pd.unique(table['window'])

    for spw in winlist:
        window = img[img['window']==spw]

        length = len(window.index)
        name = str(i + '_dirtycube_' + spw)

        file = gt.trawler(path='E:/ExtHDDBackup/Data/',kword=name) #Finding and opening relevant fits file
        if file == 'none':
            continue

        hdul = fits.open('E:/ExtHDDBackup/Data/' +  file)

        xc = hdul[0].header['CRVAL1']*u.deg
        yc = hdul[0].header['CRVAL2']*u.deg

        wref = SkyCoord(xc,yc) # skycoord of the central pixel

        for source in range(length):
            coords = [window.iloc[source,0],window.iloc[source,1]]
            clump_pos = SkyCoord(coords[0]*u.deg,coords[1]*u.deg) #Gets world coord of clump

            f = window.iloc[source,2] * 1e9 # frequency source is located at
            wavel = 3e+8/f
            fwhm = 1.13 * (wavel/11)
            seplim = Angle(fwhm/np.sqrt(2), unit=u.rad)
 
            sep = clump_pos.separation(wref)
            if sep.radian <= seplim.radian:
                y+=1
                print('sources in range:', y)          



            
