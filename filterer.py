import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.io import fits
import generictrawler as gt
from astropy.wcs import WCS

y=0

folder = str('D:/Data/')

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

        file = gt.trawler(path=folder,kword=name) #Finding and opening relevant fits file
        if file == 'none':
            continue

        hdul = fits.open(folder +  file)

        '''xc = hdul[0].header['CRVAL1']*u.deg
        yc = hdul[0].header['CRVAL2']*u.deg
'''
        w = WCS(hdul[0].header)

        for source in range(length):
            coords = [window.iloc[source,3],window.iloc[source,4]]

            f = window.iloc[source,2] * 1e9 # frequency source is located at

            wref = w.pixel_to_world([513,513,f,1],0) # skycoord of the central pixel
            clump_pos = w.pixel_to_world([coords[0],coords[1],f,1]) #Gets world coord of clump

            wavel = 3e+8/f
            fwhm = 1.22 * (wavel/12)
            seplim = Angle(fwhm/np.sqrt(2), unit=u.rad)
 
            sep = clump_pos.separation(wref)
            if sep.radian <= seplim.radian:
                y+=1
                print('sources in range:', y)       

        hdul.close()   



            
