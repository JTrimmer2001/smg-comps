import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import Angle
from astropy.io import fits
import generictrawler as gt
from astropy.wcs import WCS

y=0

######################CHANGE THE FOLDER BEFORE RUNNING#########################

folder = str('D:/Data/') #Gets the suffix for the file address

table = pd.read_csv('master_fidlim.csv')

###############################################################################

'''flim = table[table['F_kw'] > 0.6] #Imposes a fidelity per kernel width limit of 0.6
flim.to_csv('master_fidlim.csv')'''

def beamlimit():
    sauces = pd.unique(table['source']) #Makes a list of the source images in the source column of the table

    for i in sauces: #Iterates through the images
        img = table[table['source'] == i] #Gets a list of objects in the image

        winlist = pd.unique(table['window']) # makes a list of the spectral windows 

        for spw in winlist:
            window = img[img['window']==spw] #Gets a table of the objects in the provided window

            length = len(window.index) #Gets the number of objects in the current window
            name = str(i + '_dirtycube_' + spw) #Produces the name of the object so image files can be located

            file = gt.trawler(path=folder,kword=name) #Finding and opening relevant fits file
            if file == 'none':
                continue #Skips iteration if no files are found, this shouldn't occur

            hdul = fits.open(folder +  file) #Opens image file for the window

            '''xc = hdul[0].header['CRVAL1']*u.deg
            yc = hdul[0].header['CRVAL2']*u.deg
    '''
            w = WCS(hdul[0].header) #SHOULD set the wcs settings for this file

            for source in range(length):
                coords = [window.iloc[source,3],window.iloc[source,4]] #Gets the image coordinates of the object

                f = window.iloc[source,2] * 1e9 # frequency source is located at

                wref = SkyCoord(w.pixel_to_world(513,513,f,1)) # skycoord of the central pixel
                clump_pos = SkyCoord(w.pixel_to_world(coords[0],coords[1],f,1)) #Gets world coord of clump

                wavel = 3e+8/f
                fwhm = 1.22 * (wavel/12)
                seplim = Angle(fwhm/np.sqrt(2), unit=u.rad)
    
                sep = clump_pos.separation(wref)
                if sep.radian <= seplim.radian:
                    y+=1
                    print('sources in range:', y)       

            hdul.close()   


def separations():

    masterList = []

    images = pd.unique(table['source']) # List of all unique source images

    for img in images:
        currentObjects = table[table['source'] == img] #Table of objects within the current image
                                                       # All table data is preserved

        currentObjectsNum = len(currentObjects.index) #Gets the number of objects currently being
                                                      #studied
        
        for p in range(currentObjectsNum-1):

            primaryObj = {'ra':  currentObjects.iloc[p,0],
                          'dec': currentObjects.iloc[p,1],
                          'freq':currentObjects.iloc[p,2]
                         }
                #^^ Creates a list of characteristics for the current
                #primary object. Follows the format of RA, DEC, F for
                #use when studying the separation of objects

            primarySC = SkyCoord(primaryObj['ra'],primaryObj['dec'], unit=u.deg)

            for s in range(currentObjectsNum-1):

                if s == p:
                    continue #Catch case for no duplicates

            
                secondaryObj = {'ra':  currentObjects.iloc[s,0],
                                'dec': currentObjects.iloc[s,1],
                                'freq':currentObjects.iloc[s,2]
                                }
                    #^^ Creates a list of characteristics for the current
                    #primary object. Follows the format of RA, DEC, F for
                    #use when studying the separation of objects

                secondarySC = SkyCoord(secondaryObj['ra'],secondaryObj['dec'], unit=u.deg)

                ### Finding separations ###
                
                sepAngle = primarySC.separation(secondarySC)
                dFrequency = abs(secondaryObj['freq']-primaryObj['freq'])

                sublist = [sepAngle.degree, dFrequency, img]
                masterList.append(sublist)

    seps = pd.DataFrame(masterList)  #Produces a data frame and saves as a csv

    seps.to_csv('separationsMaster.csv')


separations()           
