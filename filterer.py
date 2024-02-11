#from ctypes import windll
import pandas as pd
import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import SpectralCoord
from astropy.coordinates import SpectralQuantity
from astropy.coordinates import Angle
from astropy.io import fits
#from sympy import false
import generictrawler as gt
from astropy.wcs import WCS
from astroquery.splatalogue import Splatalogue

y=0

######################CHANGE THE FOLDER BEFORE RUNNING#########################

folder = str('F:/Data/') #Gets the suffix for the file address

table = pd.read_csv('matched_err1.0_with_lines.csv')

###############################################################################
def fidlim():
    tbl = pd.read_csv('master.csv')
    flim = tbl[tbl['F'] > 0.6] #Imposes a fidelity per kernel width limit of 0.6
    flim.to_csv('master_fidlim.csv',index=False)

def beamlimit():
    y=0
    sauces = pd.unique(table['source']) #Makes a list of the source images in the source column of the table

    dflist = []

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

            with fits.open(folder + file) as hdul:
                data = hdul[0].data
                hdr = hdul[0].header

            '''hdul = fits.open(folder +  file) #Opens image file for the window'''

            '''xc = hdul[0].header['CRVAL1']*u.deg
            yc = hdul[0].header['CRVAL2']*u.deg
    '''
            w = WCS(hdr, naxis=(1,2)) #SHOULD set the wcs settings for this file

            crx = hdr['CRPIX1']
            cry = hdr['CRPIX2']

            for source in range(length):
                coords = [window.iloc[source,3],window.iloc[source,4]] #Gets the image coordinates of the object

                f = window.iloc[source,2] * 1e9 # frequency source is located at

                '''fslice = data[:,:,f,1]'''

                wref = SkyCoord(w.pixel_to_world(crx,cry)) # skycoord of the central pixel
                clump_pos = SkyCoord(w.pixel_to_world(coords[0],coords[1])) #Gets world coord of clump

                wavel = 3e+8/f
                fwhm = 1.22 * (wavel/12) #Calculating beam radius at current f
                seplim = Angle(fwhm/np.sqrt(2), unit=u.rad)
                tester = seplim.arcsec
    
                sep = clump_pos.separation(wref) #finding distance from center
                if sep.radian <= seplim.radian: #determining if object is in the beam

                    y+=1
                    row = window.iloc[[source]]#selecting the current objects data
                    dflist.append(row)
     
            hdul.close()
            
    filtered = pd.concat(dflist)
    filtered.to_csv('master_beamlim.csv',index = False)
    print('sources in range:', y)   

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

    seps.to_csv('sepMasterBeamLim.csv',index = False)

def dupeCatcher():
    count = 0
    idlist = [] #List of IDs to be added to the database (for duplicate entry checking)
    masterlist = [] #List of objects to be added, will add dicts and construct from that

    images = pd.unique(table['source']) #List of target objects

    for img in images:
        imgTable = table[table['source'] == img] #Filters table down to current image
        windows = pd.unique(imgTable['window']) #generates a list of the spws used in this image

        for window in windows:

            winTable = imgTable[imgTable['window']==window]#Narrows table to objects within current window
            ids = pd.unique(winTable['ID']) #Gets list of IDs

            for i in ids:

                currentObj = winTable[winTable['ID']==i] #Single row of window table containing current primary object data

                primaryObj = {'ra':currentObj.iloc[0,0],
                              'dec':currentObj.iloc[0,1],
                              'freq':currentObj.iloc[0,2],
                              'id':currentObj.iloc[0,11]} #Dict organising relevant data into a easily accesible spot
                
                primarySC = SkyCoord(primaryObj['ra'],primaryObj['dec'], unit=u.deg) 
                #skycoord object for primary object, allows for the separation 
                
                for s in ids:

                    if s == i: #skips same id
                        continue

                    currentObj = winTable[winTable['ID']==s]

                    secondaryObj = {'ra':currentObj.iloc[0,0],
                                    'dec':currentObj.iloc[0,1],
                                    'freq':currentObj.iloc[0,2],
                                    'id':currentObj.iloc[0,11]}
                    
                    secondarySC = SkyCoord(secondaryObj['ra'],secondaryObj['dec'], unit=u.deg)

                    sepAngle = primarySC.separation(secondarySC)

                    if sepAngle.arcsec <=3:
                        if s in idlist:
                            continue
                        else:
                            #idlist.append(s)
                            count+=1

    mask = table['ID'].isin(idlist)
    dupes = table[mask]

    dupes.to_csv('duplicates_fixed.csv',index=False)
    
           
def ClashOfTheClumps():
    Ids = []
    
    knockouts = []

    sourceFile = pd.read_csv('duplicates_fixed.csv',index_col='ID')
    sourceImgList = pd.unique(sourceFile['source'])

    for source in sourceImgList:
        sourceImg = sourceFile[sourceFile['source']==source]
        windowList = pd.unique(sourceImg['window'])

        for window in windowList:
            img = sourceImg[sourceImg['window']==window]

            windowIds = img.index.values.tolist()

            for i in windowIds:
                if i in knockouts:
                    continue

                primaryObj = {'ra':img.loc[i,'RA'],
                              'dec':img.loc[i,'DEC'],
                              'f_kw':img.loc[i,'F_kw'],
                              'snr':img.loc[i,'SNR']}
                
                primarySC = SkyCoord(primaryObj['ra'],primaryObj['dec'], unit=u.deg)
                
                for s in windowIds:

                    if s in knockouts:
                        continue
                    elif s == i:
                        continue

                    secondaryObj = {'ra':img.loc[s,'RA'],
                                    'dec':img.loc[s,'DEC'],
                                    'f_kw':img.loc[s,'F_kw'],
                                    'snr':img.loc[s,'SNR']}
                    
                    secondarySC = SkyCoord(secondaryObj['ra'],secondaryObj['dec'],unit=u.deg)

                    sep = primarySC.separation(secondarySC)

                    if sep.arcsec > 1:
                        continue
                    elif sep.arcsec <= 1:
                        if primaryObj['f_kw']>secondaryObj['f_kw']:
                            knockouts.append(s)
                        elif primaryObj['f_kw']<secondaryObj['f_kw']:
                            knockouts.append(i)
                        elif primaryObj['f_kw']==secondaryObj['f_kw']:
                            if primaryObj['snr']>secondaryObj['snr']:
                                knockouts.append(s)
                            elif primaryObj['snr']<secondaryObj['snr']:
                                knockouts.append(i)
                            else:
                                knockouts.append(s)

            for id in windowIds:

                if id in knockouts:
                    continue
                else:
                    Ids.append(id)

    print(Ids)

    mask = ~table['ID'].isin(knockouts)
    out = table[mask]

    out.to_csv('master_doubles.csv',index=False)

def lineIdentifier():
    '''Takes in a set of values with assosciated high and low redshifts 
        and converts to a rest frame wavelength, before comparing to a 
        catalogue of line emission frequencies. Built on the photometric
        data catalogue used for the masters project in 2024.'''
    ids = list(table.index.values) #Gets list of indexes for iterating through

    for i in ids:
        specinfo = {'freq_ghz':table.at[i,'FREQ_GHZ'],
                    'zhi'     :table.at[i,'zhi'],
                    'zlo'     :table.at[i,'zlo'],
                    'zphot'   :table.at[i,'zphot']} #Dictionary of hi and lo z, plus the line freq
        
        hiZcoord = SpectralCoord(value=specinfo['freq_ghz'],unit=u.GHz,redshift=specinfo['zhi'])
        hiRestCoord = hiZcoord.to_rest()
        hiRest = hiRestCoord.value*u.GHz
        loZcoord = SpectralCoord(value=specinfo['freq_ghz'],unit=u.GHz,redshift=specinfo['zlo'])
        loRestCoord = loZcoord.to_rest() #gets rest frame hi and lo estimates of emission freq
        loRest = loRestCoord.value*u.GHz
        photZcoord = SpectralCoord(value=specinfo['freq_ghz'],unit=u.GHz,redshift=specinfo['zphot'])
        midRest = photZcoord.to_rest() #Same thing but for best guess Z

        length_range = (hiRest,loRest)

        colines = Splatalogue.query_lines(min_frequency=hiRest,max_frequency=loRest,chemical_name=' HCO+ ',only_astronomically_observed=True)
        colines.keep_columns(['Species','Freq-GHz(rest frame,redshifted)','Resolved QNs','Lovas/AST Intensity'])
        
        '''wavelength_range = (hiRestCoord.value * u.GHz,loRestCoord.value * u.GHz)

        lineList = AtomicLineList.query_object(wavelength_range=wavelength_range, wavelength_type='Air',
                            wavelength_accuracy=20, element_spectrum='C II-IV',
                                                output_columns=('spec','term','prob'))'''
        
        print('Frequencies used: {}'.format(length_range))
        colines.pprint_all(max_lines=-1)
        existing = table.iat[i,61]
        
        transition = str(input())
        

        if transition == '':
            continue
        else:
            if pd.isna(table.at[i,'transition']) == True:
                table.at[i,'transition']=transition
            else:
                table.at[i,'transition']=str(transition+' or '+existing)
            


    table.to_csv('matched_err1.0_with_lines.csv',index=False)


lineIdentifier()
#ClashOfTheClumps()
#dupeCatcher()
