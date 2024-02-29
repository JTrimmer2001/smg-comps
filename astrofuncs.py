import numpy as np
import math
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

######### Cosmology details ###############
cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

def obsVolume(z1,z2,f1,f2,dish,factor=None):
    '''Function to find the observed volume of a telescope pointing.
    Assumes a flat lambda cdm model
    Uses 1/3*pi*(r2**2*h + r2**2*D - r1**2*h) = V
    
    Inputs:
        z1 (float): Lower redshift limit
        z2 (float): Upper redshift limit
        f1 (float): Lower frequency limit
        f2 (float): Upper frequency limit
        dish (float): diameter of receiving instrument
        factor (float): Factor used to modify radius of observations
        
    Returns:
        Volume (float): Total observed volume for the given freqeuency and redshifts in Mpc^-3'''
    
    ###### Getting relevant measurements ######
    wavel1 = 3e+8 / f1
    fwhm1 = 1.22 * (wavel1/dish)
    r1 = fwhm1/2
    
    wavel2 = 3e+8/f2
    fwhm2 = 1.22 * (wavel2/dish)
    r2 = fwhm2/2

    conversion1=cosmo.kpc_comoving_per_arcmin(z1).to(u.Mpc/u.arcmin)
    conversion2 = cosmo.kpc_comoving_per_arcmin(z2).to(u.Mpc/u.arcmin)

    r1 = r1 * conversion1.value
    r2 = r2 * conversion2.value # Converts angular radius to Mpc radius at z1 and z2

    if factor != None:
        A1 = math.pi*factor*(r1**2)
        A2 = math.pi*factor*(r2**2)
    if factor == None:
        A1 = math.pi*(r1**2)
        A2 = math.pi*(r2**2)

    h = cosmo.lookback_distance(z1)
    D = cosmo.lookback_distance(z2)

    volume = (1/3)*(A2 * h + A2 * D - A1*h)

    return volume.value

    
def deltaBeamsize(bmaj,bmin,rmaj,rmin):
    beamArea = math.pi * bmaj/2 * bmin/2
    regArea = math.pi * rmaj * rmin

    factor = regArea/beamArea

    return factor
    
def dexHistogram(min,max,data):
    '''Method to get a "per dex" visualisation of data
    
    Inputs: min (int)  :Minimum power of 10 to use for the histogram
            max (int)  :Maximum power of 10
            data(list) :List of data to be used for the histogram
            
    Outputs: hist (npy histogram):  values of histogram
                                    Has a values and bin edge component'''
    
    binEdges = []

    for power in range(min,max):
        binEdges.append(10**power)

    hist = np.histogram(data,bins=binEdges)

    return hist



    
