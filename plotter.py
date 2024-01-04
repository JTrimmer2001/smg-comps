import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table as tb
import generictrawler as gt
import os
import pandas as pd
import math

def fplotter():
    cats = gt.trawler('mastercats')

    for i in cats:
        parts = i.split('_')
        name = 'mastercats/'+i
        data = tb.read(name, format='ascii')
        mask = (data['F_kw']>0)
        fidelity = data[mask]

        if os.path.isdir('plots/fidelity/'+parts[0]+'/') == False:
            os.makedirs('plots/fidelity/'+parts[0]+'/')

        plt.hist(fidelity['F_kw'],bins=40)
        plt.savefig('plots/fidelity/'+parts[0]+'/'+parts[1]+'.png')
        plt.xlabel('F_kw')
        plt.clf()

def fvsnr():

    cats = gt.trawler('mastercats')
    for i in cats:
        parts = i.split('_')
        name = 'mastercats/'+i
        data = tb.read(name, format='ascii')

        if os.path.isdir('plots/fvsnr/'+parts[0]+'/') == False:
            os.makedirs('plots/fvsnr/'+parts[0]+'/')

        mask = (data['F_kw']>0)&(data['SNR']>0)
        points = data[mask]

        plt.scatter(points['SNR'],points['F_kw'])
        plt.xlabel('SNR')
        plt.ylabel('F_kw')

        plt.savefig('plots/fvsnr/'+parts[0]+'/'+parts[1]+'.png')
        plt.clf()

def snrplotter():
    cats = gt.trawler('mastercats')

    for i in cats:
        parts = i.split('_')
        name = 'mastercats/'+i
        data = tb.read(name, format='ascii')
        mask = (data['SNR']>0)
        fidelity = data[mask]

        if os.path.isdir('plots/SNRhist/logy/'+parts[0]+'/') == False:
            os.makedirs('plots/SNRhist/logy/'+parts[0]+'/')

        plt.hist(fidelity['SNR'],bins=40)
        plt.savefig('plots/SNRhist/logy/'+parts[0]+'/'+parts[1]+'.png')
        plt.xlabel('SNR')
        plt.yscale('log')
        plt.clf()
        

def radec():
    cats = gt.trawler('mastercats')

    for i in cats:
        parts = i.split('_')
        name = 'mastercats/'+i
        data = tb.read(name, format='ascii')
        #mask = (data['RA']>0)&(data['DEC']>0)
        #fidelity = data[mask]

        if os.path.isdir('plots/radec/plain/'+parts[0]+'/') == False:
            os.makedirs('plots/radec/plain/'+parts[0]+'/')

        plt.scatter(data['RA'], data['DEC'],s=1)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.savefig('plots/radec/plain/'+parts[0]+'/'+parts[1]+'.png',bbox_inches='tight')
        plt.clf()

def radecFidelity():
    cats = gt.trawler('mastercats')

    for i in cats:
        parts = i.split('_')
        name = 'mastercats/'+i
        data = tb.read(name, format='ascii')
        mask0_6 = (data['F_kw']>0.6)
        mask0_8 = (data['F_kw']>0.8)
        fidelity = data[mask0_6]
        mask2 = (fidelity['F_kw']<=0.8)
        fidelity0_6 = fidelity[mask2]
        fidelity0_8 = data[mask0_8]

        if os.path.isdir('plots/radec/betterwFkw/'+parts[0]+'/') == False:
            os.makedirs('plots/radec/betterwFkw/'+parts[0]+'/')

        fig, ax = plt.subplots()

        f00=ax.scatter(data['RA'], data['DEC'],s=0.5,c='0.8',label='F_kw > 0')
        f06=ax.scatter(fidelity0_6['RA'],fidelity0_6['DEC'],marker='2',c='r',label='F_kw > 0.6')
        f08=ax.scatter(fidelity0_8['RA'],fidelity0_8['DEC'],marker='1',c='b',label='F_kw > 0.8')

        ax.set(xlabel='RA',ylabel='DEC')
        ax.legend(handles=[f00,f06,f08])

        fig.savefig('plots/radec/betterwFkw/'+parts[0]+'/'+parts[1]+'.pdf',bbox_inches='tight')
        plt.close('all')

def special():
        data = tb.read('mastercats/aless62_spw2123_clumpsP_minSNR_3.0_cropped_wfidelity_sorted.cat', format='ascii')
        mask = (data['FREQ_GHZ']>97)&(data['FREQ_GHZ']<98)
        data = data[mask]
        mask0_6 = (data['F_kw']>0.6)
        mask0_8 = (data['F_kw']>0.8)
        fidelity0_6 = data[mask0_6]
        fidelity0_8 = data[mask0_8]

        fig, ax = plt.subplots()

        f00=ax.scatter(data['RA'], data['DEC'],s=0.5,c='0.8',label='F_kw > 0')
        f06=ax.scatter(fidelity0_6['RA'],fidelity0_6['DEC'],marker='2',c='r',label='F_kw > 0.6')
        f08=ax.scatter(fidelity0_8['RA'],fidelity0_8['DEC'],marker='1',c='b',label='F_kw > 0.8')

        ax.set(xlabel='RA',ylabel='DEC')
        ax.legend(handles=[f00,f06,f08])

        if os.path.isdir('plots/special/') == False:
            os.makedirs('plots/special/')
            
        fig.savefig('plots/special/aless62_spw2123_radec_freqlimited9798.pdf',bbox_inches='tight')
        plt.close('all')

def sepAngleHist():

    data = pd.read_csv('separationsMaster.csv') #Loads entire bank of data
    imgList = pd.unique(data['refObj']) #Gets a list of available images

    for img in imgList:
        refObjtbl = data[data['refObj']==img]

        if os.path.isdir('plots/separation/FDbins/') == False:
            os.makedirs('plots/separation/FDbins/')

        n = refObjtbl['sepAngle'].size
        q1 = refObjtbl['sepAngle'].quantile(0.25)
        q3 = refObjtbl['sepAngle'].quantile(0.75)
        iqr = q3-q1

        max = refObjtbl['sepAngle'].max()
        min = refObjtbl['sepAngle'].min()
        binw = 2*(iqr * (n)**(-1/3)) #Freedman-Diaconis rule: https://en.wikipedia.org/wiki/Freedman%E2%80%93Diaconis_rule
        bins = math.ceil((max-min)/binw)

        fig, ax = plt.subplots()

        hist1 = ax.hist(refObjtbl['sepAngle'],bins=bins)
        ax.set(xlabel='Separation angle (degrees)',ylabel='N')

        fig.savefig('plots/separation/FDbins/'+img+'.pdf',bbox_inches='tight')
        plt.close('all')


sepAngleHist()
