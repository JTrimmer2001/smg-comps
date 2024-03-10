import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
from astropy.table import Table as tb
from astropy.cosmology import FlatLambdaCDM
import generictrawler as gt
import os
import pandas as pd
import math
import astrofuncs as af

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
    '''cats = gt.trawler('mastercats')

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
        plt.close('all')'''
    
    alldata = pd.read_csv('master_beamlim.csv')
    limdata = pd.read_csv('master_doubles.csv')
    path = 'plots/radec/doubles_testing/'
    #### INIT File paths and strings ####

    imglist = pd.unique(limdata['source'])

    for img in imglist:
        imgAll = alldata[alldata['source'] == img]
        limAll = limdata[limdata['source'] == img]

        limWindow = pd.unique(limAll['window'])

        if os.path.isdir(path+img+'/') == False:
            os.makedirs(path+img+'/')

        for window in limWindow:

            imgWin = imgAll[imgAll['window'] == window]
            limWin = limAll[limAll['window'] == window]

            '''mask0_6 = (limWin['F_kw']>0.6)
            mask0_8 = (limWin['F_kw']>0.8)
            fidelity = limWin[mask0_6]
            mask2 = (fidelity['F_kw']<=0.8)
            fidelity0_6 = fidelity[mask2]
            fidelity0_8 = limWin[mask0_8]'''###Use when using a full set of data

            fig, ax = plt.subplots()

            '''f00=ax.scatter(imgWin['RA'], imgWin['DEC'],s=0.5,c='0.8',label='Raw data')
            f06=ax.scatter(fidelity0_6['RA'],fidelity0_6['DEC'],marker='2',c='r',label='0.6 < F_kw <= 0.8')
            f08=ax.scatter(fidelity0_8['RA'],fidelity0_8['DEC'],marker='1',c='b',label='F_kw > 0.8')
            ax.legend(handles=[f00,f06,f08])'''###Use this when doing plots with a full sample of raw data

            raw = ax.scatter(imgWin['RA'], imgWin['DEC'], s=0.5,c='r', label = 'All Data')
            data = ax.scatter(limWin['RA'], limWin['DEC'],marker='1',c='b',label='Duplicates')
            ax.legend(handles=[raw,data])

            ax.set(xlabel='RA',ylabel='DEC')
            

            fig.savefig(path+img+'/'+window+'.pdf',bbox_inches='tight')
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

    data = pd.read_csv('sepMasterBeamLim.csv') #Loads entire bank of data
    imgList = pd.unique(data['refObj']) #Gets a list of available images
    path = 'plots/separation/beamlim/FDbins/'

    for img in imgList:
        refObjtbl = data[data['refObj']==img]

        if os.path.isdir(path) == False:
            os.makedirs(path)

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

        fig.savefig(path+img+'.pdf',bbox_inches='tight')
        plt.close('all')

def luminosityFunction():
    table = pd.read_csv('matched_err1.0_with_lines.csv')
    index = list(table.index.values)
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    hist2_1 = np.array([0,0,0,0],dtype=np.float64)
    hist3_2 = np.array([0,0,0,0],dtype=np.float64)
    hist4_3 = np.array([0,0,0,0],dtype=np.float64)

    candidateStats = {}

    for item in index:
        z = table.at[item, 'zphot']
        D_l = cosmo.luminosity_distance(z=z).value
        fidelity = table.at[item,'F']
        line = str()

        if table.at[item,'CO2-1'] == True:
            line = '2-1'
            freq_obs = 230.538
        elif table.at[item,'CO3-2'] == True:
            line = '3-2'
            freq_obs = 345.796
        elif table.at[item,'CO4-3'] == True:
            line = '4-3'
            freq_obs = 461.041
        else:
            continue

        F = table.at[item, 'f_line']

        L = ((3.25e+7)/(1+z)**3)*F*(freq_obs**(-2))*(D_l**2)

        if line == '2-1':
            hist, bins = af.dexHistogram(7,12,L)
            hist = fidelity * hist
            hist2_1 += hist
        elif line == '3-2':
            hist, bins = af.dexHistogram(7,12,L)
            hist = fidelity * hist
            hist3_2 += hist
        elif line == '4-3':
            hist, bins = af.dexHistogram(7,12,L)
            hist = fidelity * hist
            hist4_3 += hist 
        else:
            continue
    
    v2_1 = af.obsVolume(z1=1.3954,z2=1.5125,f1=91.756e+9,f2=96.243e+9,factor=np.sqrt(2),dish=12)
    v3_2 = af.obsVolume(z1=2.593,z2=2.7686,f1=91.756e+9,f2=96.243e+9,factor=np.sqrt(2),dish=12)
    v4_3 = af.obsVolume(z1=3.79,z2=4.025,f1=91.756e+9,f2=96.243e+9,factor=np.sqrt(2),dish=12)
    vs = [v2_1,v3_2,v4_3]
    hists = [hist2_1,hist3_2,hist4_3]

    '''f, axs = plt.subplots(1,3,sharey=True)

    x, y = zip(*co2_1)
    axs[0].scatter(x,y,label='CO2-1')
    #axs[0].xlabel('luminosity')
    #axs[0].ylabel('Density Mpc^-3')

    x, y = zip(*co3_2)
    axs[1].scatter(x,y,label='CO3-2')
    #axs[1].xlabel('luminosity')
    #axs[1].ylabel('Density Mpc^-3')

    x, y = zip(*co4_3)
    axs[2].scatter(x,y,label='CO4-3')
    #axs[2].xlabel('luminosity')
    #axs[2].ylabel('Density Mpc^-3')

    for ax in axs:
        ax.set(xlabel='Luminosity',ylabel='Density Mpc^-3')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(7,12)
        ax.set_ylim(-6,-0.5)
        ax.label_outer()
        ax.legend()

    f.tight_layout()
    plt.show()'''

    f, axs = plt.subplots(1,3,sharey=True)

    for c in range(3):
        v = vs[c]
        hist = hists[c]

        for i in range(4):
            xmin = bins[i]
            xmax = bins[i+1]
            width = xmax - xmin

            value = hist[i]
            lf = (1/v)*value

            err = (1/v)*np.sqrt(value)
            ymin = lf - err
            height = 2*err
            xy = (xmin, ymin)

            rect = patches.Rectangle(xy,width,height)

            axs[c].add_patch(rect)

    for ax in axs:
        ax.set(xlabel='Luminosity',ylabel='Density Mpc^-3')
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlim(10**7,10**12)
        ax.set_ylim(10**(-6),10**(-0.5))
        ax.label_outer()
        ax.legend()

    f.tight_layout()
    plt.show()




luminosityFunction()