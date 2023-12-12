import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table as tb
import generictrawler as gt
import os

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
        

snrplotter()
