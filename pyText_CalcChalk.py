#!/usr/bin/python

import os
import glob
import matplotlib
matplotlib.use('Agg') #Allows for saving figures without producing a window
import matplotlib.pyplot as plt
from itertools import count
import numpy as np
from scipy.optimize import curve_fit
from pyTexture import guessDSpacing, getPeakCens, gauss, getHKL, fakePseudoVoigt, textResample, runQueues, feedQueue

np.seterr(divide='ignore')

if os.path.isfile('doPlots'):
    iPlot = True
else:
    iPlot = False

num_X = 2048
num_Y = 2048
num = 2000
areaScale = np.sqrt(np.log(2)/np.pi)
arcRad = 800
rhoStart, rhoEnd = -180, 180
rhoStep = 4
rhoN = 90
outDir = './'

import logging
logger = logging.getLogger('myapp')
hdlr = logging.FileHandler(outDir + '/pyText.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)


def textureFitsFile(inFilename, pF, wavelength, cake=True, doPlot=False):

    runNum = inFilename.split('_')[-1].split('.')[0]
    # print 'Reading ' + inFilename

    with open(inFilename, mode='rb') as fileobj:
        rawData = np.fromfile(fileobj, np.float32, num_X * num_Y).astype('float32')
    z = rawData.reshape((num_X, num_Y))

    dMax = 3.1
    dMin = 1.5

    x1, y1 = pF[0], pF[1]
    c = [ 2.79, 2.56, 2.45, 1.88, 1.60]
    cMult=[ 6, 2, 12, 12, 6]
    pkInd=['100', '002', '101', '102', '110']

    if doPlot:
        fig, ax = plt.subplots(1, 1)
        for mC in c:
            ax.axvline(mC + 0.05, color='k', ls='dotted')
            ax.axvline(mC - 0.05, color='k', ls='dotted')

    if cake:
        xd, yd = pF[0] - np.arange(2048), pF[1] - np.arange(2048)
        xv, yv = np.meshgrid(xd, yd)
        psi = np.degrees(np.arctan2(-np.sin(pF[3])*xv + np.cos(pF[2]) * yv, np.cos(pF[3]) * (np.cos(pF[2]) * xv + np.sin(pF[2]) * yv)))
        #arcWid = np.abs(rhoStart - rhoEnd) / rhoN / 2.
        arcWid = rhoStep
    #for angle in np.linspace(rhoStart, rhoEnd, rhoN):
    datBuff = np.zeros((np.arange(rhoStart, rhoEnd, rhoStep).size, 6))
    for angle, k in zip(np.arange(rhoStart, rhoEnd, rhoStep), count(0)):
        x2, y2 = x1 + arcRad * np.cos(np.radians(angle)), y1 + arcRad * np.sin(np.radians(angle))
        if cake:
            cond1 = (psi > angle) & (psi < angle + arcWid)
            cakeMat = np.where(cond1)
            x, y = xv[cakeMat]+pF[0], yv[cakeMat]+pF[1]
        else:
            x, y = np.linspace(x1, x2, num), np.linspace(y1, y2, num)
        dSpac = guessDSpacing(pF, x, y, wavelength)
        kKeys = np.where((dMin < dSpac) & (dSpac < dMax) & (0 < x) & (x < num_X) & (0 < y) & (y < num_Y))
        x = x[kKeys]
        y = y[kKeys]
        dSpac=dSpac[kKeys]
        if len(dSpac) is 0:
            print 'No measurement for angle of %d.' % angle
            pass
        else:
            if cake:
                binedges = np.linspace(np.min(dSpac),np.max(dSpac),256)
                dg = np.digitize(dSpac,binedges)
                ind = np.argsort(dg)
                splits = np.where(np.diff(dg[ind]))[0]
                sliceInt = np.split(z[cakeMat][kKeys][ind],splits)[1:]
                zi= np.array([np.maximum(x.mean(),1) for x in sliceInt])
                dSpac = np.array(binedges[:-1])
            else:
                zi = textResample(z,x,y)
            if np.abs(dSpac.size - zi.size) > 0:
                print 'Choose a wider rho angle!'
            if doPlot:
                ax.plot(dSpac,zi,'-')
            backG = np.mean(zi[np.where((2.0 - 0.05 < dSpac) & (2.0 + 0.05 > dSpac))])
            cbuff = np.zeros(5)
            for mC,i in zip(c,count(0)):
                thisPeak = np.where((mC - 0.05 < dSpac) & (mC + 0.05 > dSpac))
                thisArea = np.mean(zi[thisPeak]) / backG / cMult[i]
                cbuff[i] = thisArea
            datBuff[k,0] = angle + arcWid / 2.0
            datBuff[k,1:] = cbuff

    try:
        os.makedirs(outDir + 'dat/' + runNum)
    except:
        pass

    for thisInd,i in zip(pkInd,count(0)):
        datOut = outDir + 'dat/' + runNum + '/' + runNum + '.' + thisInd
        thisDat= np.zeros((datBuff[:,0].size,3))
        thisDat[:,0] = 90
        thisDat[:,1] = datBuff[:,0]
        thisDat[:,2] = datBuff[:,i+1]
        np.savetxt(datOut, thisDat, fmt='%.5f')

    if doPlot:
        ax.set_xlim(dMin, dMax)
        figOut = outDir + 'fig/' + runNum + '.png'
        plt.savefig(figOut)
        plt.close()

lam = 0.1441
ceriaDat = np.array([1.01038456e+03, 1.02344190e+03,  7.03712550e+03, -1.89784818e-02, 1.60107261e+00])
# fList = glob.glob('/home/chris/Tubes/ge/GE4Tube*sum')
fList = glob.glob('/home/chris/Tubes/ge/GE4Tube_0*sum')
a = (ceriaDat,lam)

runQueues(fList,textureFitsFile, a)

# for f, n in zip(fList,count()):
#     textureFitsFile(f, *a)
#     print '%d of %d files complete.' % (n+1, len(fList))