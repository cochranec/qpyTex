#!/usr/bin/python

import os
import glob
from itertools import count

import matplotlib
matplotlib.use('Agg') #Allows for saving figures without producing a window
import matplotlib.pyplot as plt

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
rhoStep = 5
rhoN = 90
outDir = './'

lam = 0.1441
ceriaDat = np.array([1.01038456e+03, 1.02344190e+03,  7.03712550e+03, -1.89784818e-02, 1.60107261e+00])
# fList = glob.glob('/home/chris/Tubes/ge/GE4Tube*sum')
fList = glob.glob('/home/chris/Tubes/ge/CC*sum')

import logging
logger = logging.getLogger('myapp')
hdlr = logging.FileHandler(outDir + '/pyText.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)


def textureFitsFile(inFilename, pF, wavelength, cake=True, doPlot=False):
    dMax = 3.1
    dMin = 1.5

    x1, y1 = pF[0], pF[1]
    c = [2.79, 2.56, 2.45, 1.88, 1.60]
    pkInd = ['100', '002', '101', '102', '110']
    cMult = [6, 2, 12, 12, 6]

    runNum = inFilename.split('_')[-1].split('.')[0]
    with open(inFilename, mode='rb') as fileobj:
        rawData = np.fromfile(fileobj, np.float32, num_X * num_Y).astype('float32')
    z = rawData.reshape((num_X, num_Y))

    if doPlot:
        fig, ax = plt.subplots(1, 1)
        for mC in c:
            ax.axvline(mC + 0.05, color='k', ls='dotted')
            ax.axvline(mC - 0.05, color='k', ls='dotted')

    if cake:
        xd, yd = pF[0] - np.arange(2048), pF[1] - np.arange(2048)
        xv, yv = np.meshgrid(xd, yd)
        psi = np.degrees(
            np.arctan2(
                -np.sin(pF[3])*xv + np.cos(pF[2]) * yv, np.cos(pF[3]) * (np.cos(pF[2]) * xv + np.sin(pF[2]) * yv)
            )
        )
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
        zSlice = z[cakeMat][kKeys]

        if len(dSpac) is 0:
            print 'No measurement for angle of %d.' % angle
            pass
        else:
            if cake:
                sliceStyle = True
                if sliceStyle:
                    # Split the d-spacing space into 256 equally spaced bins and take the mean.
                    # The number of pixels per bin will vary.
                    binedges = np.linspace(np.min(dSpac),np.max(dSpac),256)
                    dg = np.digitize(dSpac,binedges)
                    zi = np.zeros(binedges.size)
                    for ii in np.intersect1d(dg,np.arange(256)):
                        ind = np.where(dg == ii)
                        zi[ii] = np.sum(zSlice[ind])
                    # Fill empty slices with the average of adjacent slices
                    for ii in np.setxor1d(dg,np.arange(256)):
                        if ii < np.size(zi):
                            zi[ii] = np.mean([zi[ii-1],zi[ii+1]])
                    dSpac = np.array(binedges)
                else:
                    # Split d-Spacing array into 256 equal pixel number segments and take the mean
                    # The number of pixels in each bin is equal.
                    # NB: Throws out last N - 256 * k data points
                    ind = np.argsort(dSpac)
                    ll =  dSpac.size / 256
                    zi = np.mean(np.reshape(zSlice[ind][:ll*256], (256,-1),), axis=1)
                    dSpac = np.mean(np.reshape(dSpac[ind][:ll*256], (256,-1), ), axis=1)
            else:
                zi = textResample(z,x,y)
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

a = (ceriaDat,lam)
runQueues(fList,textureFitsFile, a)
