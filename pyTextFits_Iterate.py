#!/usr/bin/python
#import threading
from multiprocessing import Process, Queue
import os
import glob
import matplotlib
matplotlib.use('Agg') #Allows for saving figures without producing a window
import matplotlib.pyplot as plt
from itertools import count
import numpy as np
from scipy.optimize import curve_fit
from pyTexture import guessDSpacing, getPeakCens, gauss, getHKL, fakePseudoVoigt, textResample

np.seterr(divide='ignore')

if os.path.isfile('doPlots'):
    iPlot = True
else:
    iPlot = False

num_X = 2048
num_Y = 2048
num = 2000
areaScale = np.sqrt(np.log(2)/np.pi)
arcRad = 1024 * 1.42
rhoStart, rhoEnd = 0, 20
rhoN = 10
outDir = './'

import logging
logger = logging.getLogger('myapp')
hdlr = logging.FileHandler(outDir + '/pyText.log')
formatter = logging.Formatter('%(asctime)s %(levelname)s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)


def textureFitsFile(inFilename, pF, wavelength, cake=True, doPlot=False):
    # ptList = []
    xLi = []
    yLi = []
    areaLi = []
    angleList = []
    cLi = []
    pkLi = []
    runNum = inFilename.split('_')[-1].split('.')[0]
    print 'Reading ' + inFilename
    with open(inFilename, mode='rb') as fileobj:
        rawData = np.fromfile(fileobj, np.float32, num_X * num_Y).astype('float32')
    z = rawData.reshape((num_X, num_Y))

    #peakCens = getPeakCens('hcp',[3.23433391, 5.163978], 2)
    peakCens = np.array([2.78126, 2.4591, 2.585,1.90, 1.60544, 1.2384])
    #peakHKL  = getHKL('hcp',2)
    peakHKL = [(1,0,0),(1,0,1),(0,0,2),(1,0,2),(1,1,0),(2,0,2)]
    peakWids = np.repeat(0.01,len(peakCens))
    dMax = np.max(peakCens) + 0.5
    x1, y1 = pF[0], pF[1]
    if doPlot:
        fig, ax = plt.subplots(1, 1)
    if cake:
        xd, yd = pF[0] - np.arange(2048), pF[1] - np.arange(2048)
        xv, yv = np.meshgrid(xd, yd)
        psi = np.degrees(np.arctan2(-np.sin(pF[3])*xv + np.cos(pF[2]) * yv, np.cos(pF[3]) * (np.cos(pF[2]) * xv + np.sin(pF[2]) * yv)))
        arcWid = np.abs(rhoStart - rhoEnd) / rhoN / 2.
    for angle in np.linspace(rhoStart, rhoEnd, rhoN):
        x2, y2 = x1 + arcRad * np.cos(np.radians(angle)), y1 + arcRad * np.sin(np.radians(angle))
        if cake:
            cond1 = (psi > angle - arcWid) & (psi < angle + arcWid)
            cakeMat = np.where(cond1)
            x, y = xv[cakeMat]+pF[0], yv[cakeMat]+pF[1]
        else:
            x, y = np.linspace(x1, x2, num), np.linspace(y1, y2, num)
        dSpac = guessDSpacing(pF, x, y, wavelength)
        kKeys = np.where((dSpac < dMax) & (0 < x) & (x < num_X) & (0 < y) & (y < num_Y))
        x = x[kKeys]
        y = y[kKeys]
        dSpac=dSpac[kKeys]
        if len(dSpac) is 0:
            print 'No measurement for angle of %d.' % angle
            pass
        else:
            if cake:
                binedges = np.linspace(np.min(dSpac),np.max(dSpac),1024)
                dg = np.digitize(dSpac,binedges)
                ind = np.argsort(dg)
                splits = np.where(np.diff(dg[ind]))[0]
                sliceInt = np.split(z[cakeMat][kKeys][ind],splits)[1:]
                zi= np.array([x.mean() for x in sliceInt])
                dSpac = np.array(binedges[1:])
            else:
                zi = textResample(z,x,y)
            if doPlot:
                plt.cla()
                ax.plot(dSpac[:np.size(zi)],zi,'k.-')
            '''
            BGmm1=np.argmin(np.abs(dSpac-(2.17+0.01*7)))
            BGmm2=np.argmin(np.abs(dSpac-(2.17-0.01*7)))
            BGmmRange = range(np.max((np.min((BGmm1,BGmm2)),0)),np.min((np.max((BGmm1,BGmm2)),len(zi))))
            BGArea = np.sum(zi[BGmmRange])
            '''

            for c in [xx for xx in peakCens if xx > np.min(dSpac)]:
                # mmX = np.argmin(np.abs(dSpac-c))
                # print c
                mm1=np.argmin(np.abs(dSpac-(c+peakWids[peakCens==c][0]*10)))
                mm2=np.argmin(np.abs(dSpac-(c-peakWids[peakCens==c][0]*10)))
                mmRange = np.array(range(np.max((np.min((mm1,mm2)),0)),np.min((np.max((mm1,mm2)),len(zi)))))
                # if np.any(zi[mmRange] <= 0) or np.all(c > (dSpac[mmRange])):
                if False:
                    pass
                else:
                    '''
                    p0 = [np.max(zi[mmRange]) - np.min(zi[mmRange]),dSpac[mmX],peakWids[peakCens==c],np.min(zi[mmRange])]
                    p0 = [np.max(zi[mmRange]) - np.min(zi[mmRange]),dSpac[mm2+np.argmax(zi[mmRange])],peakWids[peakCens==c],np.min(zi[mmRange])]
                    try:
                        coeff, _ = curve_fit(gauss, dSpac[mmRange], zi[mmRange], p0,xtol=1e-2)
                    except:
                        coeff = np.array([0, c, 1e-3, 0])
                    '''
                    coeff    = np.array([np.sum(zi[mmRange])/np.min(zi[mmRange]),c,1e-3,0])
                    coeff[0] = np.abs(coeff[0])
                    coeff[2] = np.abs(coeff[2])

                    mRan = zi[mmRange].argsort()[:5]
                    bx = np.hstack(([mmRange[0],mmRange[-1]], mmRange[mRan]))
                    by = np.hstack(([zi[mmRange[0]],zi[mmRange[-1]]],zi[mmRange][mRan]))
                    bkfit = np.polyfit(bx,by,0)
                    newArea = np.sum(zi[mmRange]) - np.sum(np.polyval(bkfit,mmRange))

                    p0 = [np.max(zi[mmRange]) - np.min(zi[mmRange]),dSpac[mm2+np.argmax(zi[mmRange])],peakWids[peakCens==c],np.min(zi[mmRange])]
                    try:
                        coeff, pcov = curve_fit(fakePseudoVoigt, dSpac[mmRange], zi[mmRange]-np.polyval(bkfit,mmRange), p0,xtol=1e-2,method='lm')
                        perr = np.sqrt(np.diag(pcov))
                        print 'FIRST LINE'
                        print pcov
                        print ''
                        print coeff
                        print perr
                        print ''
                    except:
                        coeff = np.array([0, c, 1e-3, 0])
                                       
                    if doPlot:
                        resca = np.linspace(dSpac[mmRange[0]],dSpac[mmRange[-1]],500)
                        resca2 = np.linspace(mmRange[0],mmRange[-1],500)
                        # pk_fit = gauss(resca, *coeff)
                        # ax.plot(dSpac[mmRange],zi[mmRange],'.')
                        # ax.axvline(dSpac[mmRange[0]])
                        # ax.axvline(dSpac[mmRange[-1]])
                        ax.plot(dSpac[bx],by,'ro')
                        ax.plot(resca,fakePseudoVoigt(resca,*coeff)+np.polyval(bkfit,resca2),'-',lw=2)

                    if coeff[2] > 0.1:
                        print coeff
                        # ax.plot(resca,gauss(resca,*p0),'k--')
                        pass
                    else:
                        zi[mmRange] = zi[mmRange] - fakePseudoVoigt(dSpac[mmRange],*coeff)
                        # if doPlot: ax.plot(resca,pk_fit,'-')
                        # areaLi = np.append(areaLi, coeff[0] * coeff[2] * areaScale)
                        areaLi = np.append(areaLi, newArea)
                        if not cake:
                            xLi = np.append(xLi, np.interp(coeff[1],range(np.size(kKeys)),x))
                            yLi = np.append(yLi, np.interp(coeff[1],range(np.size(kKeys)),y))
                        angleList = np.append(angleList, angle)
                        cLi = np.append(cLi, coeff[1])
                        pkLi= np.append(pkLi,np.where(peakCens==c)[0][0])
                        if (c == 2.585):
                            pass
                        else:
                            np.place(peakCens,peakCens == c,coeff[1])
                        # np.place(peakWids,peakCens == c,coeff[2])
        if doPlot:
            ax.set_ylim(1,1e3)
            ax.set_xlim(1.2,3)
            ax.set_yscale('log')
            # print '%.2f' % angle
            ax.plot(dSpac[:np.size(zi)],zi,'r.')
            #plt.show()
            figOut = outDir + 'fig/' + runNum + '_' + str(np.round(angle)) + '.png'
            print figOut
            plt.savefig(figOut)
    theta = 90 - 0 * cLi  #- np.degrees(np.arcsin(wavelength / (2 * cLi))/2)
    try:
        os.makedirs(outDir + runNum)
    except:
        pass
    outputName = inFilename.split('/')[-1].split('.')[0]
    if cake:
        outputName += 'cake'
    for K, n in zip(peakHKL,range(len(peakHKL))):
        nKeys = np.where(pkLi == n)
        if len(nKeys[0]) > 5:
            #print '**\n%d%d%d has ' % K + str(len(nKeys[0])) + ' fits.'
            outFile = outDir + runNum + '/' + outputName + '.%d%d%d' % K
            #print 'Output to ' + outFile
            M = np.vstack((theta[nKeys], angleList[nKeys], areaLi[nKeys]))
            np.savetxt(outFile,M.T,fmt='%.6f',delimiter=' ')
    print 'Files output to ' + outFile

lam = 0.1441
ceriaDat = np.array([  1.01038456e+03, 1.02344190e+03,  7.03712550e+03, -1.89784818e-02, 1.60107261e+00])
fList = glob.glob('../tubes/GE4Tube*sum')
a = (ceriaDat,lam)

for f, n in zip(fList,count()):
    oName = outDir + f.split('_')[-1].split('.')[0] + '/'

#    if (~os.path.exists(oName)):
#     try:
    textureFitsFile(f, *a, doPlot=iPlot)
    # except:
    #     st = 'File error: ' + f.split('_')[-1].split('.')[0]
    #     logger.error(st)
    #     print st
#    else:
#        pass
    print '%d of %d files complete.' % (n+1, len(fList))
