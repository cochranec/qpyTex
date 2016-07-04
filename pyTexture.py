'''
Created on Nov 21, 2014

@author: chris
'''

from itertools import combinations_with_replacement, combinations, product, repeat,chain
# from warnings import warn
from numpy import sin, cos, sqrt, arctan2, exp, log, divide, clip
from scipy.optimize import leastsq, curve_fit
# import numpy as np
import threading
from Queue import Queue
from time import strftime
# import matplotlib.pyplot as plt

def isBCC(x):
    if sum(x) == 0:
        return 0
    else:
        return not (sum(x) % 2) == 1

def isFCC(x):
    if sum(x) == 0:
        return 0
    else:
        return not ((x[0] + x[1]) % 2) and not ((x[1]+x[2]) % 2)

def isHCP(x):
    if sum(x) == 0:
        return 0
    else:
        return not (x[2] % 2) or (x[0] + 2*x[1] % 3)

def getHKL(sys, n=3):
    if sys.lower() == "bcc":
        fxn = isBCC
    elif sys.lower() =="fcc":
        fxn = isFCC
    elif sys.lower() =="hcp":
        fxn = isHCP
        y = [[x[0][1], x[0][0], x[1]] for x in product(combinations_with_replacement(range(n+1),2),range(n+1))]
        return [tuple(x) for x in y if fxn(x)]
    else:
        warn("No crystal system given.", SyntaxWarning)
        return []
    return [x for x in combinations_with_replacement(range(n+1), 3) if fxn(x)]

def getPeakCens(sys, a, n=3):
    hkl = getHKL(sys, n)
    squares = []
    if 'hcp' in sys.lower():
        for x in hkl:
            squares.append(4/3. * (x[0]**2 + x[0] * x[1] + x[1]**2)/a[0]**2 + x[2]**2/a[1]**2)
    else:
        for x in hkl:
            squares.append((x[0]**2 + x[1]**2 + x[2]**2)/a)
    return [1 / sqrt(x) for x in squares]

def guessDSpacing(p, xv, yv, wavelength):
    [xcg,ycg,betag,phig,dg] = p
    pxC = 200e-6
    xd = (xv - xcg) * pxC
    yd = (yv - ycg) * pxC
    TA = cos(phig)**2 *  ( cos(betag) * xd + sin(betag) * yd)**2
    TB =                 (-sin(betag) * xd + cos(betag) * yd)**2
    TC = (dg + sin(phig)*( cos(betag) * xd + sin(betag) * yd))**2
    R = sqrt(TA + TB)
    D = sqrt(TC)
    delR = R * (1 - D / sqrt(TA + TB + TC))
    delD = D * (1 - D / sqrt(TA + TB + TC))
    trueR = R - delR
    trueD = D - delD
    TT = arctan2(trueR, trueD)
    return wavelength / TT

def gauss(x, *p):
    A, mu, sigma, y0 = p
    return A * exp(-(x - mu) ** 2 / (2.*sigma ** 2)) + y0

def gaussDoublet(x, *p):
    A1, mu1, sigma1, A2, mu2, sigma2, y0 = p
    return gauss(x, A1, mu1, sigma1, y0) + gauss(x, A2, mu2, sigma2, 0)

def fakePseudoVoigt(x, *p):
    height, center, hwhm, y0 = p
    shape=0.4
    return height*((1-shape)*exp(-log(2)*(divide((x-center),hwhm))**2)+shape*1/(1+(divide(x-center,hwhm))**2)) +  y0

def PseudoVoigt(x, *p):
    height, center, hwhm, sha, y0 = p
    shape = sin(sha)**2
    return height*((1-shape)*exp(-log(2)*(divide((x-center),hwhm))**2)+shape*1/(1+(divide(x-center,hwhm))**2)) +  y0
    #return A * exp(-(x - mu) ** 2 / (2.*sigma ** 2)) + y0

def textResample(z, x, y):
    x1 = x.astype('int')
    x2 = clip(x1+1,0,2047)
    y1 = y.astype('int')
    y2 = clip(y1+1,0,2047)

    EPSILON = (x1 + .5 - x) * -2
    ETA     = (y1 + .5 - y) * -2
    N1 = .25 * (1 - EPSILON) * (1 - ETA)
    N2 = .25 * (1 + EPSILON) * (1 - ETA)
    N3 = .25 * (1 + EPSILON) * (1 + ETA)
    N4 = .25 * (1 - EPSILON) * (1 + ETA)

    zo = z[x1,y1] * N1 + z[x2,y1] * N2 + z[x2,y2] * N3 + z[x1,y2] * N4
    return zo


def runQueues(files, fxn, params):
    nThread = 4
    nFiles  = len(files)
    # startTime = time.time()

    q = Queue(maxsize=0)
    for fs in files:
        q.put(fs)

    print strftime('%H:%M:%S')  + ' Starting queue with %d files.' % (len(files))

    for _ in range(nThread):
        worker = threading.Thread(target=feedQueue, args=(q, fxn, params,))
        worker.setDaemon(True)
        worker.start()

    q.join()

    # print 'Done.  Average time per file: %.3fs' % ((time.time() - startTime)/nFiles)


def feedQueue(q, fxn, param):
    while True:
        f = q.get()
        sizeOfQueue = q.qsize()
        fxn(f, *param)
        print strftime('%H:%M:%S') +' ' + str(sizeOfQueue) + ' files remain in queue.'
        # nFilesComplete = nFile - sizeOfQueue
        # timeSpent = time.time() - startT
        # timeRemaining = numpy.round(timeSpent / nFilesComplete * sizeOfQueue)
        # m, s = divmod(timeRemaining, 60)
        # h, m = divmod(m, 60)
        # print 'File ', nFilesComplete, '/', nFile, '-',
        # print "%d:%02d:%02d to completion." % (h, m, s)
        q.task_done()

#def guessDoublet(foo):

# if __name__ == '__main__':
#     bb=False
#     if bb:
#         p = [3.23, 5.15]
#         a = np.array(getPeakCens('hcp',p,2))
#         err = lambda p: [np.array(getPeakCens('hcp',p,2))[0] - 2.581989, np.array(getPeakCens('hcp',p,2))[2] - 2.462139]
#         print err(p)
#         x, cov = leastsq(err, p, args=())
#         print x
#     if ~bb:
#         xr = np.arange(-10,10,.1)
#         pCor= [1, 0, .3,
#                2, 1, .2,
#                0.1]
#         foo = gaussDoublet(xr, *pCor) + np.random.rand(len(xr))*.2
#
#         errfunc = lambda p, x, y: gauss(x, *p) - y
#
#         c1 = np.argmax(foo)
#         guess1=[1, xr[c1],.1,0]
#         f1 = gauss(xr, *guess1)
#
#         p1, _ = curve_fit(gauss, xr, foo, guess1,xtol=1e-8)
#
#         foo2 = -errfunc(p1,xr,foo)
#         c2 = np.argmax(foo2)
#         guess2=[1, xr[c2],.1,0]
#         f2 = gauss(xr, *guess2)
#         p2, _ = curve_fit(gauss, xr, foo2, guess2,xtol=1e-8)
#
#         f3 = gauss(xr, *p1) + gauss(xr, *p2)
#         fig, ax = plt.subplots(1,1)
#
#         ax.plot(xr, foo,'.')
#         ax.plot(xr, foo2,'.')
#         ax.plot(xr, f1+f2-guess1[-1],'g-', label='Initial Guess')
#         ax.plot(xr, f3,'r-', label='Fit')
#         ax.legend()
#         peakCen = [guess1[1], guess2[1]]
#         print peakCen
#         plt.show()