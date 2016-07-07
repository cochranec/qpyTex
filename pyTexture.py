'''
Created on Nov 21, 2014

@author: chris
'''

from itertools import combinations_with_replacement, combinations, product, repeat,chain
# from warnings import warn
from numpy import sin, cos, sqrt, arctan2, exp, log, divide, clip
# from scipy.optimize import leastsq, curve_fit
# import numpy as np
import threading
from Queue import Queue
from time import strftime, sleep, time
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

    q = Queue(maxsize=0)
    for fs in files:
        q.put(fs)

    monitor = threading.Thread(target=monitorQueue, args=(q, 30))
    monitor.setDaemon(True)
    monitor.start()

    for _ in range(nThread):
        worker = threading.Thread(target=feedQueue, args=(q, fxn, params,))
        worker.setDaemon(True)
        worker.start()

    q.join()


def feedQueue(q, fxn, param):
    while True:
        f = q.get()
        fxn(f, *param)
        q.task_done()


def monitorQueue(q, delay):
    start_size = q.qsize()
    start_time = time()

    print strftime('%H:%M:%S')  + ' Starting queue with %d files.' % (start_size)
    while True:
        sleep(delay)
        current_size = q.qsize()
        current_time = time()
        time_remain = ((current_time - start_time) / (start_size - current_size) * current_size) / 60
        print strftime('%H:%M:%S') + ' ' + str(current_size) + ' files remain in queue.  Process time remaining: ' + str(time_remain)
