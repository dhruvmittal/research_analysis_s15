#!/usr/bin/python

import numpy

MAX_STEPS = 22000

START_AT_NX = {40:90, 60:110}

# Format: [Npart, EK, EP, EVext, ET]

def get(filename):
    with open(filename) as f:
        for _ in xrange(START_AT_NX[60]):
            next(f)
        count_entries = 0
        Npart = 0
        EK = 0
        EP = 0
        EVext = 0
        ET =  0
        for line in f:
            l =  line.split()
            if int(l[0]) % 10 == 0:
                Npart = Npart + float(l[6])
                EK = EK + float(l[7])
                EP = EP + float(l[8])
                EVext = EVext + float(l[9])
                ET = ET + float(l[10])
                count_entries += 1
            else:
                pass
            if int(l[0]) == MAX_STEPS:
                break
        f.close()
    avg_Npart = Npart / count_entries
    avg_EK = EK / count_entries
    avg_EP = EP / count_entries
    avg_EVext = EVext / count_entries
    avg_ET = ET / count_entries
    return [avg_Npart, avg_EK, avg_EP, avg_EVext, avg_ET]

def getError(filename):
    with open(filename) as f:
        for _ in xrange(80):
            next(f)
        Npart  = []
        EK = [] 
        EP = []
        EVext = []
        for line in f:
            l = line.split()
            if int(l[0]) % 10 == 0:
                Npart.append(float(l[6]))
                EK.append(float(l[7]))
                EP.append(float(l[8]))
                EVext.append(float(l[9]))
                ET.append(float(l[10]))
            else:
                pass
            if int(l[0]) == MAX_STEPS:
                break
        f.close()
    stdev_Npart = numpy.std(Npart)
    stdev_EK = numpy.std(EK)
    stdev_EP = numpy.std(EP)
    stdev_EVext = numpy.std(EVext)
    stdev_ET = numpy.std(ET)
    return [stdev_Npart, stdev_EK, stdev_EP, stdev_EVext, stdev_ET]




def movingAverage(a, n=3):
    ret = numpy.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

