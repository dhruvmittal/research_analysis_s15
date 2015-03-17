import averageData
import os
import math
import scipy
from scipy import integrate
import numpy as np

import aggregateData as ag
import matplotlib.pyplot as plt

BETA = 1.00
Nx = 40
L = 1.00 * Nx
HBAR = 1
OMEGA_FREQ = 1

# var_num: [0:Npart, 1:EK, 2:EP, 3:EVext, 4:ET]
def prepareForPlot(path, var_num, need_coupling=False):
    all_data = ag.getDataFromSubDirs(ag.listSubDir(path))
    x = []
    y = []
    for datum in all_data:
        x.append(datum[0])
        y.append(datum[1][var_num])
    if need_coupling:
        return x, y, path.split('/')[-1].split('_')[-1]
    return x,y

def calculateOmegaNaught(beta, beta_mu):
    result = 0.0
    prev = 1
    for n in range(0, 100):
        this =  math.log(1 + math.exp(-1 * beta * (HBAR * OMEGA_FREQ * (float(n) + 1./2.)) - beta_mu/beta))
        result = result + this
        try:
            ratio = this / prev 
        except:
            break
        prev = this
    return result / (-1 * beta)

def OmegaNaughtDistribution(beta, beta_mu_list=None, beta_mu_list_path=None):
    if beta_mu_list == None:
        f = open(beta_mu_list_path)
        beta_mu_list = []
        for line in f:
            beta_mu_list.append(float(line))
    omega_naught_list = []
    for beta_mu in beta_mu_list:
        omega_naught_list.append(calculateOmegaNaught(beta, beta_mu))
    return omega_naught_list

def NumberDensityNoCoupling(omega_naught_list):
    # dOmega_0 / d mu
    return np.gradient(omega_naught_list)



def prepareBMvsDimensionlessDensity(path_no_coupling, path):
    no_coupling = prepareForPlot(path_no_coupling, 0)
    yes_coupling = prepareForPlot(path, 0)
    x = []
    y = []
    for i in range(len(no_coupling[0])):
        try:
            x.append(math.log(no_coupling[0][i]))
            y.append(yes_coupling[1][i]/no_coupling[1][i])
        except:
            pass
    return x,y

def prepareBMvsDimensionlessContact(path):
    #no_coupling = prepareForPlot(path_no_coupling, 2, need_coupling=True)
    yes_coupling = prepareForPlot(path, 2, need_coupling=True)
    x = []
    y = []
    for i in range(len(yes_coupling[0])):
            x.append(math.log(yes_coupling[0][i]))
            g = float(yes_coupling[2])
            interaction_avg_energy =  yes_coupling[1][i]
            contact = -1 * g * interaction_avg_energy
            gamma_squared = BETA * g**2
            y.append(contact * math.pi * BETA**2 / (2 * L * gamma_squared))
        #except:
            #pass
    return x,y

def prepareBMvsDimensionlessPressure(path):
    #no_coupling = ([],[]) 
    #yes_coupling = ([],[])
    x = []
    y = []
    #with open(path_no_coupling) as f:
        #for line in f:
            #l = line.split('\t')
            #no_coupling[0].append(l[0])
            #no_coupling[1].append(l[1])
    with open(path) as f:
        for line in f:
            l = line.split('\t')
            x.append(float(l[0]))
            y.append(float(l[1]))
        f.close()
    x = x[1:-1]
    new_y = averageData.movingAverage(y)
    integrated_y = [1.0] * len(y)
    for i in range(1, len(y)):
        integrated_y[i] = integrate.simps(new_y[0:i], x[0:i])
    return x, integrated_y

def plot(data_tuple):
    plt.plot(data_tuple[0], data_tuple[1], marker='.') #, linestyle='')

def plotAll(initial_path, var_num):
    for dir in os.listdir(initial_path):
        plot(prepareBMvsDimensionless(
            os.path.join(initial_path, 'Bare_coupling_0.01/'),
            os.path.join(initial_path, dir), 
            var_num)
            )

def writeAllDataToFile(initial_path, output_path, case):
    for dir in os.listdir(initial_path):
        f = open(os.path.join(output_path, dir), 'w')
        if case == 0: #Density
            data_tuple = prepareBMvsDimensionlessDensity(
                    os.path.join(initial_path, 'Bare_coupling_0.01/'),
                    os.path.join(initial_path, dir), 
                    )
        if case == 2: #Contact
            data_tuple = prepareBMvsDimensionlessContact(
                    os.path.join(initial_path, dir),
                    )
        if case == 3: #Pressure
            data_tuple = prepareBMvsDimensionlessPressure(
                    os.path.join(initial_path, dir)
                    )
        for i in range(len(data_tuple[0])):
            try:
                f.write(str(data_tuple[0][i]) + '\t' + str(data_tuple[1][i]) + '\n')
            except:
                pass


def plotFromFile(filename, mv_avg=False):
    x = []
    y = []
    with open(filename) as f:
        for line in f:
            l = line.split('\t')
            x.append(float(l[0]))
            y.append(float(l[1]))
    f.close()
    if mv_avg:
        new_y = averageData.movingAverage(y)
        error = [a - b for a, b in zip(new_y, y[1:-1])]
        plot((x[1:-1], new_y))
        plt.errorbar(x[1:-1], new_y, yerr=error, fmt='')
    else:
        plot((x,y))


def plotAllFromFile(initial_path, mv_avg=False):
    for f in os.listdir(initial_path):
        if f[0] != '.':
            plotFromFile(os.path.join(initial_path, f), mv_avg=True)

        
def getDensityBMPlot(calculate=False):
    if calculate == True:
        writeAllDataToFile('../exec_beta_2.0/beta_2.0/', '../data_Density/',
                case=0)
    plotAllFromFile('../data_Density', mv_avg=True)
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$N/N_0$')
    plt.show()


def getContactBMPlot(calculate=False) :
    if calculate == True:
        writeAllDataToFile('../exec_beta_2.0/beta_2.0/', '../data_Contact/',
                case=2)
    plotAllFromFile('../data_Contact/')
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$C\pi\beta^2/(2L\lambda^2)$')
    plt.show()


def getPressure(calculate=False):
    if calculate == True:
        writeAllDataToFile('../data_Density/', '../data_Pressure/', 3)
    plotAllFromFile('../data_Pressure/')
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$P/P_0$')
    plt.show()

    


#getDensityBMPlot()
#getContactBMPlot()
getPressure(calculate=True)
#print calculateOmegaNaught(BETA, 0)
#print NumberDensityNoCoupling(OmegaNaughtDistribution(BETA, None, '../fugacity_input.txt'))
