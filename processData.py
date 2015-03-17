import averageData
import os
import math
import scipy
from scipy import integrate
import numpy as np

import aggregateData as ag
import matplotlib.pyplot as plt

#--------- CONSTANTS ----------------#
INPUT = '../exec/beta_8.0'
BETA = 8.0
Nx = 40

L = 1.00 * Nx

HBAR = 1
OMEGA_FREQ = 1
dBetaMu = .5

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
    #prev = 1
    mu  = beta_mu / beta
    for n in range(0, 100):
        energy_n = HBAR * OMEGA_FREQ * (n + .5)
        current = math.log(1 + math.exp(-1 * beta * (energy_n - mu)))
        result = result + current
    return 2 * result 

def NumberDensityNoCoupling(omega_naught_list):
    # dOmega_0 / d mu
    return np.gradient(omega_naught_list) / dBetaMu

def OmegaNaughtDistribution(beta, beta_mu_list=None, fugacity_list_path=None):
    if beta_mu_list == None:
        f = open(fugacity_list_path)
        beta_mu_list = []
        for line in f:
            # convert fugacity to betaMu 
            beta_mu_list.append(math.log(float(line)))
    omega_naught_list = []
    for beta_mu in beta_mu_list:
        omega_naught_list.append(calculateOmegaNaught(beta, beta_mu))
    return omega_naught_list

def calc_log_z(beta, z):
	result = 0.0
	for n in range(0, 100):
		energy_n = HBAR * OMEGA_FREQ * (n + .5)
		current = math.log(1 + z * math.exp(-1 * beta * energy_n))
		result = result + current
	return 2 * result

def number_p_distribution_numerical(beta, fugacity_list_path):
	f = open(fugacity_list_path)
	fugacity_list = []
	for line in f:
		fugacity_list.append(float(line))
	f.close()
	log_z_list = []
	for z in fugacity_list:
		log_z_list.append(calc_log_z(beta, z))
	#get Number of particles now
	dy = np.gradient(log_z_list)
	dx = np.gradient(fugacity_list)
	number_p = []
	for i in range(len(dy)): 
		number_p.append(fugacity_list[i] * (dy[i]/dx[i]))
	return number_p
			
def number_p_distribution(beta, fugacity_list_path):
	f = open(fugacity_list_path)
	fugacity_list = []
	for line in f:
		fugacity_list.append(float(line))
	f.close()
	number_p_dist = []
	for z in fugacity_list:
		result = 0.0
		for n in range(0,100):
			energy_n = HBAR * OMEGA_FREQ * (n + .5)
			current = (2 * z * math.exp(-1 * beta * energy_n))/(1 + z * math.exp(-1*beta*energy_n))
			result = result + current
		number_p_dist.append(result)
	return number_p_dist
	
def read_in_number_p_distribution(beta):
	f = open("number_particles_analytical/beta_" + "%0.1f" % beta) 
	val = f.readline().split(',')
	for i in range(len(val) -1):
		val[i] = float(val[i])
	return val

def prepareBMvsDimensionlessDensity_2(path):
    #no_coupling_results = NumberDensityNoCoupling(OmegaNaughtDistribution(BETA, None, '../fugacity_input.txt'))
    #no_coupling_results = NumberDensityNoCoupling(number_p_distribution_numerical(BETA, '../fugacity_input.txt'))
    #no_coupling_results = NumberDensityNoCoupling(number_p_distribution(BETA, '../fugacity_input.txt'))
    no_coupling_results = read_in_number_p_distribution(BETA)
    yes_coupling = prepareForPlot(path, 0)
    x = []
    y = []
    for i in range(len(yes_coupling[0])):
        try:
            #print yes_coupling[0][i], yes_coupling[1][i], no_coupling_results[i]
            x.append(math.log(yes_coupling[0][i]))
            y.append(yes_coupling[1][i] / no_coupling_results[i])
            print yes_coupling[1][i] / no_coupling_results[i]
            #y.append(no_coupling_results[i])
            #y.append(yes_coupling[1][i])
        except:
	    print 'passed'
            pass
    return x,y

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
    no_coupling_omega = OmegaNaughtDistribution(BETA, fugacity_list_path='../fugacity_input_25.txt')
    yes_coupling_density = prepareForPlot(path, 0)
    #print yes_coupling_density
    x = []
    N = []
    for i in range(len(yes_coupling_density[0])):
        x.append(math.log(yes_coupling_density[0][i]))
        N.append(yes_coupling_density[1][i])
    new_N = averageData.movingAverage(N)
    integrated_N = [1.0] * len(new_N)
    x = x[1:-1]
    #print len(no_coupling_omega)
    #print len(new_N)
    #print len(x)
    no_coupling_omega = no_coupling_omega[1:-1]
    count = 0
    for i in range(1, len(new_N)):
        #print integrated_N
        integrated_N[i] = integrate.simps(new_N[0:i], x[0:i]) / no_coupling_omega[i]
        count += 1
        #print 'Number of traversals' + str(count)
    return x, integrated_N

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
            #data_tuple = prepareBMvsDimensionlessDensity(
                    #os.path.join(initial_path, 'Bare_coupling_0.01/'),
                    #os.path.join(initial_path, dir), 
                    #)
            data_tuple = prepareBMvsDimensionlessDensity_2(
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
        writeAllDataToFile(INPUT, '../data_Density/',
                case=0)
    plotAllFromFile('../data_Density', mv_avg=True)
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$N/N_0$')
    plt.show()


def getContactBMPlot(calculate=False) :
    if calculate == True:
        writeAllDataToFile(INPUT, '../data_Contact/',
                case=2)
    plotAllFromFile('../data_Contact/')
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$C\pi\beta^2/(2L\lambda^2)$')
    plt.show()


def getPressure(calculate=False):
    if calculate == True:
        writeAllDataToFile(INPUT, '../data_Pressure/', 3)
    plotAllFromFile('../data_Pressure/')
    plt.xlabel(r'$\beta\mu$')
    plt.ylabel(r'$P/P_0$')
    plt.show()

    


#getDensityBMPlot(True)
#getContactBMPlot()
getPressure(calculate=True)
#print calculateOmegaNaught(BETA, 0)
#print NumberDensityNoCoupling(OmegaNaughtDistribution(BETA, None, '../fugacity_input.txt'))


