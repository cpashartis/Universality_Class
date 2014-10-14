# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Wed May 14 15:32:02 2014

@author: christoforos
"""

import numpy as np
from subprocess import call

def analyze(combos,dat_file = '01.data'):

    data=np.loadtxt(dat_file)
    
    #axis  = 0, i.e. over columns
    var = data.var(axis = 0)
    
    #propagation of error for a log is sqrt((1/x)**2*var)
    err = ((1/data)**2*var)**(0.5)
    #storable values
    store = np.zeros((len(var),4))
    store[:,:2]=combos
    store[:,2] = -1*np.log(np.mean(data,axis = 0)) #means of entropies
    store[:,3] = (var/len(var))**(0.5) #uncertainty of means
    #reshape

    np.savetxt('R_Entropy.data', store, fmt='%.6f', header = 'x\ty\tS2\tUnc')

#list of stuff to plot
#assuming the data is all being done in the same run
    
def plot_m(dim=0, dat_file = '01.data'):
    
    from matplotlib import use 
    use('agg')

    import matplotlib.pyplot as plt
    
    data = np.loadtxt(dat_file)
    
    #axis  = 0, i.e. over columns
    var = data.var(axis = 0)
    #propagation of error for a log is sqrt((1/x)**2*var)
    err = ((1/data)**2*var)**(0.5)
    
        #plot over m, i.e. pathlength
    range(len(data))
    -1*np.log(data[:,dim])
    plt.errorbar(range(len(data)), -1*np.log(data[:,dim]), yerr = err[dim])
    
    

def prep_A():
    "Makes the regionA data file."
    #load in param.dat
    param_file = "param.dat"
    regionA = "regionA.dat"
    params = np.loadtxt(param_file)
    
    #temp_write = [(params[0]-1) * (params[1]-1) ]
    temp_write = []
    A_file = open(regionA, 'w')
    A_file.write('%i' %( (params[0]-1) * (params[1]-1) ) ) #writes number of lattice combos first
    
    hold = []
    for y in range(1,int(params[1])): #loop over n_y
        for x in range(1,int(params[0])): #loop over n_x
            sites = np.zeros( ( params[1],params[0] ) ) #make 0 matrix of sites
            sites[:y,:x] = 1
            sites = sites.astype('int')
            hold.append([x,y])
            for line in sites:
                A_file.write('\n' + ' '.join([str(i) for i in line]))
            A_file.write('\n-99')
    
    A_file.close()    
    
    return hold    
        
if __name__ == '__main__':
    list_combos = prep_A()
    #run the simulation
    exe = 'tfim.out'
    call('./' + exe)
    
    #put data in file
    analyze(list_combos)