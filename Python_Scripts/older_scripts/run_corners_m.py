# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Wed May 14 15:32:02 2014

@author: christoforos
"""
from matplotlib import use 
use('agg')

from matplotlib import rcParams, rc
rc('font',**{'family':'serif','serif':['Computer Modern'], 'size':14})
rc('text', usetex=True)
rc('xtick.major',pad=8)
rc('ytick.major',pad=8)
rcParams['axes.color_cycle'] = [(.4,0,0),(0,.4,0),(0,0,.4),(0,0,0),(.4,.4,0),(0,.4,.4),
(.4,0,.4),(.6,0,0),(0,.6,0),(0,0,.6),(0.6,.6,0),(0.6,0,.6),(0,.6,.6)]

import matplotlib.pyplot as plt
import numpy as np
from subprocess import Popen
import os
from shutil import copy
import glob

def analyze(combos, path = ''):

    data=np.loadtxt(path + '01.data')
    
    #axis  = 0, i.e. over columns
    estimated_var = data.var(axis = 0)/(len(data)-1) #sample variance divided by N-1 for ubiased estimate of the average of x
    
    #propagation of error for a log is sqrt((1/x)**2*var)
    err = ((1/np.mean(data,axis = 0))**2*estimated_var)**(0.5)
    #storable values
    store = np.zeros((len(estimated_var),4))
    store[:,:2]=combos #store x,y
    store[:,2] = -1*np.log(np.mean(data,axis = 0)) #mean entropy per configuration
    store[:,3] = err #uncertainty of means with propagation of error

    np.savetxt(path + 'R_Entropy.data', store, fmt='%.10f', header = 'x\ty\tS2_avg\tErr')

#list of stuff to plot
#assuming the data is all being done in the same run
    
def plot_n(path = '', **kwarg):
    
    data = np.loadtxt(path + '01.data')
    
    x = np.arange(len(data))
    f = plt.figure()
    ax = f.add_subplot(111)
    
    print "Plotting bin plots %s ... " %(path)
    
    for dim in range(len(data[0,:])):
        y =-1*np.log(data[:,dim])
        #plot over m, i.e. pathlength
        plt.plot(x, y,'k.')
        
        #linear fit the plot
        a,b = np.polyfit(x, y, 1)
        
        plt.plot(x, a*x + b, 'r')
        plt.text(.5,.025,"$y= %.4e x + %.4e$" %(a,b), fontsize = 12, transform=ax.transAxes)
        plt.xlabel("Bin Number")
        plt.ylabel("$S_2$")
        
        plt.savefig(path + 'plot_n%i.png' %(dim+1), bbox_inches = 'tight')
        plt.clf()
        
    print "Done"
    plt.close()

def prep_A(path,new_A = True):
    "Makes the regionA data file."
    #load in param.dat
    param_file = path+ "/param.dat"
    regionA = path + "/regionA.dat"
    params = np.loadtxt(param_file)
    
    hold = []
    if new_A == True: #only write a new regionA.dat if true
        #temp_write = [(params[0]-1) * (params[1]-1) ]
        temp_write = []
        A_file = open(regionA, 'w')
        A_file.write('%i' %( (params[0]-1) * (params[1]-1) ) ) #writes number of lattice combos first
        
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
    
    else:
        for y in range(1,int(params[1])): #loop over n_y
            for x in range(1,int(params[0])): #loop over n_x
                hold.append([x,y])
    
    return hold
    
def plot_m(path = ''):
    
    ann_loc = '../Ann_Check/'
    
    if isinstance(path, list) == False: #if not list
        path = list(path)
    
    for loc in path: #loop over all paths given, (generic handling, yay!)
        
        print 'Plotting pathlength for %s ...' %(loc)
        
        cluster = loc.rstrip('/') #5x4 etc
        entropies = []
        err = []
        m = []
        
        for m_name in glob.glob(loc+'*/R_Entropy.data'): #load in all data from entropy by searching for the files
            #append data to lists for now            
            m.append(m_name.split('/')[1].lstrip('m_'))  
            data = np.loadtxt(m_name)            
            entropies.append(data[:,2])
            err.append(data[:,3])
        
        #reformating data so it is easy to plot, i.e. connect all the m_* entropies
        #in the same cluster to be in the same row, different rows correspond to different clusters
        xy = data[:,:2]
        m = np.array(m, dtype= 'int') #sort m
        index = m.argsort() #sort m and EVERYTHING ELSE 
        m=m[index]
        #sort entropies and err according to order m_100..m_200..m_300
        entropies = np.array(entropies)[index].T #reformat
        err = np.array(err)[index].T #reformat

        
        #plot all ms with combinations on the same graph
        f = plt.figure()
        ax = f.add_subplot(111)
        f2 = plt.figure()
        ax2 = f2.add_subplot(111)
        corner_cut_all = []
        
        Ann_Entropy = np.loadtxt(ann_loc+cluster+'.dat') #load ann's data, note my combinations
        #work in the same order as ann's i corresponds to the row
        
        for i in range(len(data)):
            
            corner_cut = '%ix%i' %(xy[i,0], xy[i,1]) #.join(map(str,xy[i])) #[5,4] ->5x4
            
            #plotting all of the corner cuts on the same plot figure
            ax.errorbar(m, entropies[i,:], yerr = err[i,:], fmt = '.-', aa=True)
            ax.set_xlabel('Pathlength Parameter (m)')
            ax.set_ylabel('$S_2$')
            
            print "Plotting individual pathlength %s ..." %(corner_cut)
            #plot the individual plots
            ax2.errorbar(m,entropies[i,:], yerr= err[i,:],fmt = 'k.-', aa=True)
            ax2.plot(m, [Ann_Entropy[i,2]]*len(m), color = 'r', linewidth = 1.5 )
            ax2.set_xlabel('Pathlength Parameter (m)')
            ax2.set_ylabel('$S_2$')
            ax2.set_title(corner_cut + " Corner Cut")
            #ax2.set_xticks([1200], minor=True)            
            ax2.grid(axis = 'x', which='both')
            ax2.legend(("MC","Ann's $S_2$"),'lower right',fontsize = 'xx-small')
            f2.savefig(loc+ 'm_%s_cut%s.png' %(cluster,corner_cut), bbox_inches = 'tight', dpi = 500 )
            ax2.clear()
            print "Done individual ..."
            
            corner_cut_all.append(corner_cut) #store for large graph
            
        #final touches on the all corner cut plot
        ax.grid(axis = 'x')
        ax.legend(corner_cut_all, fontsize = 'xx-small')
        ax.set_title("All Corner Cuts")
        f.savefig('m_all_%s.png' %(cluster), bbox_inches = 'tight', dpi = 500)
        f.clf()
        
        #big plot doesn't show all the error bars due to scale, plot individually as well:
        #plt.errorbar(m, entropies[i,:], yerr = err[i,:], fmt = 'co-')
        #plt.xlabel('pathlength parameter')
        #plt.ylabel('$S_2$')
        #plt.legend(xy[i])
        #plt.savefig('_%s.png' %(loc.rstrip('/')))
        #plt.clf()
        
        print 'Done'
        
    plt.close('all')
        
def plot_all(combos, plt_n = True, plt_m = False):

    data_list = glob.glob("*/*/01.data") #find data in subfolder, dont care about order   
    
    #generic loop
    lattice_path_list = []
    for name in data_list: 
        lattice_size = map(int, name.split('/')[0].split('x'))
        ind = N.index(lattice_size)
        save_path = name.rstrip('01.data')
        analyze(combos[ind], path = save_path)
        
        #make a list of lattice paths
        if name.split('/')[0]+'/' not in lattice_path_list:
            lattice_path_list.append(name.split('/')[0]+ '/')
        
        #plot all subdirectories of n
        if plt_n == True: #plot bin vs entropy
            plot_n(path = name.rstrip('01.data'))

    if plt_m == True:
        plot_m(path = lattice_path_list)

        
if __name__ == '__main__':
    
    exe = 'tfim.out'
    #let's make directories for everything
    N = [[5,4], [7,3], [10,2]]
    m = range(100,3901, 100) + range(4000, 5001, 1000) + [1190,1210]
    path = os.getcwd() + '/'
    
    #loop over
    command_list = []
    list_combos = []
    for i in range(len(N)):
        #make directory
        base = path + "%ix%i" %(N[i][0], N[i][1])
        params = np.loadtxt("param.dat")
        params[0] = N[i][0]
        params[1] = N[i][1]
        
        #if file exists, continue
        try:
            os.mkdir(base)
             #do param.dat stuff
            np.savetxt(base+"/param.dat",params, fmt = "%i")
            list_combos.append(prep_A(base, new_A = True))
            #make region file
        except OSError: 
            list_combos.append(prep_A(base, new_A = False))       
        
        for j in m: #loop over pathlength vars
            #make dir
            m_path = base + "/m_%i" %(j)
            
            try:
                os.mkdir(m_path)
                
                params[-2] = j
                    
                #write new param file
                np.savetxt(m_path+ '/param.dat',params, fmt = '%g')
                #copy over regionA.dat, tfim.out
                copy('tfim.out', m_path + "/tfim.out")
                copy(base + "/regionA.dat", m_path + "/regionA.dat")
                copy(path + "regionX.dat", m_path+"/regionX.dat" )
                #run all processes
                os.chdir(m_path)
                command_list.append(Popen('./'+exe))
            except OSError:
                continue                   
        
    if command_list != []: #run commands if not empty
        #wait for stuff to finish
        [com.wait() for com in command_list]
    
    plot_all(list_combos, plt_n = False, plt_m = True)
    
    #run the simulation ----- pre upgrade for m tests
    
    #call('./' + exe)
    
    #put data in file
    #analyze(list_combos)
    #plot_m()--------------------

