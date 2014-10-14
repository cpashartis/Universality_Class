# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Thu Jun  5 10:11:36 2014

@author: christoforos
"""

from matplotlib import use 
use('agg')
from matplotlib import rc_file
rc_file('original.rc')

import matplotlib.pyplot as plt
import numpy as np
import MCstat


def plot_n(path = '', Ann_check = ''):
    
    if path[-1] != '/':
        path += '/'
    
    data = np.loadtxt(path + '01.data')
    x = np.arange(len(data))
    f = plt.figure()
    ax = f.add_subplot(111)
    
    configuration = path.strip('/').split("/")[-1]
    print "Plotting bin plots %s ... " %(path)
    
    try:
        for dim in range(len(data[0,:])):
            y =-1*np.log(data[:,dim])
            #plot over m, i.e. pathlength
            plt.plot(x, y,'k.')
            
            #linear fit the plot
            a,b = np.polyfit(x, y, 1)
            
            plt.plot(x, a*x + b, 'r')
            plt.text(.5,.025,"$y= %.4e x + %.4e$" %(a,b), fontsize = 12,
                     transform=ax.transAxes)
            plt.xlabel("Bin Number")
            plt.ylabel("$S_2$")
            plt.title(configuration)
            plt.savefig(path + configuration + ' plot_n%i.png' %(dim+1),
                        bbox_inches = 'tight')
            plt.clf()
            
    except IndexError:
        y =-1*np.log(data[:])
        #plot over m, i.e. pathlength
        plt.plot(x, y,'k.')
        plt.xlabel("Bin Number")
        plt.ylabel("$S_2$")
        plt.title(path.split("/")[-2])
        plt.savefig(path + configuration + ' plot_n.png',
                    bbox_inches = 'tight')
        
    print "Done"
    plt.close()

def Re_Bin(path ='', plot = True):
    """This function uses MCstat to rebin all data,
    plots and returns the data"""
    
    if path[-1] != '/':
        path +='/'
    
    bins_o = np.loadtxt(path + "01.data")

    bin_err = MCstat.bin(bins_o)
    
    if plot == True:
        #get names
        tags = open(path+ "01.data", "r").readline().lstrip("#").split()
        #replace("_", "-").split() #issue with latex and using '_'
        #begin plotting subroutine
        fig = plt.figure()
        ax = fig.add_subplot(111)
        bin_level = range(len(bin_err))
        #loop over all cut names
        
        print "Plotting ReBining: %s ..." %(path)
        
        #check for data to make sure data stored in each row
        for i in range(len(tags)):
            ax.plot(bin_level, bin_err[:,i], label = tags[i], linewidth = 2)
        
        ax.legend(loc = 'upper left',fontsize = 'xx-small')
        ax.set_title("Re-Bin of raw %s" %(path.strip('/').split('/')[-1]))
        ax.set_xlabel("Bin Level")
        ax.set_ylabel("Error")
        plt.savefig(path + 'ReBin_Error.png')
        plt.close() #removes figure/canvas 
        
        print "Done"
    
    #returns the maximum error and its associated mean (without nasty loops)
    index = bin_err.argmax(axis = 0)
        
    return bin_err[index, range(len(index))]
 
def plot_universality(path, order_type, hybrid = False):
    """This function takes the output from All_Results and Original_Results 
    and generates 3 plots. These plots are the unfitted data, the fitted data with
    the slope coefficient extracted (the error of which is a function of the covariance matrix),
    and a plot of the error as a function of the number of iterations of the gaussian noise.
    
    Notes:
    -> order_type is just for graph labels
    -> Hybrid allows for partial weighted fits
    """
    
    # import dependencies
    rc_file('latex_nice.rc')
    import Lauren_stat
    from scipy import polyfit
    
    print 'Plotting universality graphs for : ' + order_type
    
    #check if plot
    #begin plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
#------------------Import data-------------------------------------------------
    data = np.loadtxt(path + "All_Results")
    o_data = np.loadtxt(path + "Original_Results")
    nlce_avg, nlce_std = data.mean(axis = 0), data.std(axis=0)
    order = np.array(open("All_Results", "r").readline().strip('#').split(), 
                  dtype = float)
    
    #define fitting function
    fit = lambda x, a ,b : a*np.log(x) +b
###############################################################################
#FIRST PLOT, order NLCE
#--------------fit non weighted and original-----------------------------------
    coef2 = polyfit(np.log(1/order), nlce_avg, 1)
    #o_coef = polyfit(np.log(1/order), o_data,1)    
    order_lin = np.linspace(order.min(), order.max())
    ax.errorbar(1/order, nlce_avg, yerr = nlce_std, fmt = 'o',
                label = 'Gaussian Propagated Error' )
    #ax.plot(1/order, o_data, 'o', label = 'Original (Pre-Gaussian)')
    
#--------------fit weighted-----------------------------------    
    #first check if zombie, to remove exact data
    if hybrid == True :
        #update order_type
        index_cut = 0
        for element in nlce_std:
        
            if element >= 10**-9: #tolerance for zero
                break
            index_cut +=1
        #update parameters to cut out exact
        nlce_avg = nlce_avg[index_cut:]
        nlce_std = nlce_std[index_cut:]
        order = order[index_cut:]
        o_data = o_data[index_cut:]
        data = data[index_cut:]
        order_lin = np.linspace(order.min(), order.max())
    
    #now that data is fixed let us fit it
    coef, err, reduced, cov = Lauren_stat.fit_leastsq(1/order, nlce_avg, 
                    nlce_std, fit, Lauren_stat.residuals_weighted,[0,0])
    #plot weighted fit                
    ax.plot(1/order_lin, fit(1/order_lin, *coef), '-', 
            label = 'Weighted Fit, $\chi_{reduced}^2$ %.3f' %(reduced) )
    
#-------------------Add Ann's data to order plot------------------------------
    ann_data = np.loadtxt('/Users/christoforos/NSERC_Melko/QMC/Ann_Check/order_g.dat') #first col order, second is entropy
    ax.plot(1/np.sqrt(ann_data[:,0]), ann_data[:,1]/2., 'ko-',label = "Ann's Geometric Data.data")
    
    ax.legend(loc = 'lower right')
    ax.set_xlabel(r"$\frac{1}{l}$")
    ax.set_ylabel(r"$c_{2}$")
    ax.set_title("NLCE $O_{%s}$" %(order_type[0]))
    plt.savefig("Raw_Data_%s.png" %(order_type))
    ax.clear()

###############################################################################
#SECOND PLOT, corner entropy term
    #load data for the graph of all alpha
    tfim = np.loadtxt('/Users/christoforos/NSERC_Melko/QMC/Ann_Check/tfim.dat')
    
    #plot the alpha curve
    ax.plot(tfim[:,0],tfim[:,1], '#990000', label = r"$O_G$ coefficients")
    ax.plot(tfim[:,0],tfim[:,1]+tfim[:,2], '#A14444', tfim[:,0],tfim[:,1]-tfim[:,2], '#A14444', linewidth =.5)
    ax.fill_between(tfim[:,0], tfim[:,1]+tfim[:,2], tfim[:,1]-tfim[:,2], color = '#FFA8A8')
    ax.set_ylim(0.003,0.018)
    ax.set_xlim(0.85,3)
    
    ax.errorbar(2, coef[0], yerr = err[0], fmt = 'o', label = "QMC, weighted fit, error =  %.5e" %(err[0]))
    new_file = open('Store_Universality_Class', 'a')
    new_file.write(order_type + '\t%.10f\t%.10f\n' %(coef[0],err[0]) )
    new_file.close()
    #plot(2, o_coef[0], 'o', label = 'QMC, original fitted' )
    ax.plot(2, coef2[0], 'o', label = 'QMC, non-weighted fit')
    
    ax.set_ylabel(r"$-a_{\alpha}$")
    ax.set_xlabel(r"$\alpha$")
    ax.set_title("Universality Class Number $O_{%s}$" %(order_type[0]))
    ax.legend()
    plt.savefig('Fitted_Term_Alphas_%s.png' %(order_type), dpi = 500)
    ax.clear()
    
###############################################################################
#THIRD PLOT, plot of error
    nlce_error = []
    for j in range(1,len(data)):
        nlce_error.append(data[:j,:].std(axis=0))
    
    #convert to array
    nlce_error= np.array(nlce_error)

    ax.plot(range(1,len(data)), nlce_error)
    ax.set_xlabel("Iteration Number")
    ax.set_ylabel("Error $(\sigma)$")
    plt.savefig('Gaussian_Error_%s.png' %(order_type), dpi = 500)
    
    plt.close() #dont clober figures   
    
    print """Plotting done! Nothing to see here...
    Move along....
    Move along..."""

"""
def plot_m(paths = ''):
    
    
    data = []
    for direc in paths:
        data.append(np.loadtxt(direc))
    
    for row in range(len(data[0])): #loop over all paths given, (generic handling, yay!)
        
        entropies = []
        err = []
        m = []
        
        for pathlength_ind in range(len(data))
            m.append(m_name.split('/')[1].lstrip('m_'))  
            data = np.loadtxt(m_name)            
            entropies.append(data[:,2])
            err.append(data[:,3])
        for m_name in glob.glob(loc+'*/R_Entropy.data'): #load in all data from entropy by searching for the files
            #append data to lists for now            
            
        
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

if __name__ == '__main__':
    
    cwd = getcwd()
    for root, dirname, filename_list in walk(cwd): 
    #recursively finds all files in a directory
        for filename in filename_list:
            if filename == '01.data':
                #plot_n(path = root)
                #Re_Bin(path = root)
                plot_m()"""