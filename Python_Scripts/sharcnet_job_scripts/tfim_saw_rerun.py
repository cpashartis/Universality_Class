#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 16:17:54 2014

@author: christoforos

This script runs the specified directories over.

Usage:

-a, --all -> Runs all tfim.out in directories found recursively

-i, --include '3x3 5x5 6x6' -> runs the directories given

-e, --exclude '4x3 5x3 6x5' -> runs everything but the directories given
"""

import argparse
import sys
from os import walk, getcwd, system, chdir,remove, path


def edit_param(bins, path = ''):
    
    params_in = open(path + '/param.dat', 'r')
    params = params_in.readlines()
    params_in.close()
    
    remove('param.dat')
    param_save = open(path + '/param.dat', 'w')
    
    counter = 1
    for line in params:
        if counter == 3:
           param_save.write('0\n')
        elif counter == 5 :
            param_save.write(str(bins) + '\n')
        else:
            param_save.write(line)
        
        counter +=1
    param_save.close()
    
    
    
parser = argparse.ArgumentParser(description='Rerun jobs using defined satistics')
parser.add_argument("-a", "--all", default = False, action = 'store_true',
                    help="runs all the files in the directory")
parser.add_argument("-i", "--include", nargs='+', default = [],
                  help="include the given file directories such as 3x3")
parser.add_argument("-e", "--exclude", nargs='+',default = [],
                  help="exclude the the given file such as 3x3")
parser.add_argument("-b", "--bins", nargs = 1, default = 90000, 
                    help = "define the number of bins to rerun" )
parser.add_argument("-t", "--time", nargs = 1, default = '7d', 
                    help = "define the length of the qsub command" ) 
parser.add_argument("-r", "--restart", default = False, action = 'store_true', 
                    help = "restart all runs in directory" )   
options = parser.parse_args()
#options is a dictionary

print options
if raw_input("Would you like to run this setup? y") != 'y':
    sys.exit('You decided not to run')
    
#now to do the combinations
#first check for all /exclude
cwd = getcwd()
if options.all == True and options.restart == False:
    for root, dirname, filename_list in walk(cwd): 
        #recursively finds all files in a directory
        for filename in filename_list:
            if filename == '01.data':
                if root.split('/')[-3] not in options.exclude and root.split('/')[-2] not in options.exclude:
                    chdir(root)
                    edit_param(options.bins, path = root)
                    system('sqsub -r %s -q NRAP_893 --mpp 500MB -o tfim.log ./tfim.out' %(options.time))
                    chdir(cwd)
elif options.all == False:
    if options.include != []:
        for dir_inc in options.include:
            for root, dirname, filename_list in walk(cwd + '/' + dir_inc ):
                for filename in filename_list: 
                    if filename == '01.data':

                        chdir(root)
                        
                        if options.restart == True and path.isfile('00.data') == False: #restart, just rerun
                            system('sqsub -r %s -q NRAP_893 --mpp 500MB -o tfim.log ./tfim.out' %(options.time))                            
                        elif options.restart == False:
                            edit_param(options.bins, path = root)
                            system('sqsub -r %s -q NRAP_893 --mpp 500MB -o tfim.log ./tfim.out' %(options.time))
                            
                        chdir(cwd)
    else:
        sys.exit('No options found, exiting')
    
