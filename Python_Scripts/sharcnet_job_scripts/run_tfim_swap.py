# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Wed May 14 15:32:02 2014

@author: christoforos

Brief description:

This script eats dakota input and creates directories for various m, n, n_skip
runs. The script then spits out the second renyi entropy and the error for dakota use.

"""

#user stuff###
#1 - Geometric
#2 - Arithmetic

order_min = 4 #1
#order_min = 2 #1
order_max = 36
#order_max = 12
order_type = 'Geometric' #1
#order_type = 'Arithmetic' #2
#order_min = 2 #3
#order_max = 100 #3
#order_type = 'Quadratic' #3
seed_num = 1809237
user_swap = True
################

import numpy as np
import os
from sys import exit
import shutil
import mxn_order
from numpy import ceil

class run_tfim:
    
    def __init__(self, bins, seed, L_x, L_y, swap = False):
        
        self.m = L_x*L_y*125
        self.bins = bins
        self.seed = seed
        self.L_x = L_x
        self.L_y = L_y
        self.mcs = L_x*L_y*500 #find a way to generalize
        self.eql = self.bins*self.mcs/5000
        self.swap = swap
        
        #check for symmetry
        if self.L_x == self.L_y:    
            self.symmetry = True
            
            self.num_of_cuts = int(sum(range(1,L_x-1)))+ L_x-1 +int((L_x+1)/2)
            #sum of range is the number of upper triangle corner cuts
            #L_X-1 term is the number of symmetric squares from corner
            #third term is the number of vertical/horizontal symmetry cuts
        else:
            self.symmetry = False 
            self.num_of_cuts = (L_x-1)*(L_y-1) + int((L_y+1)/2) + int((L_x+1)/2)
        
###############################################################################
        
    def __print__(self):
        
        print np.array([[self.L_x],[self.L_y], [self.eql], [self.mcs],
                           [self.bins], [self.seed], [self.m], [3.044]])

###############################################################################

    def tag_gen(self,x,y, x_reduced = 0, y_reduced = 0):
        
        if y == self.L_y:
            if self.symmetry == True:
                tag = 'L%.2i%.2i_%.2i' %(self.L_x,self.L_y, x)
            else:
                tag =  'L%.2i%.2i_%.2iY' %(self.L_x,self.L_y, x)  
        elif x == self.L_x:
            if self.symmetry == True:
                tag = 'L%.2i%.2i_%.2i' %(self.L_x,self.L_y, y)
            else:
                tag =  'L%.2i%.2i_%.2iX' %(self.L_x,self.L_y, y)
        else: #non zero x and y
            tag = 'C%.2i%.2i_%.2i%.2i' %(self.L_x,self.L_y, x, y)
        
        if x_reduced != 0:
            if y_reduced !=0: #shouldn't happen as far as I know
                tag += '_red%.2i%.2i' %(x_reduced,y_reduced)
            elif y == 1:# for corner cuts of single line, re-use corner cuts of single line
                tag = 'C%.2i%.2i_%.2i%.2i' %(self.L_x,self.L_y, x_reduced, y)
            else:
                tag += '_red%.2iX' %(x_reduced)
        elif y_reduced !=0 :                
            tag += '_red%.2iY' %(y_reduced)
        
        return tag       
            
###############################################################################
    
    def cut_gen(self):
        
        try:
            from collections import OrderedDict
            tag_dic = OrderedDict({})
        except ImportError:
            tag_dic = {}
        
#-------first find all corner cuts        
        for y in range(1, self.L_y): #loops over ys
            if self.symmetry == True:
                lower = y
            else:
                lower = 1
            #loops over xs in ys, excludes same types for symmetry
            for x in range(lower,self.L_x): #doing this topleft down, shouldnt matter
                tag_dic.update({self.tag_gen(x,y):(x,y)})

#-------horizontal cut
        x = self.L_x
        for y in range(1, int((self.L_y)/2)+1):
            tag_dic.update({self.tag_gen(x,y):(x,y)})
#-------vertical cuts
        y = self.L_y
        if self.symmetry == False:  #need to write the vertical cuts as well
            for x in range(1,int((self.L_x)/2)+1):
                tag_dic.update({self.tag_gen(x,y):(x,y)})
        
        return tag_dic 
            
###############################################################################
    
    def write_generic(self,x,y, x_reduced=0, y_reduced=0):
        
        #check for reduced 0
        if x_reduced == 0:
            x_reduced = self.L_x
        if y_reduced == 0:
            y_reduced = self.L_y
        self.write_file.write("\n")
        for row in range(1,self.L_y+1):
            if row <= y:
                if row == y and x_reduced != self.L_x:#modify row for row cut
                    self.write_file.write(" ".join(map(str,[1]*(x_reduced) + [0]*(self.L_x-x_reduced)))+'\n')
                elif row > y_reduced:
                    self.write_file.write(" ".join(map(str,[1]*(x-1) + [0]*(self.L_x-x+1)))+'\n')
                else:
                    self.write_file.write(" ".join(map(str,[1]*x+[0]*(self.L_x-x)))+'\n')
            else:
                self.write_file.write(" ".join(map(str,[0]*self.L_x))+'\n')
        self.write_file.write("-99")
        
###############################################################################        
    def swap_gen(self,x,y, build_list = [], xy_list = []):
        
        #for some reason doing [:] prevents errors ...
        #buffer variables
        build_list_dummy = build_list[:]
        xy_list_dummy = xy_list[:]
        #number of allowed difference in sites
        num_sites = 8
        
        if build_list_dummy == []:
            #find original name
            build_list_dummy.append(self.tag_gen(x,y))
            xy_list_dummy.append((x,y,0,0))
        
        #first check size
        if x+y-1 > num_sites : #perimeter
            #need in between cuts
            if x >= y: #determine how many cuts are required
                num_cuts = ceil(float(x)/num_sites)
                cut = int(ceil(x/num_cuts))
                row = x
                
                if y ==1 and x>1: #check for single line, reduced to rest of corner
                    x = x-1
                    build_list_dummy.append(self.tag_gen(x,y))
                    xy_list_dummy.append((x,y,0,0))
                else:
                    while row > 0:
                        row = row - cut
                        if row <= 0 :
                            y = y-1
                            row = 0
                        build_list_dummy.append(self.tag_gen(x,y,x_reduced = row))
                        xy_list_dummy.append((x,y,row,0))
            else: #just use columns
                num_cuts = ceil(float(y)/num_sites)
                cut = int(ceil(y/num_cuts))
                col = y
                while col > 0:
                    col = col - cut
                    if col <= 0 :
                        if x == 1:
                            y = 1
                        else:
                            x = x-1
                        col = 0
                    build_list_dummy.append(self.tag_gen(x,y,y_reduced = col))
                    xy_list_dummy.append((x,y,0,col))   
        else: #just use previous cuts
            if y > x:                    
                y = y-1
            elif x > y:
                x = x-1
            elif x > 1 and y >1 : #i.e. equal and thus a symmetrical corner
                x = x-1
                y = y-1
            else: #return the final list of names
                return build_list_dummy, xy_list_dummy
                
            build_list_dummy.append(self.tag_gen(x,y))
            xy_list_dummy.append((x,y,0,0))
    
        #recursive
        return self.swap_gen(x,y,build_list = build_list_dummy, xy_list = xy_list_dummy)
        
###############################################################################
        
    def make_dirs(self, op_dir, overwrite = True, verbose = False):
        #returns true if made a file
    
        if os.path.isdir(op_dir) == False:
            os.makedirs(op_dir) #make directories(including new subfolders)
            return False
        else : #file exists
            #if statement to handle existing directories and potential files
            if overwrite == True:
                shutil.rmtree(op_dir)
                os.makedirs(op_dir)
                if verbose == True:
                    print "Deleting and overwriting ..."
                return False
            else:
                if verbose == True:
                    print "Directory not removed...continuing"
                return True
        
    def Generate(self):
    
        op_dir = "%ix%i/eql%i_mcs%i_nbins%i_m%i" %(self.L_x, self.L_y, 
                                                   self.eql, self.mcs,
                                                   self.bins, self.m)
                                    
#--------now that directory is made, let's populate it with the necessary files
        #making param.dat file
        params = np.array([[self.L_x],[self.L_y], [self.eql], [self.mcs],
                           [self.bins], [self.seed], [self.m], [3.044]])
        
        
        #generate cut names for no swap
        if self.swap == False:
            #make directory
            self.make_dirs(op_dir)
            #copy tfim.out
            shutil.copy('tfim.out', op_dir+ '/tfim.out')
            #chdir need to run ./tfim.out
            os.chdir(op_dir)
            np.savetxt("param.dat", params, fmt = '%.16g') #write as type in array and allow 16 digits
            
            self.write_file = open("regionA.dat", "w")
            dat_01 = open("01.data", "w")
            self.write_file.write("%i" %(self.num_of_cuts)) #write # of regions
            
            cuts = self.cut_gen()
            dat_01.write('#' + "\t".join(cuts.keys()) + '\n') #write cut names for later
            for x,y in cuts.itervalues():
                self.write_generic(x,y)
            
            #close files
            self.write_file.close()
            dat_01.close()
            #make regionX file
            regionX = open("regionX.dat", "w")
            regionX.write("0")
            regionX.close()
            
            print op_dir + " prepared ..."
            os.system('sqsub -r 7d -q serial --mpp 500MB -o tfim_saw.log ./tfim.out') #run command for sharcnet   
            #os.popen("./tfim.out")
            os.chdir('../..') #back to original directory
        
        else:
            original_cuts = self.cut_gen()
            
            self.make_dirs(op_dir, overwrite = False)
            swap_file = open(op_dir + "/Swap_Progression.dat", 'w')
            
            for cut_x,cut_y in original_cuts.itervalues():
                #now generate the swap boxes with the files
                swap_tag, swap_val = self.swap_gen(cut_x,cut_y, )
                #print cut_x, cut_y, swap_tag, swap_val
                counter = 0
                
                #write swap progression to file
                swap_file.write('\t'.join(swap_tag)+ '\n')
                
                for counter in range(len(swap_val)):
                    #time to make the files
                    #first add to original directory
                    new_op_dir = op_dir + '/%s' %(swap_tag[counter])
                    existing = self.make_dirs(new_op_dir, overwrite = False)
                    
                    #if file exists don't rerun it, need this since swap generator gives full list
                    if existing == False:
                        
                        shutil.copy('tfim.out', new_op_dir+ '/tfim.out')
                        #chdir need to run ./tfim.out
                        os.chdir(new_op_dir)
                        np.savetxt("param.dat", params, fmt = '%.16g') #write as type in array and allow 16 digits
                        
                        #regionA.dat                    
                        self.write_file = open("regionA.dat", "w")
                        self.write_file.write('1') #write # of regions
                        self.write_generic(swap_val[counter][0],
                                               swap_val[counter][1],
                                            x_reduced = swap_val[counter][2],
                                            y_reduced= swap_val[counter][3])
                        self.write_file.close()
                        
                        #regionX.dat
                        self.write_file = open("regionX.dat",'w')
                        if counter == len(swap_val)-1:
                            self.write_file.write('0')
                        else:
                            self.write_file = open("regionX.dat",'w')
                            self.write_file.write('1')
                            self.write_generic(swap_val[counter+1][0],
                                               swap_val[counter+1][1],
                                            x_reduced = swap_val[counter+1][2],
                                            y_reduced= swap_val[counter+1][3])
                        self.write_file.close()
                        
                        #now do the 01.data header for consistancy/compatibility
                        dat_01 = open("01.data", "w")
                        dat_01.write('#%s' %(swap_tag[counter]) + '\n') #write cut names for later
                        dat_01.close()
                        
                        print new_op_dir + " prepared ..."
                        os.system('sqsub -r 7d -q serial --mpp 500MB -o tfim.log ./tfim.out') #run command for sharcnet   
                        #os.popen("./tfim.out")
                        os.chdir('../../..') #back to original directory
                    else:
                        print new_op_dir + " already generated and run, continuing ..."
                        #continue with loop
                        
            swap_file.close()
        
if __name__ == '__main__':
    for I in range(order_min,order_max+1):
        
        if order_type == 'Arithmetic':
            order = mxn_order.Arithmetic()
        elif order_type == 'Geometric':
            order = mxn_order.Geometric()
        elif order_type == 'Quadratic':
            order = mxn_order.Quadratic()
        else:
            exit("QUIT: Order type not defined")
        for m,n in order.clusters(I):
            run_tfim(10000, seed_num, m, n, swap = user_swap).Generate()
