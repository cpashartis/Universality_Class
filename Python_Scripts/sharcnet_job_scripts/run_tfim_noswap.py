# -*- coding: utf-8 -*-
#!/usr/bin/python
"""
Created on Wed May 14 15:32:02 2014

@author: christoforos

Brief description:

This script eats dakota input and creates directories for various m, n, n_skip
runs. The script then spits out the second renyi entropy and the error for dakota use.

"""
import numpy as np
import os
import shutil


#user stuff###
order_min = 2
order_max = 10
################

class Arithmetic:
    """Pulled from Ravi's code, it will return the lattice sizes for a given order"""
    def length(self,N):
        return (2+N)/2.

    def lengthstr(self,N):
        return "%.1f"%self.length(N)

    def clusters(self,N):
        x = N
        y = 2
        res = [(x,y),]
        while x > (y+1):
            x -= 1
            y += 1
            res.append((x,y))
        return res
        
def run_all():
    
    cwd = os.getcwd()
    for root, dirname, filename_list in os.walk(cwd): #recursively finds all files in a directory
        for filename in filename_list:
            if filename == '01.data':
                os.chdir(root)
                os.system('sqsub -r 1d -q serial -o tfim_saw.log ./tfim.out')
                os.chdir(cwd)
                
def directory_manage(m, n , n_skip, n_bins, seed,cut_tag = '2DNLCE_0404' ):
    
    if '2DNLCE_' not in cut_tag: #all commands
        L_x = int(cut_tag[1:3])
        L_y = int(cut_tag[3:5])
        cut = cut_tag.split("_")[-1]
        if int(cut[:-1]) == 0:
            print "Invalid cut_tag"            
            return 0
    else:
        L_x = int(cut_tag.replace('2DNLCE_', '')[:2])
        L_y = int(cut_tag.replace('2DNLCE_', '')[2:])
        cut = '2DNLCE'
        
    #ANYWAYS
    #first thing is first, must build directory structure
    #use format 4x4/m#_n#_nskip#
    op_dir = "%ix%i/eql%i_mcs%i_nbins%i_m%i" %(L_x, L_y, n_skip, n, n_bins, m)
    try:
        os.makedirs(op_dir) #make directories(including new subfolders)
        
    except OSError: #file exists
        #check if the rest of the files are in here
        if raw_input("""The %s directory exists, would you like to \
        overwrite (y) ?""" %(op_dir)) == 'y': #if the data file exists
            #remove directory and make a new one            
            shutil.rmtree(op_dir)
            os.makedirs(op_dir)
        else:
            print "Directory not removed...continuing"
            return
            
        
    
    #now that directory is made, let's populate it with the necessary files
    
    #copy tfim.out
    shutil.copy('tfim.out', op_dir+ '/tfim.out')
    #chdir
    os.chdir(op_dir)
    
    #param.dat file
    params = np.array([[L_x],[L_y], [n_skip], [n], [n_bins], [seed], [m], [3.044]])
    np.savetxt("param.dat", params, fmt = '%.16g') #write as type in array and allow 16 digits
    
    #make regionA file
    regionA_make(L_x, L_y, cut)
    
    #make regionX file
    regionX = open("regionX.dat", "w")
    regionX.write("0")
    regionX.close()
    
    #run it
    print op_dir
    #cmd = Popen('./tfim.out') #run normally
    os.system('sqsub -r 1d -q serial -o tfim_saw.log ./tfim.out') #run command for sharcnet
    #cmd.wait()
    
    os.chdir('../..') #back to original directory
    
    return op_dir
    
#make regionA_file  only y cuts right now  
def regionA_make(L_x, L_y, cut):
    
    regionA = open("regionA.dat", "w")
    dat_01 = open("01.data", "w")
    #check for symmetry:
    if L_x == L_y:    
        symmetry = True
        
    else:
        symmetry = False
        
    
    if cut == '2DNLCE': #will want all corner cuts from bottom left as well as vertical and horizontal    
        
        if symmetry == True:
            num_of_cuts = int(sum(range(1,L_x-1)))+ L_x-1 +int((L_x+1)/2)
            
            #sum of range is the number of upper triangle corner cuts
            #L_X-1 term is the number of symmetric squares from corner
            #third term is the number of vertical/horizontal symmetry cuts
            
            #num_of_cuts = ((L_x-1)*(L_y-1)-L_x)/2 +L_x-1 + int(L_x/2) +1 
            #first term is for non symmetrical diagonals corners like 3x2 cut,
            #second is for adding symmetrical corner terms back
            #third term is for number of line cuts (removing symmetry hence /2)
        else:        
            num_of_cuts = (L_x-1)*(L_y-1) + int((L_y+1)/2) + int((L_x+1)/2) #check for symmetric and non symmetric
        
        regionA.write("%i" %(num_of_cuts)) #write # of regions
        
        #first find all corner cuts:
        tag_c = 'C%.2i%.2i_' %(L_x,L_y)
        temp = '#'
        for y in range(1, L_y):
            if symmetry == True:
                lower = y
            else:
                lower = 1
            for x in range(lower,L_x): #doing this topleft down, shouldnt matter
                temp += tag_c + '%.2i%.2i\t' %(x,y)
                row0 = [0]*L_x
                row1 = [1]*x + [0]*(L_x-x)
                regionA.write("\n")
                for row in range(L_y):
                    if row >=y:
                        regionA.write(" ".join(map(str,row0)) + '\n')
                    else:
                        regionA.write(" ".join(map(str,row1)) + '\n')
                regionA.write("-99")
        
        #for line cuts
        tag_l = 'L%.2i%.2i_' %(L_x,L_y)
        
        #vertical cut
        for x in range(1,int((L_x+1)/2)+1):
            if symmetry == True:
                temp += tag_l+ '%.2i\t' %(x)
            else:
                temp += tag_l+ '%.2iY\t' %(x)
            regionA.write("\n")
            line = [1]*x + [0]*(L_x - x )
            line = map(str, line)
            for y in range(L_y):
                regionA.write(" ".join(line) + '\n')
            regionA.write("-99")
        
        if symmetry == False:
            #horizontal cut
            for y in range(1, int((L_y+1)/2)+1): #doing the squares top left down, easier and shouldn't matter
                temp += tag_l+ '%.2iX\t' %(y)
                regionA.write("\n")
                row0 = [0]*L_x
                row1 = [1]*L_x
                for row in range(L_y):
                    if row >=y:
                        regionA.write(" ".join(map(str,row0)) + '\n')
                    else:
                        regionA.write(" ".join(map(str,row1)) + '\n')
                regionA.write("-99")
        temp = temp.rstrip("\t")
        dat_01.write(temp + '\n') #write headers for 01.data file for later use
            
########vertical cuts in
    else:           
        if cut[-1] == 'Y' : #vertical cuts
            regionA.write("1\n")
            line = [1]*int(cut[:-1]) + [0]*(L_x - int(cut[:-1]) )
            line = map(str, line)
            for y in range(L_y):
                regionA.write(" ".join(line) + '\n')
            regionA.write("-99")
        elif cut[-1] == 'X': #horizontal cuts
            regionA.write("1\n")
            col = [[0]]*(L_x - int(cut[:-1]) ) + [[1]]*int(cut[:-1])
            col = map(str, col)
            total = col *L_x
            for line in total:
                regionA.write(" ".join(line) + '\n')
            regionA.write("-99")
        regionA.close() 

#for 2DNLCE
if __name__ == '__main__':
    for I in range(order_min,order_max+1):
       
        #define order here
        order = Arithmetic()
        
        for m,n in order.clusters(I):
            mcs = m*n*5000
            pathlength = m*n*125
            directory_manage(pathlength,mcs,mcs/5, 1000, 2005967, cut_tag = '2DNLCE_%.2i%.2i' %(m,n))
        
###otherwise just cut_tag = 'L....'            
    
    
    
    
