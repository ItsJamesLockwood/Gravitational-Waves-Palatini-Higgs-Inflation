# -*- coding: utf-8 -*-
"""
Created on Tue Mar 23 18:20:32 2021

@author: James
"""

from sys import getsizeof
from plotUtils import *
from time import time,sleep

p1 = r"D:\Physics\MPhys Project\DatasetArcive\noMpl_tach_1_1024_slices_screen.log"
p2 = r"D:\Physics\MPhys Project\DatasetArcive\Remote tests\rsimp-lf4-run1_screen.log"

path = p1
file = open(path,'r+')
print("Size of file:", getsizeof(file))

def access_line(line_number,path,block_size=3,preamble=1):
    file = open(path,'r')
    if type(line_number)==int:
        max_read_out = line_number
        range_number = 1
    else:
        #Range: (first_line, last_line)
        max_read_out = line_number[0]
        range_number = line_number[1]-line_number[0]
    print("Range: ",max_read_out,range_number)
    # Skip preamble
    for _ in range(preamble):
        file.readline()
    # Skip n-1 lines
    print("Acessing range of interest...")
    for i in range(max_read_out-1):
        if i%500==0:
            print("\rSkipping... Progress: %i %%"%(i/max_read_out*100),end='',flush=True)
        for _ in range(block_size):
            file.readline()
    print("\rSkipping... Progress: 100%",end='',flush=True)
    print("\nDone.")
    #Read out range
    lines = []
    for _ in range(range_number):
        l = file.readline()
        if l=='':
            break
        lines.append(l)
    #Close file
    file.close()
    if len(lines)==0:
        print("Warning: there were no empty lines in the file. Are you sure the range is correct?")
    return lines
    
        
def reduce_slices(path,field_number=1,skip=50,mode='f',count_lines=True):
    inpath = trim_name(path) + "_slice_%s_%i.log"%(mode,field_number) 
    outpath = trim_name(path) + "_COPY%i_slice_%s_%i.log"%(skip, mode,field_number)
    
    file_in = open(inpath,'r')
    file_out = open(outpath,'w')
    
    #Copy over the header (resolution of simulation) to new file.
    l_res = file_in.readline()
    file_out.write(l_res)
    
    l1 = file_in.readline()
    l2 = file_in.readline()
    l3 = file_in.readline()
    counter = 0
    while l1!='':
        if counter%skip==0:
            print('\rProcessing line: %i'%counter,end='',flush=True)
            file_out.write(l1)
            file_out.write(l2)
            file_out.write(l3)
        l1 = file_in.readline()
        l2 = file_in.readline()
        l3 = file_in.readline()
        counter +=1
    file_in.close()
    file_out.close()
    return counter


lines = reduce_slices(path,mode='f',skip=200)
copied_file = r"D:\Physics\MPhys Project\DatasetArcive\noMpl_tach_1_1024_slices_slice_p_1.log"

#my_range = access_line([10500,10800],copied_file)
    
        