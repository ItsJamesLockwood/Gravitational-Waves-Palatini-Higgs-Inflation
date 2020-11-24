#! /usr/bin/env python
import sys
from math import ceil
def BitReverse(n,p):
    '''convert denary integer n to binary string bStr'''
    br = 0
    if n < 0:  raise ValueError, "must be a positive integer"
    for i in range(0,p):
        br = (n % 2) + (br << 1)
        n = n >> 1
    return br


p=int(sys.argv[1])
if (p>15):
    print "a cubic lattice with (2^16) grid points seems to be too large"
    sys.exit()
if (p<3):
    print "syntax: \n python genparam.py p \n here 2<p<13 is an integer"
    sys.exit()
n=2**p
istr="Integer(IB),parameter ::"
fp=open("n_params_"+str(p)+".h","w")
fp.write(istr + " n = " + str(n) + "\n")
fp.write("#if USE_STANDARD_WAVENUMBER\n")
fp.write(istr + "  fft_numk = " + str(int(ceil(n*0.86603-0.4999)))+"\n")
fp.write("#if DIS_SCHEME == LATTICEEASY\n")
fp.write(istr + "  fft_cutoff_ind = " + str(n/2) +"\n")
fp.write("#elif DIS_SCHEME == HLATTICE1\n")
fp.write(istr + "  fft_cutoff_ind = " + str(int(ceil(n*0.25))) +"\n")
fp.write("#elif DIS_SCHEME == HLATTICE2\n")
fp.write(istr + "  fft_cutoff_ind = " + str(int(ceil(n*0.29))) +"\n")
fp.write("#endif\n")
fp.write("#elif DIS_SCHEME == LATTICEEASY\n")
fp.write(istr + "  fft_numk = " + str(int(ceil(n*0.55133-0.4999)))+"\n")
fp.write(istr + "  fft_cutoff_ind = " + str(n/2) +"\n")
fp.write("#elif DIS_SCHEME == HLATTICE1\n")
fp.write(istr + "  fft_numk = " + str(int(ceil(n*0.27566-0.4999)))+"\n")
fp.write(istr + "  fft_cutoff_ind = " + str(int(ceil(n*0.25))) +"\n")
fp.write("#elif DIS_SCHEME == HLATTICE2\n")
fp.write(istr + "  fft_numk = " + str(int(ceil(n*0.37222-0.4999)))+"\n")
fp.write(istr + "  fft_cutoff_ind = " + str(int(ceil(n*0.29))) +"\n")
fp.write("#endif\n")
fp.write(istr + " sind(-"+str(2*n)+":" + str(n*3) + ") = (/ & \n")
for i in range(-2*n, 3*n) :
    fp.write(str((i-1)%n+1) + ", ")
    if( i % 17 == 0):
        fp.write(" & \n")
fp.write( str(n)+" /) \n")
fp.write(istr+" bit_reverse_array(n) = (/ & \n")
for i in range(0,n-1):
    fp.write(str(BitReverse(i,p)+1)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
fp.write(str(BitReverse(n-1,p)+1)+" /)\n")
fp.write(istr+" fft_conj_index(n) = (/ 1, & \n")
for i in range(2,n):
    fp.write(str(n+2-i)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
fp.write(" 2 /)\n")
fp.write(istr+" fft_rawk_index(n) = (/0,  & \n")
for i in range(2,n/2+2):
    fp.write(str(i-1)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
for i in range(n/2+2,n):
    fp.write(str(i-1-n)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
fp.write(" -1 /)\n")
fp.write(istr+" fft_rawk_index_abs(n) = (/0,  & \n")
for i in range(2,n/2+2):
    fp.write(str(i-1)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
for i in range(n/2+2,n):
    fp.write(str(n+1-i)+", ")
    if( (i+1) % 17 == 0):
        fp.write(" & \n")
fp.write(" 1 /)\n")


fp.close()
