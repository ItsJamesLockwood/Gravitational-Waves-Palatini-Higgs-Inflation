Integer(IB),parameter :: n = 32
#if USE_STANDARD_WAVENUMBER
Integer(IB),parameter ::  fft_numk = 28
#if DIS_SCHEME == LATTICEEASY
Integer(IB),parameter ::  fft_cutoff_ind = 16
#elif DIS_SCHEME == HLATTICE1
Integer(IB),parameter ::  fft_cutoff_ind = 8
#elif DIS_SCHEME == HLATTICE2
Integer(IB),parameter ::  fft_cutoff_ind = 10
#endif
#elif DIS_SCHEME == LATTICEEASY
Integer(IB),parameter ::  fft_numk = 18
Integer(IB),parameter ::  fft_cutoff_ind = 16
#elif DIS_SCHEME == HLATTICE1
Integer(IB),parameter ::  fft_numk = 9
Integer(IB),parameter ::  fft_cutoff_ind = 8
#elif DIS_SCHEME == HLATTICE2
Integer(IB),parameter ::  fft_numk = 12
Integer(IB),parameter ::  fft_cutoff_ind = 10
#endif
Integer(IB),parameter :: sind(-64:96) = (/ & 
32, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13,  & 
14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30,  & 
31, 32, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,  & 
16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32,  & 
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17,  & 
18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 1, 2,  & 
3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,  & 
20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 1, 2, 3, 4,  & 
5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21,  & 
22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32 /) 
Integer(IB),parameter :: bit_reverse_array(n) = (/ & 
1, 17, 9, 25, 5, 21, 13, 29, 3, 19, 11, 27, 7, 23, 15, 31, 2,  & 
18, 10, 26, 6, 22, 14, 30, 4, 20, 12, 28, 8, 24, 16, 32 /)
Integer(IB),parameter :: fft_conj_index(n) = (/ 1, & 
32, 31, 30, 29, 28, 27, 26, 25, 24, 23, 22, 21, 20, 19, 18,  & 
17, 16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3,  2 /)
Integer(IB),parameter :: fft_rawk_index(n) = (/0,  & 
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,  & 
16, -15, -14, -13, -12, -11, -10, -9, -8, -7, -6, -5, -4, -3, -2,  -1 /)
Integer(IB),parameter :: fft_rawk_index_abs(n) = (/0,  & 
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,  & 
16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2,  1 /)
