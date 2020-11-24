Integer(IB),parameter :: n = 16
#if USE_STANDARD_WAVENUMBER
Integer(IB),parameter ::  fft_numk = 14
#if DIS_SCHEME == LATTICEEASY
Integer(IB),parameter ::  fft_cutoff_ind = 8
#elif DIS_SCHEME == HLATTICE1
Integer(IB),parameter ::  fft_cutoff_ind = 4
#elif DIS_SCHEME == HLATTICE2
Integer(IB),parameter ::  fft_cutoff_ind = 5
#endif
#elif DIS_SCHEME == LATTICEEASY
Integer(IB),parameter ::  fft_numk = 9
Integer(IB),parameter ::  fft_cutoff_ind = 8
#elif DIS_SCHEME == HLATTICE1
Integer(IB),parameter ::  fft_numk = 4
Integer(IB),parameter ::  fft_cutoff_ind = 4
#elif DIS_SCHEME == HLATTICE2
Integer(IB),parameter ::  fft_numk = 6
Integer(IB),parameter ::  fft_cutoff_ind = 5
#endif
Integer(IB),parameter :: sind(-32:48) = (/ & 
16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15,  & 
16, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,  & 
1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1,  & 
2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 1, 2,  & 
3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16 /) 
Integer(IB),parameter :: bit_reverse_array(n) = (/ & 
1, 9, 5, 13, 3, 11, 7, 15, 2, 10, 6, 14, 4, 12, 8, 16 /)
Integer(IB),parameter :: fft_conj_index(n) = (/ 1, & 
16, 15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3,  2 /)
Integer(IB),parameter :: fft_rawk_index(n) = (/0,  & 
1, 2, 3, 4, 5, 6, 7, 8, -7, -6, -5, -4, -3, -2,  -1 /)
Integer(IB),parameter :: fft_rawk_index_abs(n) = (/0,  & 
1, 2, 3, 4, 5, 6, 7, 8, 7, 6, 5, 4, 3, 2,  1 /)
