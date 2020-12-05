#if USE_OMP
#ifdef KSPACE_OPERATION
#ifdef REDUCTION_VARS
    !$omp parallel do default(shared) private(i,j,k,i1,j1,k1,i2,j2,k2) reduction(+:REDUCTION_VARS)
#undef REDUCTION_VARS
#else
    !$omp parallel do default(shared) private(i,j,k,i1,j1,k1,i2,j2,k2)
#endif
    do k=2,Nby2
       call fft_conj_pair(k,k1,k2)
       do j=1,N
          call fft_conj_pair(j,j1,j2)
          do i=1,N
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
       enddo
    enddo
    !$omp end parallel do
    do k=1,Nby2+1,Nby2
       call fft_conj_pair(k,k1,k2)
       do j=2,Nby2
          call fft_conj_pair(j,j1,j2)
          do i=1,N
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
       enddo
       do j=1,Nby2+1,Nby2
          call fft_conj_pair(j,j1,j2)
          do i=2,Nby2
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
          do i=1,Nby2+1,Nby2 !!self conjugate
             call fft_conj_pair(i,i1,i2)
	     i2=0
             call KSPACE_OPERATION
          enddo
       enddo
    enddo
#undef KSPACE_OPERATION
#else
    stop "k space operation not defined"
#endif

!!If USE_OMP is set to NO:
#else

#ifdef KSPACE_OPERATION
#ifdef REDUCTION_VARS
    !$omp parallel do default(shared) private(i,j,k,i1,j1,k1,i2,j2,k2) reduction(+:REDUCTION_VARS)
#undef REDUCTION_VARS
#else
    !$omp parallel do default(shared) private(i,j,k,i1,j1,k1,i2,j2,k2)
#endif
    do k=2,Nby2
       call fft_conj_pair(k,k1,k2)
       do j=1,N
          call fft_conj_pair(j,j1,j2)
          do i=1,N
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
       enddo
    enddo
    !$omp end parallel do
    do k=1,Nby2+1,Nby2
       call fft_conj_pair(k,k1,k2)
       do j=2,Nby2
          call fft_conj_pair(j,j1,j2)
          do i=1,N
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
       enddo
       do j=1,Nby2+1,Nby2
          call fft_conj_pair(j,j1,j2)
          do i=2,Nby2
             call fft_conj_pair(i,i1,i2)
             call KSPACE_OPERATION
          enddo
          do i=1,Nby2+1,Nby2 !!self conjugate
             call fft_conj_pair(i,i1,i2)
       i2=0
             call KSPACE_OPERATION
          enddo
       enddo
    enddo
#undef KSPACE_OPERATION
#else
    stop "k space operation not defined"
#endif
#endif