module LatticeFFT
  use Utils
  implicit none

  interface FFT3D
     module procedure CubicFFT, CubicMassiveFFT
  end interface
  

#include "configure.h"


  !!constants 
  !!define the FFT direction
  integer(IB),parameter::FFT_FORWARD = -1 
  integer(IB),parameter::FFT_BACKWARD=-FFT_FORWARD
  
  integer(IB)::fft_ind_re=1
  integer(IB)::fft_ind_im=2
  !!work space
  complex(DLC),dimension(n)::FFT_eit
  real(dl),dimension(n)::fft_sin,fft_cos,fft_sinhalf,fft_coshalf,fft_effkdx, fft_effk2dx2
  real(dl),dimension(fft_numk)::fft_vol

  logical::fft_initialized = .false.


contains

  subroutine FFT_Init()
    integer(IB) i,j,k,id
    real(dl)::keff,w
    fft_initialized = .true.
#if USE_OMP
    !$omp parallel do
    do i=1,N
       fft_cos(i)=cos(const_2pi*(i-1)*FFT_FORWARD/n)
       fft_sin(i)=sin(const_2pi*(i-1)*FFT_FORWARD/n)
       fft_eit(i)=cmplx(fft_cos(i),fft_sin(i))
       fft_coshalf(i) = cos(const_pi*(i-1)*FFT_FORWARD/n)
       fft_sinhalf(i) = sin(const_pi*(i-1)*FFT_FORWARD/n)
#if USE_STANDARD_WAVENUMBER
       !! stop "By default this feature is switched off"
       if(i.le.n/2+1)then
          fft_effkdx(i) = const_2pi*(i-1)*FFT_BACKWARD/n
       else
          fft_effkdx(i) = const_2pi*(i-1-n)*FFT_BACKWARD/n
       endif
#elif DIS_SCHEME ==  LATTICEEASY
       if(i.le.n/2+1) then
          fft_effkdx(i) = -2.*fft_sinhalf(i)
       else
          fft_effkdx(i) = 2.*fft_sinhalf(i)
       endif
#elif DIS_SCHEME == HLATTICE1
       fft_effkdx(i) = -fft_sin(i)
#elif DIS_SCHEME == HLATTICE2
       fft_effkdx(i) = -fft_sin(i)/3.*(4.-fft_cos(i))
#else
       stop "Discretization scheme not defined in configure.h"
#endif

       fft_effk2dx2(i) = fft_effkdx(i)**2
    enddo
    !$omp end parallel do
#else
    do i=1,N
       fft_cos(i)=cos(const_2pi*(i-1)*FFT_FORWARD/n)
       fft_sin(i)=sin(const_2pi*(i-1)*FFT_FORWARD/n)
       fft_eit(i)=cmplx(fft_cos(i),fft_sin(i))
       fft_coshalf(i) = cos(const_pi*(i-1)*FFT_FORWARD/n)
       fft_sinhalf(i) = sin(const_pi*(i-1)*FFT_FORWARD/n)
#if USE_STANDARD_WAVENUMBER
       !! stop "By default this feature is switched off"
       if(i.le.n/2+1)then
          fft_effkdx(i) = const_2pi*(i-1)*FFT_BACKWARD/n
       else
          fft_effkdx(i) = const_2pi*(i-1-n)*FFT_BACKWARD/n
       endif
#elif DIS_SCHEME ==  LATTICEEASY
       if(i.le.n/2+1) then
          fft_effkdx(i) = -2.*fft_sinhalf(i)
       else
          fft_effkdx(i) = 2.*fft_sinhalf(i)
       endif
#elif DIS_SCHEME == HLATTICE1
       fft_effkdx(i) = -fft_sin(i)
#elif DIS_SCHEME == HLATTICE2
       fft_effkdx(i) = -fft_sin(i)/3.*(4.-fft_cos(i))
#else
       stop "Discretization scheme not defined in configure.h"
#endif

       fft_effk2dx2(i) = fft_effkdx(i)**2
    enddo
#endif

#if MATCH_CONFIGURATION_SPACE_FLUC
    if(fft_cutoff_ind.lt.n/2) then
       write(*,*) "Since V2.0 a cutoff k is used."
       write(*,*) "with a cut off there is no way to match the configuration space fluctuation to the Fourier space integral (<f^2>=\int |f_k|^2 vol(k) d \ln k)"
       stop
    endif
    fft_vol = 0._dl
#if USE_OMP
    !$omp parallel do reduction(+:fft_vol) private(i,j,k,id,keff,w) default(shared)
    LOOP
    keff=fft_absk_ind(i,j,k)
    id = nint(keff)
    if(id.gt.0 .and. id.le.fft_numk)fft_vol(id) = fft_vol(id) + 1.
    ENDLOOP
    !$omp end parallel do
#else
    LOOP
    keff=fft_absk_ind(i,j,k)
    id = nint(keff)
    if(id.gt.0 .and. id.le.fft_numk)fft_vol(id) = fft_vol(id) + 1.
    ENDLOOP
#endif
#else
    do i=1,fft_numk
       fft_vol(i) = const_pi * 4.* real(i,dl)**2
    enddo
#endif

  end subroutine FFT_Init


  subroutine fft_conj_pair(i,i1,i2)
    integer(IB),intent(in)::i
    integer(IB),intent(out)::i1,i2
    i1=bit_reverse_array(i)
    i2=bit_reverse_array(fft_conj_index(i))
  end subroutine fft_conj_pair


  function fft_nabla2dx2(i,j,k) !!return -k_{eff}^2*(dx^2) 
    integer(IB) i,j,k
    real(dl) fft_nabla2dx2
    fft_nabla2dx2 = - fft_effk2dx2(i) - fft_effk2dx2(j) - fft_effk2dx2(k)
  end function fft_nabla2dx2


  function fft_absk_ind(i,j,k) !! return |k_eff| in unit of 2\pi/L;
    integer(IB) i,j,k
    real(dl) fft_absk_ind
    real(dl),parameter::c = (n/const_2pi)
    if(fft_rawk_index_abs(i).gt.fft_cutoff_ind .or. fft_rawk_index_abs(j) .gt. fft_cutoff_ind .or. fft_rawk_index_abs(k) .gt. fft_cutoff_ind)then
       fft_absk_ind = 0. !!ignore this mode
       return
    endif
    fft_absk_ind = c * sqrt(- fft_nabla2dx2(i,j,k) )
  end function fft_absk_ind

  Subroutine CubicFFT(ref,imf,direction) 
    real(dl),dimension(:,:,:,:),intent(inout)::ref,imf
    integer(IB),intent(in)::direction
    integer(IB) mmax,istep,m,i,j,k,indx,j1,k1
    if(.not. fft_initialized) call fft_init()
    if(size(ref,2).ne.N .or. size(ref,3).ne. n .or. size(ref,4).ne.N) stop "CubicFFT: wrong size of input array"
    if(size(imf,2).ne.N .or. size(imf,3).ne. n .or. size(imf,4).ne.N) stop "CubicFFT: wrong size of input array"
    !!FFT in x and y direction
#if USE_OMP
    !$omp parallel do private(k) default(shared)
    do k=1,N 
       call FFT_2D(ref(:,1:n,1:n,k),imf(:,1:n,1:n,k),direction)
    enddo
    !$omp end parallel do
#else
    do k=1,N 
       call FFT_2D(ref(:,1:n,1:n,k),imf(:,1:n,1:n,k),direction)
    enddo
#endif
    mmax=1
    indx=Nby2
    !!shared memory, output is scrambled format (i,j,k)-> bit_reverse_array(i),bit_reverse_array(j),bit_reverse_array(k)
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
#if USE_OMP
             !$omp parallel do 
             do i=m+1,N,istep
                call FFT_xstep(ref,imf,i,mmax,m,indx)
             enddo
             !$omp end parallel do
#else 
             do i=m+1,N,istep
                call FFT_xstep(ref,imf,i,mmax,m,indx)
             enddo
#endif
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
#if USE_OMP
             !$omp parallel do 
             do i=m+1,N,istep
                call FFT_xstep_backward(ref,imf,i,mmax,m,indx)
             enddo
             !$omp end parallel do
#else
             do i=m+1,N,istep
                call FFT_xstep_backward(ref,imf,i,mmax,m,indx)
             enddo
#endif
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
    if(direction.ne.FFT_FORWARD)then
       ref(fft_ind_re,:,:,:) = ref(fft_ind_re,:,:,:)/ncube !! FFTW does not do this
       imf(fft_ind_im,:,:,:) = imf(fft_ind_im,:,:,:)/ncube 
    endif
  End Subroutine CubicFFT


  Subroutine FFT_RandomGaussian(ref,imf,repk,impk,realspace)
    !! Set up random Gaussian ref and imf in configuration space, if realspace is omitted or set to be true. If realspace is set to be false, the result will be in scrambled momentum space. 
    !!repk is the power spectrum of ref at 2pi/L, 4pi/L, ..., 2*fft_num/L
    !!impk is the power spectrum of imf at 2pi/L, 4pi/L, ..., 2*fft_num/L

    real(dl),intent(in):: repk(fft_numk),impk(fft_numk)
    real(dl),dimension(:,:,:,:)::ref,imf
    real(dl) refk(fft_numk),imfk(fft_numk)
    integer(IB) i,j,k,i1,j1,k1,i2,j2,k2
    logical,optional,intent(in)::realspace
    if(.not. fft_initialized) call fft_init()
    if(any(repk.lt.0.) .or. any(impk.lt.0.)) stop "Error in FFT_randomGaussian :negative power spectrum."
    if(size(ref,2).ne.N .or. size(ref,3).ne. n .or. size(ref,4).ne.N) stop "CubicFFT: wrong size of input array"
    if(size(imf,2).ne.N .or. size(imf,3).ne. n .or. size(imf,4).ne.N) stop "CubicFFT: wrong size of input array"
    if(all(repk.eq.0.d0) .and. all(impk.eq.0.d0))then
       ref(fft_ind_re,:,:,:)= 0._dl
       imf(fft_ind_im,:,:,:)=0._dl
       return
    endif

    refk=sqrt(repk)
    imfk=sqrt(impk)
#define KSPACE_OPERATION assign_grid_value(i,j,k,i1,j1,k1,i2,j2,k2)
#include "headfiles/fft_scan_kspace.h"
    if(present(realspace))then
       if(realspace)then
          call CubicFFT(ref,imf,FFT_BACKWARD)
       endif
    else
       call CubicFFT(ref,imf,FFT_BACKWARD)
    endif

  contains 

    subroutine assign_grid_value(i,j,k,i1,j1,k1,i2,j2,k2)
      integer(IB) loc,i,j,k,i1,j1,k1,i2,j2,k2
      complex(DLC) tmp1,tmp2
      loc=nint(fft_absk_ind(i,j,k))
      if(loc.ge.1)then
         if(i2.gt.0)then
            tmp1 = Gaussian_cmplx()
            tmp2 = Gaussian_cmplx()
         else
            tmp1 = Gaussian_Random()
            tmp2 = cmplx(0._dl, Gaussian_Random())
         endif
         tmp1=tmp1*refk(loc)
         tmp2=tmp2*imfk(loc)
      else
         tmp1=0.
         tmp2=0.
      endif
      ref(fft_ind_re,i1,j1,k1)=real(tmp1)+real(tmp2)
      imf(fft_ind_im,i1,j1,k1)=aimag(tmp1)+aimag(tmp2)
      if(i2.gt.0)then
         ref(fft_ind_re,i2,j2,k2)=real(tmp1)-real(tmp2)
         imf(fft_ind_im,i2,j2,k2)=-aimag(tmp1)+aimag(tmp2)
      endif
    end subroutine assign_grid_value
  End Subroutine FFT_RandomGaussian  

  subroutine FFT_GetPowerSpectrum(ref,imf,repk,impk)
    real(dl),dimension(:,:,:,:)::ref,imf
    real(DL)::repk(fft_numk),impk(fft_numk)
    real(dl) mult(fft_numk),absk,w,re2,im2
    integer(IB) i,j,k,i1,j1,k1,i2,j2,k2
    if(.not. fft_initialized) call fft_init()
    mult=0.
    repk=0.
    impk=0.
#define KSPACE_OPERATION add_power(i,j,k,i1,j1,k1,i2,j2,k2,repk,impk,mult)
#define REDUCTION_VARS repk,impk,mult
#include "headfiles/fft_scan_kspace.h"
    do i=1,fft_numk
       if(mult(i).gt.0.)then
          repk(i)=repk(i)/mult(i)
          impk(i)=impk(i)/mult(i)
       endif
    enddo
  contains
    subroutine add_power(i,j,k,i1,j1,k1,i2,j2,k2,repk,impk,mult)
      integer(IB) loc,i,j,k,i1,j1,k1,i2,j2,k2
      real(dl),dimension(fft_numk)::repk,impk,mult
      loc=nint(fft_absk_ind(i,j,k))
      if(loc.ge.1)then
         if(i2.gt.0)then
            mult(loc)=mult(loc)+4.
            repk(loc)=repk(loc)+(ref(fft_ind_re,i1,j1,k1)+ref(fft_ind_re,i2,j2,k2))**2 + (imf(fft_ind_im,i1,j1,k1)-imf(fft_ind_im,i2,j2,k2))**2
            impk(loc)=impk(loc)+(ref(fft_ind_re,i1,j1,k1)-ref(fft_ind_re,i2,j2,k2))**2 + (imf(fft_ind_im,i1,j1,k1)+imf(fft_ind_im,i2,j2,k2))**2
         else
            mult(loc)=mult(loc)+2.
            repk(loc)=repk(loc)+ref(fft_ind_re,i1,j1,k1)**2*2.
            impk(loc)=impk(loc)+imf(fft_ind_im,i1,j1,k1)**2*2.
         endif
      endif
    end subroutine add_power
  end subroutine FFT_GetPowerSpectrum



  subroutine FFT_1D(ref,imf,direction)
    real(dl),dimension(:,:)::ref,imf
    integer(IB),intent(in)::direction
    integer(IB) mmax,istep,m,i,j,indx
    real(dl) retmp,imtmp
    mmax=1
    indx=Nby2
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(fft_ind_re,bit_reverse_array(i+mmax)) - fft_sin(m*indx+1)*imf(fft_ind_im,bit_reverse_array(i+mmax))
                imtmp=fft_sin(m*indx+1)*ref(fft_ind_re,bit_reverse_array(i+mmax)) + fft_cos(m*indx+1)*imf(fft_ind_im,bit_reverse_array(i+mmax))
                ref(fft_ind_re,bit_reverse_array(i+mmax))=ref(fft_ind_re,bit_reverse_array(i))-retmp
                imf(fft_ind_im,bit_reverse_array(i+mmax))=imf(fft_ind_im,bit_reverse_array(i))-imtmp
                ref(fft_ind_re,bit_reverse_array(i))=ref(fft_ind_re,bit_reverse_array(i))+retmp
                imf(fft_ind_im,bit_reverse_array(i))=imf(fft_ind_im,bit_reverse_array(i))+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(fft_ind_re,i+mmax)+fft_sin(m*indx+1)*imf(fft_ind_im,i+mmax) 
                imtmp=-fft_sin(m*indx+1)*ref(fft_ind_re,i+mmax)+fft_cos(m*indx+1)*imf(fft_ind_im,i+mmax) 
                ref(fft_ind_re,i+mmax)=ref(fft_ind_re,i)-retmp
                imf(fft_ind_im,i+mmax)=imf(fft_ind_im,i)-imtmp
                ref(fft_ind_re,i)=ref(fft_ind_re,i)+retmp
                imf(fft_ind_im,i)=imf(fft_ind_im,i)+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
  end subroutine FFT_1D



  subroutine FFT_2D(ref,imf,direction)
    real(dl),dimension(:,:,:)::ref,imf
    integer(IB),intent(in)::direction
    integer(IB) mmax,istep,m,i,j,indx
    real(dl) retmp(N),imtmp(N)
    do i=1,N
       call FFT_1D(ref(:,:,i),imf(:,:,i),direction)
    enddo
    mmax=1
    indx=Nby2
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(fft_ind_re,1:N,bit_reverse_array(i+mmax))-fft_sin(m*indx+1)*imf(fft_ind_im,1:N,bit_reverse_array(i+mmax))
                imtmp=fft_sin(m*indx+1)*ref(fft_ind_re,1:N,bit_reverse_array(i+mmax))+fft_cos(m*indx+1)*imf(fft_ind_im,1:N,bit_reverse_array(i+mmax))
                ref(fft_ind_re,1:N,bit_reverse_array(i+mmax))=ref(fft_ind_re,1:N,bit_reverse_array(i))-retmp
                imf(fft_ind_im,1:N,bit_reverse_array(i+mmax))=imf(fft_ind_im,1:N,bit_reverse_array(i))-imtmp
                ref(fft_ind_re,1:N,bit_reverse_array(i))=ref(fft_ind_re,1:N,bit_reverse_array(i))+retmp
                imf(fft_ind_im,1:N,bit_reverse_array(i))=imf(fft_ind_im,1:N,bit_reverse_array(i))+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp = fft_cos(m*indx+1)*ref(fft_ind_re,1:N,i+mmax) + fft_sin(m*indx+1)*imf(fft_ind_im,1:N,i+mmax)
                imtmp = -fft_sin(m*indx+1)*ref(fft_ind_re,1:N,i+mmax) + fft_cos(m*indx+1)*imf(fft_ind_im,1:N,i+mmax)
                ref(fft_ind_re,1:N,i+mmax)=ref(fft_ind_re,1:N,i)-retmp
                imf(fft_ind_im,1:N,i+mmax)=imf(fft_ind_im,1:N,i)-imtmp
                ref(fft_ind_re,1:N,i)=ref(fft_ind_re,1:N,i)+retmp
                imf(fft_ind_im,1:N,i)=imf(fft_ind_im,1:N,i)+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
  end subroutine FFT_2D


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine FFT_xstep(ref,imf,i,mmax,m,indx)
    integer(IB) i,mmax,m,indx
    real(dl),dimension(:,:,:,:)::ref,imf
    real(dl):: retmp(N,N),imtmp(N,N)
    retmp=fft_cos(m*indx+1)*ref(fft_ind_re,1:N,1:N,bit_reverse_array(i+mmax))-fft_sin(m*indx+1)*imf(fft_ind_im,1:N,1:N,bit_reverse_array(i+mmax))
    imtmp=fft_sin(m*indx+1)*ref(fft_ind_re,1:N,1:N,bit_reverse_array(i+mmax))+fft_cos(m*indx+1)*imf(fft_ind_im,1:N,1:N,bit_reverse_array(i+mmax))
    ref(fft_ind_re,1:N,1:N,bit_reverse_array(i+mmax))=ref(fft_ind_re,1:N,1:N,bit_reverse_array(i))-retmp
    imf(fft_ind_im,1:N,1:N,bit_reverse_array(i+mmax))=imf(fft_ind_im,1:N,1:N,bit_reverse_array(i))-imtmp
    ref(fft_ind_re,1:N,1:N,bit_reverse_array(i))=ref(fft_ind_re,1:N,1:N,bit_reverse_array(i))+retmp
    imf(fft_ind_im,1:N,1:N,bit_reverse_array(i))=imf(fft_ind_im,1:N,1:N,bit_reverse_array(i))+imtmp
  end subroutine FFT_xstep


  subroutine FFT_xstep_backward(ref,imf,i,mmax,m,indx)
    integer(IB) i,mmax,m,indx
    real(dl),dimension(:,:,:,:)::ref,imf
    real(dl):: retmp(N,N),imtmp(N,N)
    retmp=fft_cos(m*indx+1)* ref(fft_ind_re,1:N,1:N,i+mmax)+fft_sin(m*indx+1)* imf(fft_ind_im,1:N,1:N,i+mmax)
    imtmp=-fft_sin(m*indx+1)* ref(fft_ind_re,1:N,1:N,i+mmax)+fft_cos(m*indx+1)* imf(fft_ind_im,1:N,1:N,i+mmax)
    ref(fft_ind_re,1:N,1:N,i+mmax)=ref(fft_ind_re,1:N,1:N,i)-retmp
    imf(fft_ind_im,1:N,1:N,i+mmax)=imf(fft_ind_im,1:N,1:N,i)-imtmp
    ref(fft_ind_re,1:N,1:N,i)=ref(fft_ind_re,1:N,1:N,i)+retmp
    imf(fft_ind_im,1:N,1:N,i)=imf(fft_ind_im,1:N,1:N,i)+imtmp
  end subroutine FFT_xstep_backward

!!******************** massive FFT *****************************

  Subroutine CubicMassiveFFT(msize,ref,imf,direction) 
    real(dl),dimension(:,:,:,:),intent(inout)::ref,imf
    integer(IB),intent(in)::direction,msize
    integer(IB) mmax,istep,m,i,j,k,indx,j1,k1
    if(.not. fft_initialized) call fft_init()
    if(size(ref,2).ne.N .or. size(ref,3).ne. n .or. size(ref,4).ne.N .or. size(ref,1) .ne. msize .or. size(imf,1).ne.msize .or. size(imf,2).ne.N .or. size(imf,3).ne. n .or. size(imf,4).ne.N ) stop "CubicMassiveFFT: wrong size of input array"
    !!FFT in x and y direction

#if USE_OMP
    !$omp parallel do private(k) default(shared)
    do k=1,N 
       call mfft_2D(msize,ref(:,1:N,1:N,k),imf(:,1:n,1:n,k),direction)
    enddo
    !$omp end parallel do
#else
    do k=1,N 
       call mfft_2D(msize,ref(:,1:N,1:N,k),imf(:,1:n,1:n,k),direction)
    enddo
#endif
    mmax=1
    indx=Nby2
    !!shared memory, output is scrambled format (i,j,k)-> bit_reverse_array(i),bit_reverse_array(j),bit_reverse_array(k)
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
#if USE_OMP
             !$omp parallel do 
             do i=m+1,N,istep
                call mfft_xstep(msize,ref,imf,i,mmax,m,indx)
             enddo
             !$omp end parallel do
#else
             do i=m+1,N,istep
                call mfft_xstep(msize,ref,imf,i,mmax,m,indx)
             enddo
#endif
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
#if USE_OMP
             !$omp parallel do 
             do i=m+1,N,istep
                call mfft_xstep_backward(msize,ref,imf,i,mmax,m,indx)
             enddo
             !$omp end parallel do
#else
             do i=m+1,N,istep
                call mfft_xstep_backward(msize,ref,imf,i,mmax,m,indx)
             enddo
#endif
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
    if(direction.ne.FFT_FORWARD)then
       ref(:,:,:,:) = ref(:,:,:,:)/ncube !! FFTW does not do this
       imf(:,:,:,:) = imf(:,:,:,:)/ncube 
    endif
  End Subroutine CubicMassiveFFT


  subroutine mFFT_1D(msize,ref,imf,direction)
    integer(IB) msize
    real(dl),dimension(:,:)::ref,imf
    integer(IB),intent(in)::direction
    integer(IB) mmax,istep,m,i,j,indx
    real(dl) retmp(msize),imtmp(msize)
    mmax=1
    indx=Nby2
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(:,bit_reverse_array(i+mmax)) - fft_sin(m*indx+1)*imf(:,bit_reverse_array(i+mmax))
                imtmp=fft_sin(m*indx+1)*ref(:,bit_reverse_array(i+mmax)) + fft_cos(m*indx+1)*imf(:,bit_reverse_array(i+mmax))
                ref(:,bit_reverse_array(i+mmax))=ref(:,bit_reverse_array(i))-retmp
                imf(:,bit_reverse_array(i+mmax))=imf(:,bit_reverse_array(i))-imtmp
                ref(:,bit_reverse_array(i))=ref(:,bit_reverse_array(i))+retmp
                imf(:,bit_reverse_array(i))=imf(:,bit_reverse_array(i))+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(:,i+mmax)+fft_sin(m*indx+1)*imf(:,i+mmax) 
                imtmp=-fft_sin(m*indx+1)*ref(:,i+mmax)+fft_cos(m*indx+1)*imf(:,i+mmax) 
                ref(:,i+mmax)=ref(:,i)-retmp
                imf(:,i+mmax)=imf(:,i)-imtmp
                ref(:,i)=ref(:,i)+retmp
                imf(:,i)=imf(:,i)+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
  end subroutine MFFT_1D



  subroutine mFFT_2D(msize,ref,imf,direction)
    integer(IB):: msize
    real(dl),dimension(:,:,:)::ref,imf
    integer(IB),intent(in)::direction
    integer(IB) mmax,istep,m,i,j,indx
    real(dl) retmp(msize,N),imtmp(msize,N)
    do i=1,N
       call mFFT_1D(msize,ref(:,:,i),imf(:,:,i),direction)
    enddo

    mmax=1
    indx=Nby2
    if(direction.eq.FFT_FORWARD)then
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp=fft_cos(m*indx+1)*ref(:,1:N,bit_reverse_array(i+mmax))-fft_sin(m*indx+1)*imf(:,1:N,bit_reverse_array(i+mmax))
                imtmp=fft_sin(m*indx+1)*ref(:,1:N,bit_reverse_array(i+mmax))+fft_cos(m*indx+1)*imf(:,1:N,bit_reverse_array(i+mmax))
                ref(:,1:N,bit_reverse_array(i+mmax))=ref(:,1:N,bit_reverse_array(i))-retmp
                imf(:,1:N,bit_reverse_array(i+mmax))=imf(:,1:N,bit_reverse_array(i))-imtmp
                ref(:,1:N,bit_reverse_array(i))=ref(:,1:N,bit_reverse_array(i))+retmp
                imf(:,1:N,bit_reverse_array(i))=imf(:,1:N,bit_reverse_array(i))+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    else
       do while(N.gt.mmax)
          istep=2*mmax
          do m=0,mmax-1
             do i=m+1,N,istep
                retmp = fft_cos(m*indx+1)*ref(:,1:N,i+mmax) + fft_sin(m*indx+1)*imf(:,1:N,i+mmax)
                imtmp = -fft_sin(m*indx+1)*ref(:,1:N,i+mmax) + fft_cos(m*indx+1)*imf(:,1:N,i+mmax)
                ref(:,1:N,i+mmax)=ref(:,1:N,i)-retmp
                imf(:,1:N,i+mmax)=imf(:,1:N,i)-imtmp
                ref(:,1:N,i)=ref(:,1:N,i)+retmp
                imf(:,1:N,i)=imf(:,1:N,i)+imtmp
             enddo
          enddo
          mmax=istep
          indx=indx/2
       enddo
    endif
  end subroutine MFFT_2D


  !!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  subroutine mfft_xstep(msize,ref,imf,i,mmax,m,indx)
    integer(IB) msize
    integer(IB) i,mmax,m,indx
    real(dl),dimension(:,:,:,:)::ref,imf
    real(dl):: retmp(msize,N,N),imtmp(msize,N,N)
    retmp=fft_cos(m*indx+1)*ref(:,1:N,1:N,bit_reverse_array(i+mmax))-fft_sin(m*indx+1)*imf(:,1:N,1:N,bit_reverse_array(i+mmax))
    imtmp=fft_sin(m*indx+1)*ref(:,1:N,1:N,bit_reverse_array(i+mmax))+fft_cos(m*indx+1)*imf(:,1:N,1:N,bit_reverse_array(i+mmax))
    ref(:,1:N,1:N,bit_reverse_array(i+mmax))=ref(:,1:N,1:N,bit_reverse_array(i))-retmp
    imf(:,1:N,1:N,bit_reverse_array(i+mmax))=imf(:,1:N,1:N,bit_reverse_array(i))-imtmp
    ref(:,1:N,1:N,bit_reverse_array(i))=ref(:,1:N,1:N,bit_reverse_array(i))+retmp
    imf(:,1:N,1:N,bit_reverse_array(i))=imf(:,1:N,1:N,bit_reverse_array(i))+imtmp
  end subroutine mfft_xstep


  subroutine mfft_xstep_backward(msize,ref,imf,i,mmax,m,indx)
    integer(IB) msize
    integer(IB) i,mmax,m,indx
    real(dl),dimension(:,:,:,:)::ref,imf
    real(dl):: retmp(msize,N,N),imtmp(msize,N,N)
    retmp=fft_cos(m*indx+1)* ref(:,1:N,1:N,i+mmax)+fft_sin(m*indx+1)* imf(:,1:N,1:N,i+mmax)
    imtmp=-fft_sin(m*indx+1)* ref(:,1:N,1:N,i+mmax)+fft_cos(m*indx+1)* imf(:,1:N,1:N,i+mmax)
    ref(:,1:N,1:N,i+mmax)=ref(:,1:N,1:N,i)-retmp
    imf(:,1:N,1:N,i+mmax)=imf(:,1:N,1:N,i)-imtmp
    ref(:,1:N,1:N,i)=ref(:,1:N,1:N,i)+retmp
    imf(:,1:N,1:N,i)=imf(:,1:N,1:N,i)+imtmp
  end subroutine Mfft_xstep_backward



!!the input h and p should be in Fourier space,
!!this subroutine removes the scalar and vector component in h and p
!!and it gets the power spectrum of TT part in h and p
!!output_option can be: UNCHANGED (keep h and p unchanged), OUTPUT_TT (remove the scalar and vector part in h and p), OUTPUT_ADAPTIVE (perform a coordinate transformation x^i -> x^i - zeta^i; change h, and write zeta^i into p; In this case the original data stored in p is lost. You have to reload p from checkpoint file).
  subroutine fft_get_TT_power(h,p,hpk,ppk,output_option)
    integer(IB)::output_option
    real(dl) h(6,n,n,n),p(6,n,n,n),hpk(fft_numk),ppk(fft_numk),mult(fft_numk)
    integer(IB) i,j,k,i1,j1,k1,i2,j2,k2
    hpk = 0._dl
    ppk = 0._dl
    mult = 0._dl
#define KSPACE_OPERATION add_TT_power(i,j,k,i1,j1,k1,i2,j2,k2,hpk,ppk,mult)
#define REDUCTION_VARS hpk,ppk,mult
#include "headfiles/fft_scan_kspace.h"
    do i=1,fft_numk
       if(mult(i).gt.0)then
          hpk(i) = hpk(i)/mult(i)
          ppk(i) = ppk(i)/mult(i)
       endif
    enddo
    contains
      subroutine add_TT_power(i,j,k,i1,j1,k1,i2,j2,k2,hpk,ppk,mult)
        integer(IB) i,j,k,i1,j1,k1,i2,j2,k2,loc
        real(dl) hpk(fft_numk),ppk(fft_numk),mult(fft_numk)
        complex(dlc) two_hk(6), two_pk(6),lambda(2),A(3),Ap(3),trh1(3),trh2(3)
        real(dl) ksq, didj(6)
        ksq = fft_effk2dx2(i)+fft_effk2dx2(j)+fft_effk2dx2(k)
        if(ksq.le.(1.e-8/ncube))then
           if(output_option .eq. OUTPUT_TT)then
              h(:,i1,j1,k1) = 0.
              p(:,i1,j1,k1)=0.
              if(i2.gt.0)then
                 h(:,i2,j2,k2) = 0.
                 p(:,i2,j2,k2)=0.
              endif
           elseif(output_option .eq. OUTPUT_ADAPTIVE)then
              p(:,i1,j1,k1)=0.
              if(i2.gt.0) p(:,i2,j2,k2)=0.
           endif
           return
        endif
        loc = nint(fft_absk_ind(i,j,k))
        if(loc.le.0 .and. output_option.eq. UNCHANGED) return
        if(i2.gt.0)then
           two_hk = cmplx(h(:,i1,j1,k1)+h(:,i2,j2,k2),p(:,i1,j1,k1)-p(:,i2,j2,k2))
           two_pk = cmplx(p(:,i1,j1,k1)+p(:,i2,j2,k2), -h(:,i1,j1,k1)+h(:,i2,j2,k2))
           if(loc.gt.0)mult(loc) = mult(loc)+4.
        else
           two_hk = h(:,i1,j1,k1)*const_sqrt2
           two_pk = p(:,i1,j1,k1)*const_sqrt2 !!here I am actually calculating sqrt(2)*h_k and sqrt(2)*p_k, the extra factor sqrt(2)**2 will be cancelled by the weight
           if(loc.gt.0)mult(loc) = mult(loc) + 2.
        endif
        !!remove the trace
        trh1 = sum(two_hk(1:3))/3.
        two_hk(1:3) = two_hk(1:3) - trh1 
        two_pk(1:3) = two_pk(1:3) - sum(two_pk(1:3))/3.
        !!remove the other scalar part 
        didj= (/ fft_effk2dx2(i), fft_effk2dx2(j), fft_effk2dx2(k), &
             fft_effkdx(j)*fft_effkdx(k), fft_effkdx(i)*fft_effkdx(k), fft_effkdx(i)* fft_effkdx(j) /) 
        lambda = 1.5/ksq**2 * (/ sum(didj(1:3)*two_hk(1:3))+2.*sum(didj(4:6)*two_hk(4:6)), sum(didj(1:3)*two_pk(1:3))+2.*sum(didj(4:6)*two_pk(4:6)) /)
        two_hk = two_hk - lambda(1) * didj
        two_pk = two_pk - lambda(2) * didj
        trh2 =  sum(two_hk(1:3))/3.
        two_hk(1:3) = two_hk(1:3) - trh2
        two_pk(1:3) = two_pk(1:3) - sum(two_pk(1:3))/3.

        !!remove the vector part
        A(1) = (fft_effkdx(i) * two_hk(1) + fft_effkdx(j) * two_hk(6) + fft_effkdx(k) * two_hk(5))/ksq
        A(2) = (fft_effkdx(i) * two_hk(6) + fft_effkdx(j) * two_hk(2) + fft_effkdx(k) * two_hk(4))/ksq
        A(3) = (fft_effkdx(i) * two_hk(5) + fft_effkdx(j) * two_hk(4) + fft_effkdx(k) * two_hk(3))/ksq
        two_hk = two_hk - (/ 2.*fft_effkdx(i) * A(1), 2.*fft_effkdx(j) * A(2), 2.*fft_effkdx(k) *A(3), &
             fft_effkdx(j)*A(3) + fft_effkdx(k)*A(2), fft_effkdx(i) * A(3)+ fft_effkdx(k)*A(1) , &
             fft_effkdx(i)*A(2) + fft_effkdx(j)*A(1) /)

        Ap(1) = (fft_effkdx(i) * two_pk(1) + fft_effkdx(j) * two_pk(6) + fft_effkdx(k) * two_pk(5))/ksq
        Ap(2) = (fft_effkdx(i) * two_pk(6) + fft_effkdx(j) * two_pk(2) + fft_effkdx(k) * two_pk(4))/ksq
        Ap(3) = (fft_effkdx(i) * two_pk(5) + fft_effkdx(j) * two_pk(4) + fft_effkdx(k) * two_pk(3))/ksq
        two_pk = two_pk - (/ 2.*fft_effkdx(i) * Ap(1), 2.*fft_effkdx(j) * Ap(2), 2.*fft_effkdx(k) *Ap(3), &
             fft_effkdx(j)*Ap(3) + fft_effkdx(k)*Ap(2), fft_effkdx(i) * Ap(3)+fft_effkdx(k)*Ap(1) , &
             fft_effkdx(i)*Ap(2) + fft_effkdx(j)*Ap(1) /)

        !! now add the TT power
        if(loc.gt.0)then
           hpk(loc) = hpk(loc) + sum(two_hk(1:3)*conjg(two_hk(1:3))) + 2.* sum(two_hk(4:6)*conjg(two_hk(4:6)))
           ppk(loc) = ppk(loc) + sum(two_pk(1:3)*conjg(two_pk(1:3))) + 2.* sum(two_pk(4:6)*conjg(two_pk(4:6)))
        endif
        select case(output_option)
        case(UNCHANGED)
           return
        case(OUTPUT_TT)
           if(i2.gt.0)then
              h(:,i1,j1,k1) = (real(two_hk)-aimag(two_pk))/2.
              p(:,i1,j1,k1) = (aimag(two_hk)+real(two_pk))/2.
              h(:,i2,j2,k2) = (real(two_hk)+aimag(two_pk))/2.
              p(:,i2,j2,k2) = (-aimag(two_hk)+real(two_pk))/2.
           else
              h(:,i1,j1,k1) = (real(two_hk)-aimag(two_pk))/const_sqrt2
              p(:,i1,j1,k1) = (aimag(two_hk)+real(two_pk))/const_sqrt2
           endif
        case(OUTPUT_ADAPTIVE)
           two_hk(1:3) = two_hk(1:3)+ trh1 + trh2 !!restore the trace (physical mode)
           two_pk (1:3) = (A + (/ fft_effkdx(i), fft_effkdx(j), fft_effkdx(k) /) * lambda(1))/cmplx(0._dl,1._dl)
           if(i2.gt.0)then
              h(:,i1,j1,k1) = (real(two_hk)-aimag(two_pk))/2.
              p(:,i1,j1,k1) = (aimag(two_hk)+real(two_pk))/2.
              h(:,i2,j2,k2) = (real(two_hk)+aimag(two_pk))/2.
              p(:,i2,j2,k2) = (-aimag(two_hk)+real(two_pk))/2.
           else
              h(:,i1,j1,k1) = (real(two_hk)-aimag(two_pk))/const_sqrt2
              p(:,i1,j1,k1) = (aimag(two_hk)+real(two_pk))/const_sqrt2
           endif
        case default
           stop "unknown output option in fft_remove_scalar_vector"
        end select
      end subroutine add_TT_power
    end subroutine fft_get_TT_power


end module LatticeFFT




