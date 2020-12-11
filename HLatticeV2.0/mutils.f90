!!this module provide functions to calculate the potential/kinetic/gradient energy of the fields/gravity.
module mutils
  use model
  use latticefft
  implicit none
#include "configure.h"
  
  logical::need_recalculate=.true.
  real(dl)::last_potential_energy, last_fields_kinetic_energy, last_fields_gradient_energy, last_gravity_kinetic_energy, last_gravity_gradient_energy, last_total_detg, last_effective_Hubble
  real(dl)::cur_ainc,cur_total_detg
contains

  !!the total potential energy of the scalar fields
  function potential_energy()
    real(dl) potential_energy,cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       potential_energy=last_potential_energy
       return
    endif
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + potential(fields_f(:,i,j,k))*DETG(i,j,k)
    ENDLOOP      
    !$omp end parallel do
    potential_energy = sum(cach)/ncube
    return
  end function potential_energy

  !!the total kinetic energy of the scalar fields
  function fields_kinetic_energy()
    real(dl) fields_kinetic_energy,cachj(n),cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       fields_kinetic_energy = last_fields_kinetic_energy
       return
    endif
#if METRIC_PERTURB
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + sum(fields_p(:,i,j,k)**2)/DETG(i,j,k)
    ENDLOOP
    !$omp end parallel do
    fields_kinetic_energy = sum(cach) / 2._dl /ncube
#else
    fields_kinetic_energy = sum(fields_p**2)/2._dl/ncube
#if METRIC_OPTION == FRW_BACKGROUND
    fields_kinetic_energy = fields_kinetic_energy/scale_factor()**6
#endif
#endif
    return
  end function fields_kinetic_energy

  !!the total gradient energy of the scalar fields
  function fields_gradient_energy()
    real(dl) fields_gradient_energy,cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       fields_gradient_energy = last_fields_gradient_energy
       return
    endif
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) +   &
#if DIS_SCHEME == LATTICEEASY
         sum((GRID_FLD(fields_f,:,i+1,j,k)-GRID_FLD(fields_f,:,i,j,k))**2 &
         +(GRID_FLD(fields_f,:,i,j+1,k)-GRID_FLD(fields_f,:,i,j,k))**2 &
         +(GRID_FLD(fields_f,:,i,j,k+1)-GRID_FLD(fields_f,:,i,j,k))**2)
#elif DIS_SCHEME ==  HLATTICE1
#if METRIC_PERTURB
         DETG(i,j,k) * ( &
         GUP11(i,j,k) * sum(TWO_DFLD_X(i,j,k)**2) &
         + GUP22(i,j,k) * sum(TWO_DFLD_Y(i,j,k)**2) &
         + GUP33(i,j,k) * sum(TWO_DFLD_Z(i,j,k)**2) &
         + 2._dl*( &
         GUP23(i,j,k) * sum(TWO_DFLD_Y(i,j,k)*TWO_DFLD_Z(i,j,k)) &
         + GUP13(i,j,k) * sum(TWO_DFLD_X(i,j,k)*TWO_DFLD_Z(i,j,k)) &
         + GUP12(i,j,k) * sum(TWO_DFLD_X(i,j,k)*TWO_DFLD_Y(i,j,k)) &
         ))
#else
         sum(TWO_DFLD_X(i,j,k)**2) &
         + sum(TWO_DFLD_Y(i,j,k)**2) &
         + sum(TWO_DFLD_Z(i,j,k)**2)
#endif
#elif DIS_SCHEME == HLATTICE2
#if METRIC_PERTURB
         DETG(i,j,k) * ( &
         GUP11(i,j,k) * sum(TWELVE_DFLD_X(i,j,k)**2) &
         + GUP22(i,j,k) * sum(TWELVE_DFLD_Y(i,j,k)**2) &
         + GUP33(i,j,k) * sum(TWELVE_DFLD_Z(i,j,k)**2) &
         + 2._dl*( &
         GUP23(i,j,k) * sum(TWELVE_DFLD_Y(i,j,k)*TWELVE_DFLD_Z(i,j,k)) &
         + GUP13(i,j,k) * sum(TWELVE_DFLD_X(i,j,k)*TWELVE_DFLD_Z(i,j,k)) &
         + GUP12(i,j,k) * sum(TWELVE_DFLD_X(i,j,k)*TWELVE_DFLD_Y(i,j,k)) &
         ))
#else
         sum(TWELVE_DFLD_X(i,j,k)**2) &
         + sum(TWELVE_DFLD_Y(i,j,k)**2) &
         + sum(TWELVE_DFLD_Z(i,j,k)**2)
#endif
#endif
    ENDLOOP
    !$omp end parallel do
#if DIS_SCHEME == LATTICEEASY
    fields_gradient_energy = sum(cach) /ncube/ (metric%physdx)**2 / 2._dl
#elif DIS_SCHEME == HLATTICE1
    fields_gradient_energy = sum(cach) /ncube/ (metric%physdx)**2 / 8._dl
#elif DIS_SCHEME == HLATTICE2
    fields_gradient_energy = sum(cach) /ncube/ (metric%physdx)**2 / 288._dl
#endif
  end function fields_gradient_energy


!!For metric
  function gravity_kinetic_energy()
    real(dl) gravity_kinetic_energy,cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       gravity_kinetic_energy = last_gravity_kinetic_energy
       return
    endif
#if METRIC_PERTURB
    cach=0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + ( sum(MEPIV(i,j,k)**2) + 2._dl*sum(MEPIU(i,j,k)**2) - sum(MEPIU(i,j,k))**2 )/DETG(i,j,k)
    ENDLOOP
    !$omp end parallel do
    gravity_kinetic_energy = sum(cach)/Mplsq/ncube
#else
    gravity_kinetic_energy = -3.*effective_Hubble()**2*Mplsq
#endif
    return
  end function gravity_kinetic_energy


  function gravity_gradient_energy(dg)
    !!only keep O(h^2) terms, but I keep the g^{ij}'s to keep the formula covariant. Keeping all terms covariant is important when we do a coordinate transformation (no leaking of energy).
    real(dl) gravity_gradient_energy,cach(n)
    real(dl),optional::dg
    DEFINE_IND
#if ! METRIC_PERTURB
    gravity_gradient_energy = 0.
    return
#else
    if(.not.need_recalculate)then
       gravity_gradient_energy=last_gravity_gradient_energy
       return
    endif
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + ( &
#if DIS_SCHEME == HLATTICE1
         TWO_DV_X(1,i,j,k)**2 + TWO_DV_Y(2,i,j,k)**2 + TWO_DV_Z(3,i,j,k)**2  &
         - TWO_DU_Z(1,i,j,k)*TWO_DU_Z(2,i,j,k)- TWO_DU_X(2,i,j,k)*TWO_DU_X(3,i,j,k) - TWO_DU_Y(1,i,j,k)*TWO_DU_Y(3,i,j,k)  &
         )/2. &
         - (TWO_DV_X(1,i,j,k)+TWO_DV_Z(3,i,j,k))*TWO_DV_Y(2,i,j,k) - TWO_DV_X(1,i,j,k)*TWO_DV_Z(3,i,j,k) &
         + TWO_DV_Y(1,i,j,k)*TWO_DU_Z(1,i,j,k) + TWO_DV_Z(2,i,j,k)*TWO_DU_X(2,i,j,k) + TWO_DV_X(3,i,j,k)*TWO_DU_Y(3,i,j,k)
#elif DIS_SCHEME == HLATTICE2
         TWELVE_DV_X(1,i,j,k)**2 + TWELVE_DV_Y(2,i,j,k)**2 + TWELVE_DV_Z(3,i,j,k)**2  &
         - TWELVE_DU_Z(1,i,j,k)*TWELVE_DU_Z(2,i,j,k)- TWELVE_DU_X(2,i,j,k)*TWELVE_DU_X(3,i,j,k) - TWELVE_DU_Y(1,i,j,k)*TWELVE_DU_Y(3,i,j,k)  &
         )/2. &
         - (TWELVE_DV_X(1,i,j,k)+TWELVE_DV_Z(3,i,j,k))*TWELVE_DV_Y(2,i,j,k) - TWELVE_DV_X(1,i,j,k)*TWELVE_DV_Z(3,i,j,k) &
         + TWELVE_DV_Y(1,i,j,k)*TWELVE_DU_Z(1,i,j,k) + TWELVE_DV_Z(2,i,j,k)*TWELVE_DU_X(2,i,j,k) + TWELVE_DV_X(3,i,j,k)*TWELVE_DU_Y(3,i,j,k)
#endif
    ENDLOOP
    !$omp end parallel do

    if(present(dg))then
       gravity_gradient_energy = sum(cach)/ncube * dg**(1./3._dl)*Mplsq/metric%physdx**2 
    else
       gravity_gradient_energy = sum(cach)/ncube * total_detg()**(1./3._dl)*Mplsq/metric%physdx**2 
    endif
#if DIS_SCHEME == HLATTICE1
    gravity_gradient_energy = gravity_gradient_energy/8._dl
#elif DIS_SCHEME == HLATTICE2
    gravity_gradient_energy = gravity_gradient_energy/288._dl
#endif
    return
#endif
  end function gravity_gradient_energy

  function total_energy()
    real(dl) total_energy
    total_energy=potential_energy()+ fields_gradient_energy() + fields_kinetic_energy() + gravity_kinetic_energy() + gravity_gradient_energy()
  end function total_energy

  function total_detg()
    real(dl) total_detg,cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       total_detg=last_total_detg
       return
    endif
#if METRIC_PERTURB
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + DETG(i,j,k)
    ENDLOOP
    !$omp end parallel do
    total_detg = sum(cach)/ncube
#else
    total_detg = 1._dl
#endif
    return
  end function total_detg


   function scale_factor()
     real(dl) scale_factor
#if METRIC_PERTURB
     scale_factor = (total_detg())**(1._dl/3._dl)
#else
#if USE_CONFORMAL_TIME
     scale_factor = metric%y
#else
     scale_factor = metric%y**(2./3._dl)
#endif
#endif     
   end function scale_factor


  function effective_Hubble(totdg)
    real(dl) effective_Hubble
    real(dl),optional::totdg
    real(dl) cach(n)
    DEFINE_IND
    if(.not. need_recalculate)then
       effective_Hubble = last_effective_Hubble
       return
    endif
#if METRIC_PERTURB
    cach = 0._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    cach(k) = cach(k) + sum(MEPIU(i,j,k))*DETG(i,j,k)
    ENDLOOP
    !$omp end parallel do
    if(present(totdg))then
       effective_Hubble = -sum(cach)/ncube/totdg/3._dl/Mplsq
    else
       effective_Hubble =  -sum(cach)/ncube/total_detg()/3._dl/Mplsq
    endif
#else
#if USE_CONFORMAL_TIME
    effective_Hubble = -(metric%piy/metric%y**2/6._dl/Mplsq)
#else
    effective_Hubble = -(metric%piy/metric%y/4._dl/Mplsq)
#endif
#endif
  end function effective_Hubble

  subroutine coor_trans()
    real(dl) ainc
    !!================================================
    !!coordinate transformation: redefine dx-> a dx
#if METRIC_PERTURB 
    ainc = (total_detg())**(1._dl/3._dl)
    cur_ainc = ainc
    cur_total_detg = total_detg()
    write(*,*) "Ainc:",cur_ainc, "Total Detg:",cur_total_detg
    metric%a = metric%a*ainc
    metric%physdx = metric%physdx * ainc
    metric_h(1:3,:,:,:) = metric_h(1:3,:,:,:) - 2._dl*log(ainc)
    ainc = ainc**3
    metric_p = metric_p/ainc
    fields_p = fields_p/ainc
#else
    metric%a = scale_factor()
    metric%physdx = metric%dx*metric%a
#endif
    need_recalculate = .true.
    !!================================================
  end subroutine coor_trans


  subroutine get_dgup(i, j , k, dg, dgup, ddgup_dh)
    integer(IB) i,j,k
    real(dl) dg, dgup(6),ddgup_dh(6,6)
    real(dl) eh6,eh3,seii(3),eii(3), eii3(3), eii6(3), eii23(3), ss(3)
#if METRIC_PERTURB
    eii6(1) = fast_exp(-H11(i,j,k)/6._dl)
    eii6(2) = fast_exp(-H22(i,j,k)/6._dl)
    eii6(3) = fast_exp(-H33(i,j,k)/6._dl)

    eh6 = 1./(eii6(1)*eii6(2)*eii6(3))
    eh3 = eh6**2
    dg = eh6*eh3  !!sqrt(|det g|) = e^{h/2}

    eii3 = eii6**2
    eii23 = eii3 **2
    seii = eii6*eii3
    eii = seii**2

    ss = eii23*(/ H12(i,j,k)**2*eii3(2) + H13(i,j,k)**2*eii3(3), &
         H21(i,j,k)**2*eii3(1) + H23(i,j,k)**2*eii3(3), &
         H31(i,j,k)**2*eii3(1) + H32(i,j,k)**2*eii3(2) /)
    dgup(1:3) = ( eii + ss/2.)*dg
    dgup(4:6) = -MEV(i,j,k)/seii + 0.5_dl*eh6* (/ H12(i,j,k)*H13(i,j,k), &
         H21(i,j,k)*H23(i,j,k) , H31(i,j,k)*H32(i,j,k) /)

    ddgup_dh(1:3,4) = (/ 0._dl, eii3(3), eii3(2) /) * eii23 * (dg*H23(i,j,k))
    ddgup_dh(1:3,5) = (/ eii3(3), 0._dl, eii3(1) /) * eii23 * (dg*H13(i,j,k))
    ddgup_dh(1:3,6) = (/ eii3(2), eii3(1), 0._dl /) * eii23 * (dg*H12(i,j,k))
    ddgup_dh(4:6,4) = (/ -1._dl/seii(1), 0._dl, 0._dl /) &
         + 0.5_dl*eh6* (/ 0._dl, H12(i,j,k), H13(i,j,k) /)
    ddgup_dh(4:6,5) = (/ 0._dl, -1._dl/seii(2), 0._dl /) &
         + 0.5_dl*eh6* (/ H12(i,j,k), 0._dl, H23(i,j,k) /)
    ddgup_dh(4:6,6) = (/ 0._dl, 0._dl, -1._dl/seii(3) /) &
         + 0.5_dl*eh6* (/ H13(i,j,k), H23(i,j,k), 0._dl /)

    ddgup_dh(1:3,1) =  dgup(1:3)/2._dl - dg *(/ &
          eii(1)+ss(1)/3._dl, H21(i,j,k)**2*eii3(1)*eii23(2)/6._dl, H31(i,j,k)**2*eii3(1)*eii23(3)/6._dl /)
    ddgup_dh(1:3,2) =  dgup(1:3)/2._dl - dg * (/ &
         H21(i,j,k)**2 * eii3(2)*eii23(1)/6._dl, eii(2)+ss(2)/3._dl, H32(i,j,k)**2*eii3(2)*eii23(3) /)
    ddgup_dh(1:3,3) =  dgup(1:3)/2._dl - dg *(/ &
         H31(i,j,k)**2*eii3(3)*eii23(1), H32(i,j,k)**2*eii3(3)*eii23(2), eii(3)+ss(3)/3._dl /)
    ddgup_dh(4:6,1) =  eh6/12._dl*  (/ H12(i,j,k)*H13(i,j,k), H21(i,j,k)*H23(i,j,k) , H31(i,j,k)*H32(i,j,k) /)
    ddgup_dh(4:6,2) = ddgup_dh(4:6,1)
    ddgup_dh(4:6,3) = ddgup_dh(4:6,1)
    ddgup_dh(4,1) = ddgup_dh(4,1) - H23(i,j,k)/seii(1)/2._dl
    ddgup_dh(5,2) = ddgup_dh(5,2) - H13(i,j,k)/seii(2)/2._dl
    ddgup_dh(6,3) = ddgup_dh(6,3) - H12(i,j,k)/seii(3)/2._dl
#endif
    return
  end subroutine get_dgup


!!IMPORTANT: This is a single-thread subroutine
!!DO NOT CALL IT IN A PARALLEL LOOP
  subroutine get_cach(k)  
    integer(IB) i,j,k,ii
#if METRIC_PERTURB
#if DIS_SCHEME == HLATTICE1
    if(k.eq.1)then !!initialize slice #1 and slice #0
       do ii=0,1
          !$omp parallel do default(shared) private(i,j)
          do j=1,n; do i=1,n
             call get_dgup(i,j,ii,sipar%cach_dg(i,j,ii),sipar%cach_dgup(:,i,j,ii),sipar%cach_ddgup(:,:,i,j,ii))
          enddo; enddo
          !$omp end parallel do
       enddo
       sipar%cach_ind(-1) = 2 !!previous slice
       sipar%cach_ind(0) = 0  !!the current slice
       sipar%cach_ind(1) = 1 !!next slice
    endif

    sipar%cach_ind = mod(sipar%cach_ind+1, 3)  !!rotate the indices
    !!calculate the slice #k+1
    !$omp parallel do default(shared) private(i,j)
    do j=1,n; do i=1,n
       call get_dgup(i,j,k+1,sipar%cach_dg(i,j,sipar%cach_ind(1)), sipar%cach_dgup(:,i,j,sipar%cach_ind(1)),sipar%cach_ddgup(:,:,i,j,sipar%cach_ind(1)))
    enddo; enddo
    !$omp end parallel do
#elif DIS_SCHEME == HLATTICE2
    if(k.eq.1)then !!initialize slice #0, #1, #2, #3
       do ii=0,3
          !$omp parallel do default(shared) private(i,j)
          do j=1,n; do i=1,n
             call get_dgup(i,j,ii-1,sipar%cach_dg(i,j,ii),sipar%cach_dgup(:,i,j,ii),sipar%cach_ddgup(:,:,i,j,ii))
          enddo; enddo
          !$omp end parallel do
       enddo
       sipar%cach_ind(-2) = 4
       sipar%cach_ind(-1) = 0 
       sipar%cach_ind(0) = 1  !!the current slice
       sipar%cach_ind(1) = 2
       sipar%cach_ind(2) = 3
    endif

    sipar%cach_ind = mod(sipar%cach_ind+1, 5)  !!rotate the indices
    !!calculate the slice #k+1
    !$omp parallel do default(shared) private(i,j)
    do j=1,n; do i=1,n
       call get_dgup(i,j,k+2,sipar%cach_dg(i,j,sipar%cach_ind(2)), sipar%cach_dgup(:,i,j,sipar%cach_ind(2)),sipar%cach_ddgup(:,:,i,j,sipar%cach_ind(2)))
    enddo; enddo
    !$omp end parallel do
#endif
#endif
    return
  end subroutine get_cach

  subroutine get_energy()
    need_recalculate = .true.
    last_total_detg = total_detg()
    last_potential_energy = potential_energy()
    last_fields_kinetic_energy = fields_kinetic_energy()
    last_fields_gradient_energy = fields_gradient_energy()
    last_gravity_kinetic_energy = gravity_kinetic_energy()
    last_gravity_gradient_energy = gravity_gradient_energy(last_total_detg)
    last_effective_Hubble = effective_Hubble(last_total_detg)
    need_recalculate = .false.
  end subroutine get_energy

  function total_fields_energy()
    real(dl) total_fields_energy
    total_fields_energy = potential_energy() + fields_kinetic_energy() + fields_gradient_energy()
  end function total_fields_energy

  subroutine get_GW_spectrum(S_k)
    real(dl) S_k(fft_numk)
    real(dl),parameter::amp=const_pi/2.*Mplsq/ncube**2
    DEFINE_IND
    real(dl) hpk(fft_numk),ppk(fft_numk),tmp
    S_k=0._dl
#if HIJ_DEFINED
#if METRIC_PERTURB
    !!get the h_{ij}' 
    metric_p = metric_p*(2./Mplsq*metric%physdx)
    !$omp parallel do
    LOOP
    MEPIU(i,j,k) = (2.*MEPIU(i,j,k)-sum(MEPIU(i,j,k)))
    ENDLOOP
    !$omp end parallel do
#else
    metric_p = metric_p*metric%physdx
#endif

#if METRIC_OPTION == FRW_PERTURB_ADAPTIVE
    print*,"coordinate transformation..."
#endif
    call CubicMassiveFFT(6,metric_h,metric_p,FFT_FORWARD)

#ifdef DEBUG_MODE
    call fft_Get_TT_Power(metric_h,metric_p,hpk,ppk,OUTPUT_TT)
    call CubicMassiveFFT(6,metric_h,metric_p,FFT_BACKWARD)
    print*,sqrt(sum(metric_h**2)/ncube)
    tmp=0.
    !$omp parallel do reduction(+:tmp)
    LOOP
    tmp = tmp + sum(metric_h(1:3,i,j,k))**2
    ENDLOOP
    !$omp end parallel do
    print*,sqrt(tmp)
    tmp=0.
    !$omp parallel do reduction(+:tmp)
    LOOP
    tmp = tmp + (TWELVE_DU_X(1,i,j,k)+TWELVE_DV_Y(3,i,j,k)+TWELVE_DV_Z(2,i,j,k))**2 & 
         +(TWELVE_DU_Y(2,i,j,k)+TWELVE_DV_X(3,i,j,k)+TWELVE_DV_Z(1,i,j,k))**2 &
         +(TWELVE_DU_Z(3,i,j,k)+TWELVE_DV_Y(1,i,j,k)+TWELVE_DV_X(2,i,j,k))**2
    ENDLOOP
    !$omp end parallel do
    print*,sqrt(tmp)/12.
    stop
#endif
#if METRIC_OPTION == FRW_PERTURB_ADAPTIVE
    call fft_Get_TT_Power(metric_h,metric_p,hpk,ppk,OUTPUT_ADAPTIVE)
#else
    call fft_Get_TT_Power(metric_h,metric_p,hpk,ppk,UNCHANGED)
#endif
    do i=1,fft_numk
       S_k(i)=((i*const_2pi/n)**2*hpk(i) + ppk(i))/metric%physdx**2*amp*fft_vol(i)*real(i,dl)/(4._dl*const_pi)
    enddo
    call CubicMassiveFFT(6,metric_h,metric_p,FFT_BACKWARD)
#if METRIC_OPTION ==  FRW_PERTURB_ADAPTIVE
!    print*,"After coordinate transformation <h_{ij}^2>^{1/2}~", sqrt(sum(metric_h**2)/ncube)
    LOOP
#if DIS_SCHEME == HLATTICE1
    fields_f(:,i,j,k) =  fields_f(:,i,j,k) - (TWO_DFLD_X(i,j,k)*metric_p(1,i,j,k) + TWO_DFLD_Y(i,j,k)*metric_p(2,i,j,k) + TWO_DFLD_Z(i,j,k)*metric_p(3,i,j,k))/2.
    fields_p(:,i,j,k) =  fields_p(:,i,j,k) - (TWO_DF_X_FLD(fields_p,:,i,j,k)*metric_p(1,i,j,k) + TWO_DF_Y_FLD(fields_p,:,i,j,k)*metric_p(2,i,j,k) + TWO_DF_Z_FLD(fields_p,:,i,j,k)*metric_p(3,i,j,k))/2.
#elif DIS_SCHEME == HLATTICE2
    fields_f(:,i,j,k) =  fields_f(:,i,j,k) - (TWELVE_DFLD_X(i,j,k)*metric_p(1,i,j,k) + TWELVE_DFLD_Y(i,j,k)*metric_p(2,i,j,k) + TWELVE_DFLD_Z(i,j,k)*metric_p(3,i,j,k))/12.
    fields_p(:,i,j,k) =  fields_p(:,i,j,k) - (TWELVE_DF_X_FLD(fields_p,:,i,j,k)*metric_p(1,i,j,k) + TWELVE_DF_Y_FLD(fields_p,:,i,j,k)*metric_p(2,i,j,k) + TWELVE_DF_Z_FLD(fields_p,:,i,j,k)*metric_p(3,i,j,k))/12.
#else
    stop "For METRIC_OPTION=FRW_PERTURB_ADPATIVE you must use DIS_SCHEME = HLATTICE1 or HLATTICE2"
#endif
    ENDLOOP
    !!reload the metric
    call metric_load_term("data/"//trim(run_name)//"_metric_pij.log",2)
#else
#if METRIC_PERTURB
    metric_p = metric_p*(Mplsq/2./ metric%physdx)
    !$omp parallel do
    LOOP
    MEPIU(i,j,k) = (MEPIU(i,j,k)-sum(MEPIU(i,j,k)))/2.
    ENDLOOP
    !$omp end parallel do
#else
    metric_p = metric_p/ metric%physdx
#endif
#endif    
#endif
  end subroutine get_GW_spectrum


end module mutils


