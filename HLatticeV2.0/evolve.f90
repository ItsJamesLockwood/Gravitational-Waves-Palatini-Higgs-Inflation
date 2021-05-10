!!the kernal of the code
!!subroutines that evolve the system
module evolution
  use io_utils
  use mutils
  implicit none
#include "configure.h"

contains

  subroutine evolve()
    logical save_slices, save_fields, save_metric
    save_fields = .false.
    save_slices = .false.
    save_metric = .false.
#if WANTFIELDS
  save_fields= .true.
#endif
#if WANTSLICES
  save_slices = .true.
#endif
#if WANTMETRIC && (WANTGW || METRIC_PERTURB)
  save_metric = .true.
#endif
    do
      call output_to_screen()
      call output_energy()
      if(use_checkpoint .and. mod(sipar%nsteps,checkpoint_steps).eq.0 .and. (sipar%nsteps.gt.0 .or. write_check_at_step0))then
        call write_check()
        if (save_fields .and. mod(sipar%nsteps,checkpoint_steps*save_field_interval).eq.0 .and. (checkpoint_steps*save_field_interval*field_number_cutoff).ge.sipar%nsteps ) then
          write(*,*) "Saving fields at step ",sipar%nsteps, "..."
          call output_fields()
          call output_momenta()
        endif 
        if (save_slices .and. mod(sipar%nsteps,checkpoint_steps*save_slice_interval) .eq. 0) then 
          write(*,*) "Saving slices at step ",sipar%nsteps, "..."
          call output_field_slice()
          call output_momentum_slice()
        endif
        if (save_metric .and.  mod(sipar%nsteps,checkpoint_steps*save_field_interval).eq.0 .and. (checkpoint_steps*save_field_interval*field_number_cutoff).ge.sipar%nsteps ) then
          write(*,*) "Saving metric at step ",sipar%nsteps, "..."
          call output_metric_h()
        endif        
      endif
      call step(n_feedback)
      
      !! Troubleshooting statements: comment out if not testing initial values.
       if(sipar%nsteps .ge. stop_at_step .or. metric%a .ge. stop_at_a .or. isnan(metric%a))then
          if(isnan(metric%a)) then
            write(*,*) "Metric a is NaN. Exiting..."
          else if(metric%a .ge. stop_at_a)then
#if METRIC_PERTURB
            write(*,*) "Reached max a. Exiting..."
#else
            write(*,*) "Reached max a at a=",metric%a," Exiting..."
            write(*,*) "Mode: no perturbations"
#endif
          else
            write(*,*) "Reached max number of steps. Exiting..."
          end if
          call final_output()
          exit
       endif
       call coor_trans()
       !!write(*,*) "Exiting after first lattice step due to TS-mode..."
       !!exit
    enddo
    return
  end subroutine evolve

  subroutine step(nsteps)
    integer(IB) nsteps
    need_recalculate = .true.
    sipar%dt = model_dt_per_step()

#if INTEGRATOR == SYMPLECTIC_2ND
    call symp2(sipar%dt,nsteps)
#elif INTEGRATOR == SYMPLECTIC_4TH
    call symp4(sipar%dt,nsteps)
#elif INTEGRATOR == SYMPLECTIC_6TH
    call symp6(sipar%dt,nsteps)
#else
    stop "unknown integrator"
#endif
#if WANTGW && ! METRIC_PERTURB
    call NoPerturbEvolveGW(sipar%dt*nsteps)
#endif
    sipar%nsteps = sipar%nsteps + nsteps
    sipar%time = sipar%time + sipar%dt*nsteps
  end subroutine step

  subroutine symp_o2step(dt,c1,c2)
    real(dl) dt,c1,c2
    integer i
    do i=2,n_Hamiltonian_terms-1
       call Hamiltonian_Split(c1*dt/2._dl,i)
    enddo
    call Hamiltonian_Split(c1*dt, n_Hamiltonian_terms)
    do i=n_Hamiltonian_terms-1,2,-1
       call Hamiltonian_Split(c1*dt/2._dl, i)
    enddo
    call Hamiltonian_Split((c1+c2)*dt/2._dl,1)
    return
  end subroutine symp_o2step

  subroutine symp2(dt,nsteps)
    real(dl) dt
    integer(IB) nsteps,i,j
    call Hamiltonian_Split(dt/2._dl, 1)
    do j=1,nsteps-1
       call symp_o2step(dt,1._dl, 1._dl)
    enddo
    call symp_o2step(dt,1._dl, 0._dl)
  end subroutine symp2

  subroutine symp4(dt,nsteps)
    real(dl) dt
    integer(IB) nsteps,i,j
    real(dl),parameter:: c1 = 1._dl/(2._dl - 2._dl**(1._dl/3._dl))
    real(dl),parameter:: c0 = 1._dl - 2._dl*c1
    call Hamiltonian_Split(c1*dt/2._dl,1)
    do j=1,nsteps
       call symp_o2step(dt, c1, c0)
       call symp_o2step(dt, c0, c1)
       if(j.eq.nsteps)then
          call symp_o2step(dt, c1, 0._dl)
       else
          call symp_o2step(dt, c1, c1)
       endif
    enddo
  end subroutine symp4

  subroutine symp6(dt,nsteps)
    real(dl) dt
    integer(IB) nsteps,i,j
    real(dl),parameter:: c1 = -1.17767998417887_dl, c2 = 0.235573213359357_dl, c3 = 0.784513610477560_dl
    real(dl),parameter:: c0 = 1._dl-2._dl*(c1+c2+c3)

    call Hamiltonian_Split(c3*dt/2._dl,1)

    do j=1,nsteps
       call symp_o2step(dt, c3, c2)
       call symp_o2step(dt, c2, c1)
       call symp_o2step(dt, c1, c0)
       call symp_o2step(dt, c0, c1)
       call symp_o2step(dt, c1, c2)
       call symp_o2step(dt, c2, c3)
       if(j.eq.nsteps)then
          call symp_o2step(dt, c3, 0._dl)
       else
          call symp_o2step(dt, c3, c3)
       endif
    enddo
  end subroutine symp6

!!===================== Hamiltonian splitting ===============================
  subroutine Hamiltonian_Split(dt, iterm)
    integer(IB) iterm
    real(dl) dt
    select case(iterm)
    case(1)
       call Hamiltonian_gravity_kinetic_diagonal(dt)
    case(2)
       call Hamiltonian_fields_kinetic(dt)
       call Hamiltonian_gravity_kinetic_offdiagonal(dt)
    case(3)
       call Hamiltonian_potential_and_gradient(dt)
    case default
       write(*,*) "Unknown Hamiltonian term#: "//trim(int2str(iterm))
       stop
    end select
    
  end subroutine Hamiltonian_Split

  subroutine Hamiltonian_fields_kinetic(dt)
    real(dl) dt, ss, ss2
    DEFINE_IND
    !!fields kinetic term :  pi_phi ^2/(2 sqrt(|g|))
    ss = dt
    ss2 = dt/4._dl
#if METRIC_OPTION == FRW_BACKGROUND 
    ss = ss/metric%y**2
#if USE_CONFORMAL_TIME
    metric%piy = metric%piy + (2._dl*metric%y**3)* fields_kinetic_energy() * dt
#else
    metric%piy = metric%piy + (2._dl*metric%y)* fields_kinetic_energy() * dt
#endif
#endif
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    fields_f(:,i,j,k) = fields_f(:,i,j,k) + fields_p(:,i,j,k) * ss / DETG(i,j,k)
#if METRIC_PERTURB
    MEPIU(i,j,k) = MEPIU(i,j,k) + sum(fields_p(:,i,j,k)**2) * ss2 / DETG(i,j,k)
#endif
    ENDLOOP
    !$omp end parallel do
    return
  end subroutine Hamiltonian_fields_kinetic

  subroutine Hamiltonian_gravity_kinetic_diagonal(dt)
    real(dl) dt, ss, ss2
    real(dl) tv1(6),tv2(6),tv3(6),tv4(6) !!for Runge-Kutta 
    DEFINE_IND
    integer(IB) irk
    !!gravitaty kinetic term
#if METRIC_PERTURB
    !!pi_u terms, noncanonical, use Runge-Kutta
    ss = dt/Mplsq/noncanonical_Runge_Kutta_steps
    ss2 = ss/2._dl
    do irk = 1,noncanonical_Runge_Kutta_steps
       !$omp parallel do default(shared) private(i,j,k,tv1,tv2,tv3,tv4)
       LOOP
       tv1(1:3) = (2._dl*MEPIU(i,j,k) - sum(MEPIU(i,j,k)))/DETG(i,j,k)*ss
       tv1(4:6) =  (sum(MEPIU(i,j,k)**2) - 0.5_dl*sum(MEPIU(i,j,k))**2 )/DETG(i,j,k)*ss2
       tv2(1:3) = (2._dl*(MEPIU(i,j,k)+tv1(4:6)) - sum(MEPIU(i,j,k)+tv1(4:6)))*fast_exp(-sum(MEU(i,j,k)+tv1(1:3))/2._dl)*ss
       tv2(4:6) = (sum((MEPIU(i,j,k)+tv1(4:6))**2) - 0.5_dl*sum(MEPIU(i,j,k)+tv1(4:6))**2 )*fast_exp(-sum(MEU(i,j,k)+tv1(1:3))/2._dl)* ss2 
       tv3(1:3) = (4._dl*(MEPIU(i,j,k)+tv2(4:6)) - 2._dl*sum(MEPIU(i,j,k)+tv2(4:6)))*fast_exp(-sum(MEU(i,j,k)+tv2(1:3))/2._dl)*ss
       tv3(4:6) = (sum((MEPIU(i,j,k)+tv2(4:6))**2) - 0.5_dl*sum(MEPIU(i,j,k)+tv2(4:6))**2 )*fast_exp(-sum(MEU(i,j,k)+tv2(1:3))/2._dl)*ss
       tv4(1:3) = (2._dl*(MEPIU(i,j,k)+tv3(4:6)) - sum(MEPIU(i,j,k)+tv3(4:6)))*fast_exp(-sum(MEU(i,j,k)+tv3(1:3))/2._dl)*ss
       tv4(4:6) = (sum((MEPIU(i,j,k)+tv3(4:6))**2) - 0.5_dl*sum(MEPIU(i,j,k)+tv3(4:6))**2 )*fast_exp(-sum(MEU(i,j,k)+tv3(1:3))/2._dl)*ss2
       MEU(i,j,k) = MEU(i,j,k) + (tv1(1:3)+2.*tv2(1:3)+tv3(1:3)+tv4(1:3))/3._dl
       MEPIU(i,j,k) = MEPIU(i,j,k) + (tv1(4:6)+2.*tv2(4:6)+tv3(4:6)+tv4(4:6))/3._dl
       ENDLOOP
       !$omp end parallel do
    enddo
    return
#endif
#if METRIC_OPTION == FRW_BACKGROUND
#if USE_CONFORMAL_TIME
    metric%y = metric%y - metric%piy/Mplsq/6._dl *dt
    !write(*,*) "metric%y:",metric%y, "; K_f:", fields_kinetic_energy()

#else
    metric%y = metric%y - (3._dl/8._dl)*metric%piy/Mplsq *dt
#endif
#endif
  end subroutine Hamiltonian_gravity_kinetic_diagonal

  subroutine Hamiltonian_gravity_kinetic_offdiagonal(dt)
    real(dl) dt, ss, ss2
    DEFINE_IND
#if METRIC_PERTURB
    !!the pi_v terms
    ss = 2._dl*dt/Mplsq 
    ss2 = ss/4._dl
    !$omp parallel do default(shared) private(i,j,k)
    LOOP
    MEV(i,j,k) = MEV(i,j,k) + ss * MEPIV(i,j,k)/DETG(i,j,k)
    MEPIU(i,j,k) = MEPIU(i,j,k) + ss2 * sum(MEPIV(i,j,k)**2)/DETG(i,j,k)
    ENDLOOP
    !$omp end parallel do
#endif
    return
  end subroutine Hamiltonian_gravity_kinetic_offdiagonal

  !!=====================================================================
  !!gradient energy and potential energy terms, dt, the most complicated  part
  subroutine Hamiltonian_potential_and_gradient(dt)
    real(dl) dt, ss, ss2,ss4,ss3, ss5,gr(6), totdetg,crdetg, ge
    DEFINE_IND
#if METRIC_OPTION == FRW_BACKGROUND
    metric%a=scale_factor()
    metric%physdx=metric%dx*metric%a
#if USE_CONFORMAL_TIME
    metric%piy = metric%piy - (4._dl*potential_energy() + 2.*fields_gradient_energy())*metric%y**3 * dt
#else
    metric%piy = metric%piy - (2._dl*potential_energy() + (2._dl/3._dl)*fields_gradient_energy())*metric%y * dt
#endif
#endif
    totdetg= total_detg()
    crdetg=totdetg**(1./3._dl)
    ss4 = Mplsq * dt/ 16._dl / metric%physdx**2 * crdetg
    ss3=ss4 * 2._dl
    ge = gravity_gradient_energy(totdetg)
    ss5 = ge/totdetg/6._dl *dt

    ss = dt
#if METRIC_OPTION == FRW_BACKGROUND
#if USE_CONFORMAL_TIME
    ss = ss*metric%y**4
#else
    ss = ss*metric%y**2
#endif
#endif

#if DIS_SCHEME == LATTICEEASY
    ss2 = ss/metric%physdx**2
#elif DIS_SCHEME == HLATTICE1
    ss2 = ss/metric%physdx**2/4.
#elif DIS_SCHEME == HLATTICE2
    ss2 = ss/metric%physdx**2/144.
    ss4 = ss4/36.
    ss3 = ss3/36.
#endif

    do k=1,n
       call get_cach(k)
       !$omp parallel do default(shared) private(i,j,gr)
       do j=1,n; do i=1,n
          fields_p(:,i,j,k) = fields_p(:,i,j,k) + ss*( &
               - dVdf(fields_f(:,i,j,k))*CACH_DG(i,j,0)) &
               + ss2* &
#if METRIC_PERTURB
#if DIS_SCHEME == HLATTICE1
               (TWO_DFLD_X(i+1,j,k) *CACH_DGUP11(i+1,j,0) &
               -TWO_DFLD_X(i-1,j,k) *CACH_DGUP11(i-1,j,0)   &
               +TWO_DFLD_Y(i,j+1,k) *CACH_DGUP22(i,j+1,0) &
               -TWO_DFLD_Y(i,j-1,k) *CACH_DGUP22(i,j-1,0)   &
               +TWO_DFLD_Z(i,j,k+1) *CACH_DGUP33(i,j,1)   &
               -TWO_DFLD_Z(i,j,k-1) *CACH_DGUP33(i,j,-1)   &
               +TWO_DFLD_Y(i+1,j,k) *CACH_DGUP12(i+1,j,0) &
               -TWO_DFLD_Y(i-1,j,k) *CACH_DGUP12(i-1,j,0) &
               +TWO_DFLD_X(i,j+1,k) *CACH_DGUP12(i,j+1,0) &
               -TWO_DFLD_X(i,j-1,k) *CACH_DGUP12(i,j-1,0) &
               +TWO_DFLD_Z(i+1,j,k) *CACH_DGUP13(i+1,j,0) &
               -TWO_DFLD_Z(i-1,j,k) *CACH_DGUP13(i-1,j,0) &
               +TWO_DFLD_X(i,j,k+1) *CACH_DGUP13(i,j,1) &
               -TWO_DFLD_X(i,j,k-1) *CACH_DGUP13(i,j,-1) &
               +TWO_DFLD_Z(i,j+1,k) *CACH_DGUP23(i,j+1,0) &
               -TWO_DFLD_Z(i,j-1,k) *CACH_DGUP23(i,j-1,0) &
               +TWO_DFLD_Y(i,j,k+1) *CACH_DGUP23(i,j,1) &
               -TWO_DFLD_Y(i,j,k-1) *CACH_DGUP23(i,j,-1) &
               )
#elif DIS_SCHEME == HLATTICE2
               ((TWELVE_DFLD_X(i+1,j,k) *CACH_DGUP11(i+1,j,0) &
               -TWELVE_DFLD_X(i-1,j,k) *CACH_DGUP11(i-1,j,0)   &
               +TWELVE_DFLD_Y(i,j+1,k) *CACH_DGUP22(i,j+1,0) &
               -TWELVE_DFLD_Y(i,j-1,k) *CACH_DGUP22(i,j-1,0)   &
               +TWELVE_DFLD_Z(i,j,k+1) *CACH_DGUP33(i,j,1)   &
               -TWELVE_DFLD_Z(i,j,k-1) *CACH_DGUP33(i,j,-1)   &
               +TWELVE_DFLD_Y(i+1,j,k) *CACH_DGUP12(i+1,j,0) &
               -TWELVE_DFLD_Y(i-1,j,k) *CACH_DGUP12(i-1,j,0) &
               +TWELVE_DFLD_X(i,j+1,k) *CACH_DGUP12(i,j+1,0) &
               -TWELVE_DFLD_X(i,j-1,k) *CACH_DGUP12(i,j-1,0) &
               +TWELVE_DFLD_Z(i+1,j,k) *CACH_DGUP13(i+1,j,0) &
               -TWELVE_DFLD_Z(i-1,j,k) *CACH_DGUP13(i-1,j,0) &
               +TWELVE_DFLD_X(i,j,k+1) *CACH_DGUP13(i,j,1) &
               -TWELVE_DFLD_X(i,j,k-1) *CACH_DGUP13(i,j,-1) &
               +TWELVE_DFLD_Z(i,j+1,k) *CACH_DGUP23(i,j+1,0) &
               -TWELVE_DFLD_Z(i,j-1,k) *CACH_DGUP23(i,j-1,0) &
               +TWELVE_DFLD_Y(i,j,k+1) *CACH_DGUP23(i,j,1) &
               -TWELVE_DFLD_Y(i,j,k-1) *CACH_DGUP23(i,j,-1) &
               )*8.- ( &
               TWELVE_DFLD_X(i+2,j,k) *CACH_DGUP11(i+2,j,0) &
               -TWELVE_DFLD_X(i-2,j,k) *CACH_DGUP11(i-2,j,0)   &
               +TWELVE_DFLD_Y(i,j+2,k) *CACH_DGUP22(i,j+2,0) &
               -TWELVE_DFLD_Y(i,j-2,k) *CACH_DGUP22(i,j-2,0)   &
               +TWELVE_DFLD_Z(i,j,k+2) *CACH_DGUP33(i,j,2)   &
               -TWELVE_DFLD_Z(i,j,k-2) *CACH_DGUP33(i,j,-2)   &
               +TWELVE_DFLD_Y(i+2,j,k) *CACH_DGUP12(i+2,j,0) &
               -TWELVE_DFLD_Y(i-2,j,k) *CACH_DGUP12(i-2,j,0) &
               +TWELVE_DFLD_X(i,j+2,k) *CACH_DGUP12(i,j+2,0) &
               -TWELVE_DFLD_X(i,j-2,k) *CACH_DGUP12(i,j-2,0) &
               +TWELVE_DFLD_Z(i+2,j,k) *CACH_DGUP13(i+2,j,0) &
               -TWELVE_DFLD_Z(i-2,j,k) *CACH_DGUP13(i-2,j,0) &
               +TWELVE_DFLD_X(i,j,k+2) *CACH_DGUP13(i,j,2) &
               -TWELVE_DFLD_X(i,j,k-2) *CACH_DGUP13(i,j,-2) &
               +TWELVE_DFLD_Z(i,j+2,k) *CACH_DGUP23(i,j+2,0) &
               -TWELVE_DFLD_Z(i,j-2,k) *CACH_DGUP23(i,j-2,0) &
               +TWELVE_DFLD_Y(i,j,k+2) *CACH_DGUP23(i,j,2) &
               -TWELVE_DFLD_Y(i,j,k-2) *CACH_DGUP23(i,j,-2) &
               ))
#endif
#else
#if DIS_SCHEME == LATTICEEASY
               LAPLACIAN(i,j,k)
#elif DIS_SCHEME == HLATTICE1
               FOUR_LAP(i,j,k)
#elif DIS_SCHEME == HLATTICE2
               I44_LAP(i,j,k)
#endif
#endif


          !!gravity waves are generated here
#if METRIC_PERTURB
#if DIS_SCHEME == HLATTICE1
          gr(1) = sum(TWO_DFLD_X(i,j,k)**2)/2.
          gr(2) = sum(TWO_DFLD_Y(i,j,k)**2)/2.
          gr(3) = sum(TWO_DFLD_Z(i,j,k)**2)/2.
          gr(4) = sum(TWO_DFLD_Y(i,j,k)*TWO_DFLD_Z(i,j,k))
          gr(5) = sum(TWO_DFLD_X(i,j,k)*TWO_DFLD_Z(i,j,k))
          gr(6) = sum(TWO_DFLD_X(i,j,k)*TWO_DFLD_Y(i,j,k))
#elif DIS_SCHEME == HLATTICE2
          gr(1) = sum(TWELVE_DFLD_X(i,j,k)**2)/2.
          gr(2) = sum(TWELVE_DFLD_Y(i,j,k)**2)/2.
          gr(3) = sum(TWELVE_DFLD_Z(i,j,k)**2)/2.
          gr(4) = sum(TWELVE_DFLD_Y(i,j,k)*TWELVE_DFLD_Z(i,j,k))
          gr(5) = sum(TWELVE_DFLD_X(i,j,k)*TWELVE_DFLD_Z(i,j,k))
          gr(6) = sum(TWELVE_DFLD_X(i,j,k)*TWELVE_DFLD_Y(i,j,k))
#endif
          MEPIU(i,j,k) = MEPIU(i,j,k) - 0.5_dl*potential(fields_f(:,i,j,k))*CACH_DG(i,j,0)*ss 
          
          metric_p(:,i,j,k) = metric_p(:,i,j,k) - matmul(gr,sipar%cach_ddgup(:,:,i,j,sipar%cach_ind(0)))*ss2

#if DIS_SCHEME == HLATTICE1
          metric_p(1,i,j,k) = metric_p(1,i,j,k) - ss4*( &
               FOUR_D2U_ZZ(2,i,j,k)  + FOUR_D2U_YY(3,i,j,k) - 2.*FOUR_D2V_YZ(1,i,j,k))
          metric_p(2,i,j,k) = metric_p(2,i,j,k) - ss4*( &
                FOUR_D2U_ZZ(1,i,j,k) + FOUR_D2U_XX(3,i,j,k) - 2.*FOUR_D2V_XZ(2,i,j,k))
          metric_p(3,i,j,k) = metric_p(3,i,j,k) - ss4*( &
               FOUR_D2U_XX(2,i,j,k) + FOUR_D2U_YY(1,i,j,k) - 2.*FOUR_D2V_XY(3,i,j,k))
          metric_p(4,i,j,k) = metric_p(4,i,j,k) - ss3*( &
               -FOUR_D2V_XX(1,i,j,k) + FOUR_D2V_XY(2,i,j,k) + FOUR_D2V_XZ(3,i,j,k)- FOUR_D2U_YZ(1,i,j,k))
          metric_p(5,i,j,k) = metric_p(5,i,j,k) - ss3*( &
               -FOUR_D2V_YY(2,i,j,k) + FOUR_D2V_XY(1,i,j,k) + FOUR_D2V_YZ(3,i,j,k)- FOUR_D2U_XZ(2,i,j,k))
          metric_p(6,i,j,k) = metric_p(6,i,j,k) - ss3*( &
               -FOUR_D2V_ZZ(3,i,j,k) +FOUR_D2V_XZ(1,i,j,k) + FOUR_D2V_YZ(2,i,j,k)- FOUR_D2U_XY(3,i,j,k))
#elif DIS_SCHEME == HLATTICE2
          metric_p(1,i,j,k) = metric_p(1,i,j,k) - ss4*( &
               I44_D2U_ZZ(2,i,j,k)  + I44_D2U_YY(3,i,j,k) - 2.*I44_D2V_YZ(1,i,j,k))
          metric_p(2,i,j,k) = metric_p(2,i,j,k) - ss4*( &
                I44_D2U_ZZ(1,i,j,k) + I44_D2U_XX(3,i,j,k) - 2.*I44_D2V_XZ(2,i,j,k))
          metric_p(3,i,j,k) = metric_p(3,i,j,k) - ss4*( &
               I44_D2U_XX(2,i,j,k) + I44_D2U_YY(1,i,j,k) - 2.*I44_D2V_XY(3,i,j,k))
          metric_p(4,i,j,k) = metric_p(4,i,j,k) - ss3*( &
               -I44_D2V_XX(1,i,j,k) + I44_D2V_XY(2,i,j,k) + I44_D2V_XZ(3,i,j,k)- I44_D2U_YZ(1,i,j,k))
          metric_p(5,i,j,k) = metric_p(5,i,j,k) - ss3*( &
               -I44_D2V_YY(2,i,j,k) + I44_D2V_XY(1,i,j,k) + I44_D2V_YZ(3,i,j,k)- I44_D2U_XZ(2,i,j,k))
          metric_p(6,i,j,k) = metric_p(6,i,j,k) - ss3*( &
               -I44_D2V_ZZ(3,i,j,k) +I44_D2V_XZ(1,i,j,k) + I44_D2V_YZ(2,i,j,k)- I44_D2U_XY(3,i,j,k))
#endif

          MEPIU(i,j,k) = MEPIU(i,j,k) - ss5*CACH_DG(i,j,0)
#endif
       enddo;enddo
       !$omp end parallel do
    enddo
  end subroutine Hamiltonian_potential_and_gradient

!!for no metric feedback 
!!evolve h and \dot h
  subroutine NoPerturbEvolveGW(dt)
    real(dl) dt
#if WANTGW && ! METRIC_PERTURB
    real(dl) H3, a,a2,physdx2,m2,m4,dphyst,fourphysdx2,I44Physdx2
    DEFINE_IND
    a = scale_factor()
    a2 = a*a
    physdx2 = metric%physdx**2
    H3 = 3.*effective_Hubble()
    m2 = 2./Mplsq/physdx2
#if DIS_SCHEME == HLATTICE1
    fourphysdx2 = physdx2*4.
    m2 = m2/4.
#elif DIS_SCHEME == HLATTICE2
    I44physdx2 = physdx2*144.
    m2 = m2/144.
#endif

#if USE_CONFORMAL_TIME
    dphyst = dt*a
#else
    dphyst = dt
#endif

    !$omp parallel do private(i,j,k) default(shared)
    LOOP
    metric_h(:,i,j,k) = metric_h(:,i,j,k) + dphyst * metric_p(:,i,j,k)
    ENDLOOP
    !$omp end parallel do

    !$omp parallel do private(i,j,k) default(shared)
    LOOP
#if DIS_SCHEME == LATTICEEASY
    metric_p(:,i,j,k) = metric_p(:,i,j,k) + dphyst * ( -metric_p(:,i,j,k)*H3 + (metric_h(:,sind(i+1),j,k)+metric_h(:,sind(i-1),j,k)+metric_h(:,i,sind(j+1),k)+metric_h(:,i,sind(j-1),k)+metric_h(:,i,j,sind(k+1))+metric_h(:,i,j,sind(k-1))-6.*metric_h(:,i,j,k))/physdx2 + m2* (/  &
         sum((fields_f(:,sind(i+1),j,k)-fields_f(:,i,j,k))**2), & 
         sum((fields_f(:,i,sind(j+1),k)-fields_f(:,i,j,k))**2), &
         sum((fields_f(:,i,j,sind(k+1))-fields_f(:,i,j,k))**2), &
         sum((fields_f(:,i, sind(j+1),k)-fields_f(:,i,j,k))*(fields_f(:,i,j,sind(k+1))-fields_f(:,i,j,k))) , &
         sum((fields_f(:,sind(i+1),j,k)-fields_f(:,i,j,k))*(fields_f(:,i,j,sind(k+1))-fields_f(:,i,j,k))) , &
         sum((fields_f(:,sind(i+1),j,k)-fields_f(:,i,j,k))*(fields_f(:,i,sind(j+1),k)-fields_f(:,i,j,k))) &
        /) )
#elif DIS_SCHEME == HLATTICE1
    metric_p(:,i,j,k) = metric_p(:,i,j,k) + dphyst * ( -metric_p(:,i,j,k)*H3 +FOUR_NABLA2_FLD(metric_h,:,i,j,k)/fourphysdx2 + m2* (/  &
         sum(TWO_DFLD_X(i,j,k)**2), & 
         sum(TWO_DFLD_Y(i,j,k)**2), &
         sum(TWO_DFLD_Z(i,j,k)**2), &
         sum(TWO_DFLD_Y(i,j,k)*TWO_DFLD_Z(i,j,k)) , &
         sum(TWO_DFLD_Z(i,j,k)*TWO_DFLD_X(i,j,k)) , &
         sum(TWO_DFLD_X(i,j,k)*TWO_DFLD_Y(i,j,k)) &
        /) )
#elif DIS_SCHEME == HLATTICE2

    metric_p(:,i,j,k) = metric_p(:,i,j,k) + dphyst * ( -metric_p(:,i,j,k)*H3 +I44_NABLA2_FLD(metric_h,:,i,j,k)/I44physdx2 + m2* (/  &
         sum(TWELVE_DFLD_X(i,j,k)**2), & 
         sum(TWELVE_DFLD_Y(i,j,k)**2), &
         sum(TWELVE_DFLD_Z(i,j,k)**2), &
         sum(TWELVE_DFLD_Y(i,j,k)*TWELVE_DFLD_Z(i,j,k)) , &
         sum(TWELVE_DFLD_Z(i,j,k)*TWELVE_DFLD_X(i,j,k)) , &
         sum(TWELVE_DFLD_X(i,j,k)*TWELVE_DFLD_Y(i,j,k)) &
        /) )
#endif
    ENDLOOP
    !$omp end parallel do

    !!the following operation does not change h_{ij} in Fourier space (for
    !nonzero k), but it significantly reduces numerical noises
    do i=1,6
       metric_p(i,:,:,:) = metric_p(i,:,:,:) - sum(metric_p(i,:,:,:))/ncube
    enddo

#endif
    return
  end subroutine NoPerturbEvolveGW

end module evolution

