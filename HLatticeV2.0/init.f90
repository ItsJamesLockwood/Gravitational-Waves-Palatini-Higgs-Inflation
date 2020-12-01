module init
  use LatticeFFT
  use mutils
  use io_utils
  use model
  implicit none

#include "configure.h"


contains

#define HUBBLE yp(2*ns+1)
#define LNA y(2*ns+1)
#define THE_FIELDS y(1:ns)
#define THE_MOMENTA y(ns+1:2*ns)  
  subroutine initialize()
    real(dl) y(2*ns+1)
    real(dl) t, tnext
    logical ierror
    DEFINE_IND

    call fft_init()

    call allocate_fields()    
    call allocate_metric()

    if(use_checkpoint)then
       call read_check(ierror)
       if(.not.ierror) return
    endif

    call model_output(1)
    call begin_output()
    t = 0._dl
    LNA = 0._dl

    THE_FIELDS = init_fields
    if(all(init_momenta.ge.Mplsq))then
       THE_MOMENTA = - dVdf(THE_FIELDS) * Mpl/sqrt(3._dl * (potential(THE_FIELDS)) ) 
    else
       THE_MOMENTA = init_momenta
    endif

    do while(.not.Start_Box_Simulation(THE_FIELDS, THE_MOMENTA))
       tnext = t + 0.0005_dl * Mpl/sqrt((potential(THE_FIELDS)+sum(THE_MOMENTA**2)/2.)/3._dl)
       call dverk(ode_background,t,y,tnext,1.e-7_dl)
    enddo
    init_fields = THE_FIELDS
    init_momenta = THE_MOMENTA
    call print_settings_file()
    call print_info()

    sipar%time = 0
    sipar%nsteps = 0
    metric%y = 1._dl
    metric%a = 1._dl
    init_Hubble = sqrt((potential(THE_FIELDS)+sum(THE_MOMENTA**2)/2.)/3.)/Mpl
    metric%dx = boxsize_H/init_Hubble/n
    metric%physdx=metric%dx*metric%a
    sipar%dt = model_dt_per_step()
    do fld=1,ns
       if(do_init(fld))call InitSpectrum(fld) 
       write(*,*) "field # "//trim(int2str(fld))//" initialized"
       fields_f(fld,:,:,:) = fields_f(fld,:,:,:) + init_fields(fld) 
       fields_p(fld,:,:,:) = fields_p(fld,:,:,:) + init_momenta(fld)
    enddo

#if METRIC_PERTURB
    if(.not. allocated(metric_h)) stop "h not defined"
    metric_h = 0._dl
#else
#if WANTGW
    if(.not. allocated(metric_h)) stop "h not defined"
    metric_h = 0._dl
    if(.not. allocated(metric_p)) stop "p not defined"
    metric_p = 0._dl
#endif
#endif

    print*,"Initializing the metric..."

    Init_Hubble = sqrt((potential_energy() + fields_kinetic_energy() +fields_gradient_energy()+gravity_gradient_energy())/3._dl)/Mpl
    write(*,*)"Initial Hubble (10^{-8}M_p):", Init_Hubble*1.e8_dl/Mpl
#if METRIC_PERTURB
    metric_p(1:3,:,:,:) = -Mplsq*Init_Hubble
    metric_p(4:6,:,:,:) = 0._dl
#else
#if USE_CONFORMAL_TIME
    metric%piy = -6._dl * Mplsq*metric%y**2 * Init_Hubble
#else
    metric%piy = -4._dl * Mplsq*metric%y * Init_Hubble
#endif
#endif
    
    call model_output(2)

  end subroutine initialize


  subroutine InitSpectrum(fld) 
    !!initialization in k space
    real(dl) repk(fft_numk),impk(fft_numk),amp(2)
    integer i,fld
    do i=1,fft_numk
       amp=Model_Power(fld,i*k_unit/metric%physdx)
       repk(i)=amp(1)
       impk(i)=amp(2)
    enddo     
    repk=repk*(n/metric%dx)**3 
    impk=impk*(n/metric%dx)**3*metric%physdx**2
    fft_ind_re =fld
    fft_ind_im =fld
    call FFT_RandomGaussian(fields_f,fields_p,repk,impk)
    fields_p(fld,:,:,:)=fields_p(fld,:,:,:)/metric%physdx
#ifdef DEBUG_MODE
    print*,"In debug mode: the initial fields fluctuations are amplified by a factor of 1000"
    fields_f(fld,:,:,:)=fields_f(fld,:,:,:)*1000.
#endif
  end subroutine InitSpectrum

  subroutine ode_background(ny,t,y,yp)
    integer(IB) ny
    real(dl),dimension(2*ns+1)::y,yp
    real(dl) t
    yp(1:ns) = THE_MOMENTA
    HUBBLE = sqrt((potential(THE_FIELDS) + sum(THE_MOMENTA**2)/2._dl)/3._dl)/Mpl
    yp(ns+1:2*ns) = - dVdf(THE_FIELDS) - (3._dl * HUBBLE) * THE_MOMENTA 
  end subroutine ode_background

  subroutine print_info()
    call print_bar()
    write(*,*) "Start simulations when the field values (M_p) are"
    call print_fields(init_fields/Mpl)
    write(*,*) " and the field momenta (10^{-8} M_p^2) are"
    call print_fields(init_momenta*1.e8_dl/Mpl**2)
    call print_bar()
    write(*,*) "Simulation resolution: ",n
    call print_bar()
#if INTEGRATOR == SYMPLECTIC_2ND
       write(*,*) "Using 2nd order symplectic integrator"
#elif INTEGRATOR == SYMPLECTIC_4TH
       write(*,*) "Using 4th order symplectic integrator"
#elif INTEGRATOR == SYMPLECTIC_6TH
       write(*,*) "Using 6th order symplectic integrator"
#else
       stop "INTEGRATOR must be SYMPLECTIC_2ND, SYMPLECTIC_4TH, or SYMPLECTIC_6TH"
#endif
#if DIS_SCHEME == LATTICEEASY
       write(*,*) "Using LatticeEasy discretization scheme"
#elif DIS_SCHEME == HLATTICE1
       write(*,*) "Using HLattice 1st discretization scheme"
#elif DIS_SCHEME == HLATTICE2
       write(*,*) "Using HLattice 2nd discretization scheme"
#else
       stop "Unknown DIS_SCHEME"
#endif
#if USE_STANDARD_WAVENUMBER
       write(*,*) "Using the standard wave vector"
#else
       write(*,*) "Using the effective wave vector"
#endif
#if METRIC_PERTURB
    write(*,*) "Metric perturbations in synchronous gauge are included."
#if DIS_SCHEME == LATTICEEASY
    stop "LATTICEEASY discretization cannot be used for metric perturbations"
#endif
    write(*,*) "Using physical time"
#else
    write(*,*) "Metric perturbations are ignored."
#if METRIC_OPTION == FRW_BACKGROUND
    write(*,*) "Will evolve fields in FRW background."
#elif METRIC_OPTION == MINKOWSKI_BACKGROUND
    write(*,*) "Will evolve fields in Minkowski background."
#endif    
#if USE_CONFORMAL_TIME
    write(*,*) "Using conformal time"
#else
    write(*,*) "Using physical time"
#endif
#endif
  end subroutine print_info 

end module init
