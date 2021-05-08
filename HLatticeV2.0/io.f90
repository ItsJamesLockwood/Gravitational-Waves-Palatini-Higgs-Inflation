 module io_utils
   use model
   use mutils
   implicit none
#include "configure.h"
   Type(File_Pointer)::gw_file,screen_file
contains

  subroutine write_check()
    type(file_pointer) chkfile
    DEFINE_IND
    write(*,*) "Writing checkpoint files."
    chkfile = open_file("data/"//trim(run_name)//"_chk.log","u")
    write(chkfile%unit) n
    write(chkfile%unit) ns
    write(chkfile%unit) METRIC_OPTION
    write(chkfile%unit) DIS_SCHEME
    write(chkfile%unit) USE_CONFORMAL_TIME
    write(chkfile%unit) WANTGW
    write(chkfile%unit) sipar%nsteps
    write(chkfile%unit) sipar%time
    call close_file(chkfile)
    call fields_dump()
    call output_pw()
    call metric_dump()

#if ! METRIC_PERTURB && WANTGW
#if USE_CONFORMAL_TIME
    metric_h = metric_h + metric_p * (scale_factor()*sipar%dt/2.*n_feedback)
#else
    metric_h = metric_h + metric_p * (sipar%dt/2.*n_feedback)
#endif
#endif
    call output_GW()

#if ! METRIC_PERTURB && WANTGW
#if USE_CONFORMAL_TIME
    metric_h = metric_h - metric_p * (scale_factor()*sipar%dt/2.*n_feedback)
#else
    metric_h = metric_h - metric_p * (sipar%dt/2.*n_feedback)
#endif
#endif
  end subroutine write_check

  subroutine read_check(ierror)
    type(file_pointer) chkfile
    logical ierror
    integer(IB) tmp
    If(FileExists("data/"//trim(run_name)//"_fields.log") .and. FileExists("data/"//trim(run_name)//"_metric_homo.log"))then
       write(*,*) "Reading checkpoint files."
       sipar%dt = model_dt_per_step()
       chkfile = open_file("data/"//trim(run_name)//"_chk.log","u")
       read(chkfile%unit) tmp
       if(n.ne.tmp) stop "Loading checkpoint files failed: SIMU_RESOLUTION do not agree."
       read(chkfile%unit) tmp
       if(ns.ne.tmp) stop "Loading checkpoint files failed: NUM_SCALAR_FIELDS do not agree."
       read(chkfile%unit) tmp 
       if(tmp .ne. METRIC_OPTION) stop "Loading checkpoint files failed: METRIC_OPTION do not agree."
       read(chkfile%unit)tmp
       if(tmp .ne. DIS_SCHEME) stop "Loading checkpoint files failed: DIS_SCHEME do not agree."
       read(chkfile%unit)tmp
       if(tmp .ne.  USE_CONFORMAL_TIME) stop "Loading checkpoint files failed:  USE_CONFORMAL_TIME do not agree."
       read(chkfile%unit)tmp
       if(tmp .ne.  WANTGW) stop "Loading checkpoint files failed: WANTGW_TIME do not agree."
       call close_file(chkfile)
       call fields_load()
       call metric_load()
       ierror=.false.
    else
       write(*,*) "Can not find checkpoint files. Doing normal initialization."
       ierror =.true.
    endif
    return
  end subroutine read_check

  subroutine begin_output()
    integer(IB) i
    call prtsystime(.true.)
  end subroutine begin_output

  subroutine output_to_screen()
    DEFINE_IND
    character (len=20):: formatString
    character (len=100):: totalFormat
    character (len=10):: parentheses
    character (len=20):: lastTerm
    real(dl) avef(ns),flucf(ns)
    formatString =',4x,G16.7'
    totalFormat ='(F16.7'
    lastTerm=',4x,E16.7'
    parentheses=')'
    screen_file = open_file("data/"//trim(run_name)//"_screen.log","a")
    call get_energy()
    do i=1,ns
       avef(i)=sum(fields_f(i,:,:,:))/ncube
    enddo
    do i=1,ns
       flucf(i)=sqrt(sum((fields_f(i,:,:,:)-avef(i))**2)/ncube)
    enddo
    do i=1, 2*ns+4
      totalFormat=trim(totalFormat)//trim(formatString)
    enddo
    totalFormat=trim(totalFormat)//trim(lastTerm)//trim(parentheses)
    avef=avef/Mpl
    flucf=flucf/Mpl
#if FEEDBACK_ONLYAH
    if(.not. (use_checkpoint .and. mod(sipar%nsteps,checkpoint_steps).eq.0 .and. (sipar%nsteps.gt.0)))then
       write(*,'(2G16.5)') metric%a, effective_Hubble()
       write(screen_file%unit,'(2G16.5)') metric%a, effective_Hubble()
       goto 100
    endif
#endif
#if METRIC_PERTURB
    if(sipar%nsteps.eq.0)then
       write(*,'(A9, 9A11, 2A'//trim(int2str(ns*11))//')') 'a   ','H  ','E_f ','P_f/E_f', 'K_f/E_f', 'G_f/E_f', 'K_g/E_f', 'G_g/E_f', 'E_tot/E_f', 'rms_h ', 'mean_fields',  'rms_fields'
       write(screen_file%unit,'(A9, 9A11, 2A'//trim(int2str(ns*11))//')') 'a  ','H  ','E_f  ','P_f/E_f', 'K_f/E_f', 'G_f/E_f', 'K_g/E_f', 'G_g/E_f', 'E_tot/E_f', 'rms_h ', 'mean_fields',  'rms_fields'
    endif
    write(*,'(F9.5,'//trim(Int2str(2*ns+9))//'G11.3)') metric%a, effective_Hubble(), total_fields_energy(), potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), gravity_kinetic_energy()/total_fields_energy(),gravity_gradient_energy()/total_fields_energy(), total_energy()/total_fields_energy(),sqrt(sum(metric_h(:,1:n,1:n,1:n)**2)/ncube/6._dl), avef, flucf
    write(screen_file%unit,'(F10.5,'//trim(Int2str(2*ns+9))//'G11.3)') metric%a, effective_Hubble(), total_fields_energy(), potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), gravity_kinetic_energy()/total_fields_energy(),gravity_gradient_energy()/total_fields_energy(), total_energy()/total_fields_energy(),sqrt(sum(metric_h(:,1:n,1:n,1:n)**2)/ncube/6._dl), avef, flucf
#else
    if(sipar%nsteps.eq.0)then
       write(*,'(A9, 5A11, 2A'//trim(int2str(ns*11))//')') 'a  ','H  ','rho/3H^2-1','P_f/E_f   ', 'K_f/E_f   ', 'G_f/E_f   ', 'mean_fields ',  'rms_fields'
       write(screen_file%unit,'(A9, 5A11, 2A'//trim(int2str(ns*11))//')') 'a  ','H  ','rho/3H^2 -1','P_f/E_f   ', 'K_f/E_f   ', 'G_f/E_f   ', 'mean_fields ',  'rms_fields'
    endif
    write(*,(totalFormat)) metric%a,effective_Hubble(),total_fields_energy()/effective_Hubble()**2/3/Mplsq-1._dl, potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), avef,flucf
    !!write(*,*) "Time:", sipar%time , sipar%dt * n_feedback
    write(screen_file%unit,(totalFormat)) metric%a,effective_Hubble(),total_fields_energy()/effective_Hubble()**2/3/Mplsq-1._dl, potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), avef,flucf
    !!!write(*,'(F10.5,'//trim(Int2str(2*ns+5))//'G11.3)')metric%a,effective_Hubble(),total_fields_energy()/effective_Hubble()**2/3/Mplsq-1._dl, potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), avef,flucf
    !! Use the line below only during troubleshooting of y, P_f and piy.
    !write(*,*) "metric%y:",metric%y, "; P_f:", potential_energy(),"; piy:",metric%piy
    !!!write(screen_file%unit,'(F10.5,'//trim(Int2str(2*ns+5))//'G11.3)')metric%a,effective_Hubble(),total_fields_energy()/effective_Hubble()**2/3/Mplsq-1._dl, potential_energy()/total_fields_energy(), fields_kinetic_energy()/total_fields_energy(), fields_gradient_energy()/total_fields_energy(), avef,flucf
#endif
100    call model_output(3)
    call close_file(screen_file)
  end subroutine output_to_screen

  subroutine final_output()
    call print_bar()
    write(*,*)"Simulation stops at t =",sipar%time
    call model_output(4)
    call print_bar()
  end subroutine final_output
  
  subroutine output_GW()
    real(dl) tote,pk(fft_numk),freq(fft_numk)
    logical ierror
    type(file_pointer) fph
    integer(IB) i
    character(LEN=128)::fmt
#if HIJ_DEFINED
    gw_file = open_file("data/"//trim(run_name)//"_GW.log", "a")
    tote=potential_energy()+fields_kinetic_energy()+fields_gradient_energy()
    fmt='('//trim(int2str(fft_numk))//'G16.7)'
    call get_GW_spectrum(pk)
    do i=1,fft_numk
       freq(i)=4.e10 * i * k_unit / metric%physdx / (tote)**0.25
    enddo
    write(gw_file%unit,'(G16.7)') metric%a
    write(gw_file%unit,trim(fmt)) freq
    write(gw_file%unit,trim(fmt)) pk/ tote * 9.3e-6_dl
    call close_file(gw_file)
#endif
  end subroutine output_GW


  subroutine output_pw()
    real(dl) fpk(fft_numk),pipk(fft_numk)
    integer(IB) i,j
    type(file_pointer) pw_file
    character(LEN=128) fmt
    real(dl) tote,meanf(ns),meanp(ns),ge
    fmt='('//trim(int2str(fft_numk))//'E16.7)'
    ge=fields_gradient_energy()
    tote=potential_energy()+fields_kinetic_energy() + ge
    if(ge/tote .lt. 1.e-4)then
       do i=1,ns
          meanf(i) = sum(fields_f(i,:,:,:))/ncube
          fields_f(i,:,:,:)= fields_f(i,:,:,:) - meanf(i)
          meanp(i) = sum(fields_p(i,:,:,:))/ncube
          fields_p(i,:,:,:)= fields_p(i,:,:,:) - meanp(i)
       enddo
    endif
    fields_p=fields_p*(metric%physdx/metric%a**3)
    call CubicMassiveFFT(ns,fields_f,fields_p, FFT_FORWARD)
    do i=1,ns
       pw_file = open_file("data/"//trim(run_name)//"_pw_"//trim(int2str(i))//".log","a")
       write(pw_file%unit,'(G16.7)') metric%a
       fft_ind_re = i
       fft_ind_im = i
       call fft_getPowerSpectrum(fields_f,fields_p,fpk, pipk)
       do j=1,fft_numk
          fpk(j) = fpk(j)*fft_vol(j)*real(j,dl)**3
          pipk(j) = pipk(j)*fft_vol(j)*real(j,dl)
       enddo
       write(pw_file%unit,trim(fmt)) fpk/(2.*tote*ncube**2)*(const_2pi/n/metric%physdx)**2
       write(pw_file%unit,trim(fmt)) pipk/(2.*tote*ncube**2)/metric%physdx**2
       call close_file(pw_file)
    enddo
    call CubicMassiveFFT(ns,fields_f,fields_p, FFT_BACKWARD)
    fields_p=fields_p/(metric%physdx/metric%a**3)
    if(ge/tote .lt. 1.e-4)then
       do i=1,ns
          fields_f(i,:,:,:)= fields_f(i,:,:,:) + meanf(i)
          fields_p(i,:,:,:)= fields_p(i,:,:,:) + meanp(i)
       enddo
    endif
  end subroutine output_pw

  subroutine print_settings_file()
    type(file_pointer) settings_file
    settings_file = open_file("data/"//trim(run_name)//"_sim_settings.log", "w")
    write(settings_file%unit,*) "Number of fields: ", ns
    write(settings_file%unit,*) "Initial boxsize times H: ", boxsize_H
    write(settings_file%unit,*) "Discretisatin scheme [1-3]: ", DIS_SCHEME 
    write(settings_file%unit,*) "GW turned on? [0/1]: ", WANTGW 
    write(settings_file%unit,*) "Metric [1-4]: ", METRIC_OPTION
    write(settings_file%unit,*) "Integrator [1-3]: ", INTEGRATOR
    write(settings_file%unit,*) "Lattice points on side: ", n
    write(settings_file%unit,*) "Conformal time used [0/1]: ", USE_CONFORMAL_TIME
    write(settings_file%unit,*) "Max scale factor: ", stop_at_a
    write(settings_file%unit,*) "Max steps: ", MAX_STEPS
    write(settings_file%unit,*) "Feedback interval (for evolving h_ij): ", FEEDBACK_INTERVAL
    write(settings_file%unit,*) "Checkpoints (for gauge scaling): ", CHECKPOINT_AND_GW_INTERVAL
    write(settings_file%unit,*) "\n<Non-standard settings>"
    write(settings_file%unit,*) "Sub integrator steps: ", SUB_INTEGRATOR_STEPS
    write(settings_file%unit,*) "Reduced Planck mass: ", REDUCED_PLANCK_MASS
    write(settings_file%unit,*) "Use standard wavenumber (i.e. not k_eff): ",USE_STANDARD_WAVENUMBER
    write(settings_file%unit,*) "Only print a and H [0/1]: ", FEEDBACK_ONLYAH 
    call close_file(settings_file)
  end subroutine print_settings_file

   subroutine output_fields()
      type(file_pointer) fp
      integer(IB) nf,i,j
      do nf=1, ns
         fp = open_file("data/" // trim(run_name) // "_whole_field_"//trim(int2str(nf))//".log","a")
         if (sipar%nsteps.eq.0) then
            write(fp%unit,*) n
         endif
         write(fp%unit,*) "SEPARATOR"
         write(fp%unit,*) metric%a
         ! do i=1,n 
         !    do j=1,n
         !       write(fp%unit,*) i,j
         !       write(fp%unit,*) fields_f(nf,i,j,:)               
         !    enddo
         ! enddo
         write(fp%unit,*) fields_f(nf,:,:,:)
         call close_file(fp)
      enddo
      
   end subroutine output_fields

   subroutine output_momenta()
      type(file_pointer) fp
      integer(IB) nf,i,j
      do nf=1, ns
         fp = open_file("data/" // trim(run_name) // "_momenta_"//trim(int2str(nf))//".log","a")
         if (sipar%nsteps.eq.0) then
            write(fp%unit,*) n
         endif
         write(fp%unit,*) "SEPARATOR"
         write(fp%unit,*) metric%a
         ! do i=1,n 
         !    do j=1,n
         !       write(fp%unit,*) i,j
         !       write(fp%unit,*) fields_f(nf,i,j,:)               
         !    enddo
         ! enddo
         write(fp%unit,*) fields_p(nf,:,:,:)
         call close_file(fp)
      enddo
      
   end subroutine output_momenta

   subroutine output_field_slice()
      type(file_pointer) fp
      integer(IB) nf,i,j
      do nf=1, ns
         fp = open_file("data/" // trim(run_name) // "_slice_f_"//trim(int2str(nf))//".log","a")
         if (sipar%nsteps.eq.0) then
            write(fp%unit,*) n
         endif
         write(fp%unit,*) "SEPARATOR"
         write(fp%unit,*) metric%a
         ! do i=1,n 
         !    do j=1,n
         !       write(fp%unit,*) i,j
         !       write(fp%unit,*) fields_f(nf,i,j,:)               
         !    enddo
         ! enddo
         if (which_lattice_cut .eq. 1) then
            write(fp%unit,*) fields_f(nf,yz_lattice_slice,:,:)
         else if (which_lattice_cut .eq. 2) then
            write(fp%unit,*) fields_f(nf,:,:,xy_lattice_slice)
         else if (which_lattice_cut .eq. 3) then
            write(fp%unit,*) fields_f(nf,:,xz_lattice_slice,:)
         else 
            write(*,*) "Skipping field slices: slice direction not recognised"
         endif
         call close_file(fp)
      enddo
      
   end subroutine output_field_slice

   subroutine output_momentum_slice()
      type(file_pointer) fp
      integer(IB) nf,i,j
      do nf=1, ns
         fp = open_file("data/" // trim(run_name) // "_slice_p_"//trim(int2str(nf))//".log","a")
         if (sipar%nsteps.eq.0) then
            write(fp%unit,*) n
         endif
         write(fp%unit,*) "SEPARATOR"
         write(fp%unit,*) metric%a
         ! do i=1,n 
         !    do j=1,n
         !       write(fp%unit,*) i,j
         !       write(fp%unit,*) fields_f(nf,i,j,:)               
         !    enddo
         ! enddo
         if (which_lattice_cut .eq. 1) then
            write(fp%unit,*) fields_p(nf,yz_lattice_slice,:,:)
         else if (which_lattice_cut .eq. 2) then
            write(fp%unit,*) fields_p(nf,:,:,xy_lattice_slice)
         else if (which_lattice_cut .eq. 3) then
            write(fp%unit,*) fields_p(nf,:,xz_lattice_slice,:)
         else 
            write(*,*) "Skipping momentum slices: slice direction not recognised"
         endif
         call close_file(fp)
      enddo
      
   end subroutine output_momentum_slice

   subroutine output_energy()
      type(file_pointer) fp
      integer(IB) nf,i,j
      do nf=1, ns
         fp = open_file("data/" // trim(run_name) // "_FPE_"//trim(int2str(nf))//".log","a")
         if (sipar%nsteps.eq.0) then
            write(fp%unit,*) n
         endif
         write(fp%unit,*) "SEPARATOR"
         write(fp%unit,*) metric%a
         write(fp%unit,*) FPE_mesh
         write(fp%unit,*) FKE_mesh
         write(fp%unit,*) FGE_mesh
         call close_file(fp)
      enddo  
   end subroutine output_energy

end module io_utils
