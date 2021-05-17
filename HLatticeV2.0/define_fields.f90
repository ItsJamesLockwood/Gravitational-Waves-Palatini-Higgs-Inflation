module define_fields
  use parameters
  use utils
  implicit none

#include "configure.h"
  
  type FRWmetric
     !!define the metric in synchronous gauge
     real(dl) a !!scale factor, a function of time
     real(dl) dx !!grid size, comoving dx
     real(dl) physdx !!physical dx
     real(dl) y,piy
 !!if using physical time: y = a^(3/2), piy = -8/3 dy/dt
 !!if using conformal time: y=a, piy = -6 dy/dt
  end type FRWmetric

  type simulation_parameters
     integer(IB) nsteps
     real(dl) time
     real(dl) dt,boxsize
     !!fields/gravity kinetic energy, gradient energy, and fields potential energy
     real(dl) fke,gke,fge,gge,pe
     !!the cach to temporarily store calculated metric g^{ij} and the derivatives of metric on three slices on x-y plane (in z direction indices=k-1,k,k+1)
#if METRIC_PERTURB
#if DIS_SCHEME == HLATTICE1
     real(dl) cach_dg(n,n,0:2)  !!\sqrt(g)
     real(dl) cach_dgup(6,n,n,0:2)  !!\sqrt(g)g^{ij}
     real(dl) cach_ddgup(6,6,n,n,0:2) !!\frac{\paritial \sqrt(g) g^{ij}}{\partial h_{ij}}
     integer cach_ind(-1:1)
#elif DIS_SCHEME == HLATTICE2
     real(dl) cach_dg(n,n,0:4)  !!\sqrt(g)
     real(dl) cach_dgup(6,n,n,0:4)  !!\sqrt(g)g^{ij}
     real(dl) cach_ddgup(6,6,n,n,0:4) !!\frac{\paritial \sqrt(g) g^{ij}}{\partial h_{ij}}
     integer cach_ind(-2:2)
#endif
#endif
  end type simulation_parameters


  Type(FRWmetric) metric
  type(simulation_parameters)::sipar
  real(dl),dimension(:,:,:,:),allocatable::fields_f
  real(dl),dimension(:,:,:,:),allocatable::fields_p
  real(dl),dimension(:,:,:,:),allocatable::metric_h
  real(dl),dimension(:,:,:,:),allocatable::metric_p
#if METRIC_OPTION == FRW_PERTURB_ADAPTIVE
  real(dl),dimension(:,:,:,:),allocatable::gauge_lambda
#endif

   real(dl),dimension(:,:,:),allocatable::pressure_field
   real(dl),dimension(:,:,:),allocatable::density_field
contains

!! input/output utils   
   subroutine fields_dump() 
     type(file_pointer) fp
     DEFINE_IND
     fp = open_file("data/"//trim(run_name)//"_fields.log","u")
     write(fp%unit) n
     write(fp%unit) ns
     do i=1,n
        write(fp%unit) fields_f(:,:,:,i)
        write(fp%unit) fields_p(:,:,:,i)
     enddo
     call close_file(fp)
   end subroutine fields_dump

   subroutine fields_load()
     type(file_pointer) fp
     DEFINE_IND
     fp = open_file("data/"//trim(run_name)//"_fields.log","u")
     read(fp%unit) i
     read(fp%unit) j
     if(i.ne.n .or. j.ne.ns) goto 100
     do i=1,n
        read(fp%unit) fields_f(:,:,:,i)
        read(fp%unit) fields_p(:,:,:,i)        
     enddo
     call close_file(fp)
     write(*,*) "Fields configuration loaded."
     return
100  write(*,*) "Configuration does not match. Can not load the fields."
     stop
   end subroutine fields_load

   subroutine metric_dump_term(filename,term)
     character(LEN=*) filename
     type(file_pointer) fp
     integer(IB)::term
     DEFINE_IND
     fp = open_file(filename, 'u')
     select case(term)
     case(0)
        write(fp%unit) n
        write(fp%unit) metric%a
        write(fp%unit) metric%dx
        write(fp%unit) metric%physdx
        write(fp%unit) metric%y
        write(fp%unit) metric%piy
     case(1)
        write(fp%unit) n
        do i=1,n
           write(fp%unit) metric_h(:,:,:,i)
        enddo
     case(2)
        write(fp%unit) n
        do i=1,n
           write(fp%unit) metric_p(:,:,:,i)
        enddo
     end select
     call close_file(fp)
     return
   end subroutine metric_dump_term

   subroutine metric_load_term(filename,term)
     character(LEN=*) filename
     type(file_pointer) fp
     integer(IB) term
     DEFINE_IND
     fp = open_file(filename, 'u')
     read(fp%unit) i
     if(i.ne. n) goto 100
     select case(term)
     case(0)
        read(fp%unit) metric%a
        read(fp%unit) metric%dx
        read(fp%unit) metric%physdx
        read(fp%unit,END=100) metric%y
        read(fp%unit,END=100) metric%piy
     case(1)
        do i=1,n
           read(fp%unit) metric_h(:,:,:,i)
        enddo
     case(2)
        do i=1,n
           read(fp%unit) metric_p(:,:,:,i)
        enddo
     end select
     call close_file(fp)
     return
100  write(*,*) "Configuration does not match. Can not load the fields."
     stop
   end subroutine metric_load_term


   subroutine metric_dump()
     call metric_dump_term("data/"//trim(run_name)//"_metric_homo.log",0)
#if HIJ_DEFINED
     call metric_dump_term("data/"//trim(run_name)//"_metric_hij.log",1)
     call metric_dump_term("data/"//trim(run_name)//"_metric_pij.log",2)
#endif
   end subroutine metric_dump

   subroutine metric_load()
     call metric_load_term("data/"//trim(run_name)//"_metric_homo.log",0)
#if HIJ_DEFINED
     call metric_load_term("data/"//trim(run_name)//"_metric_hij.log",1)
     call metric_load_term("data/"//trim(run_name)//"_metric_pij.log",2)
     write(*,*) "metric loaded"
#endif
   end subroutine metric_load


   subroutine print_fields(f)
     real(dl) f(ns)
     character(LEN=32) fmt
     fmt = "("//trim(int2str(ns))//"G12.6)"
     write(*,fmt) f
   end subroutine print_fields

   subroutine deallocate_fields()
     if(allocated(fields_f)) deallocate(fields_f)
     if(allocated(fields_p)) deallocate(fields_p)
   end subroutine deallocate_fields

   subroutine allocate_fields()
     call deallocate_fields()
     allocate(fields_f(ns,n,n,n))
     allocate(fields_p(ns,n,n,n))
   end subroutine allocate_fields

   subroutine deallocate_metric_h()
     if(allocated(metric_h)) deallocate(metric_h)
   end subroutine deallocate_metric_h

   subroutine deallocate_metric_p()
     if(allocated(metric_p)) deallocate(metric_p)
   end subroutine deallocate_metric_p

   subroutine deallocate_metric()
     call deallocate_metric_h()
     call deallocate_metric_p()
   end subroutine deallocate_metric

   subroutine allocate_metric_h()
#if HIJ_DEFINED
     allocate(metric_h(6,n,n,n))
#endif
   end subroutine allocate_metric_h

   subroutine allocate_metric_p()
#if HIJ_DEFINED
     allocate(metric_p(6,n,n,n))
#endif
   end subroutine allocate_metric_p


   subroutine allocate_metric()
     call deallocate_metric()
#if HIJ_DEFINED
     call allocate_metric_h()
     call allocate_metric_p()
#endif     
   end subroutine allocate_metric

!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!!! Modifications to source code %%%%
!!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   subroutine deallocate_density_fields()
      if(allocated(density_field)) deallocate(density_field)
      if(allocated(pressure_field)) deallocate(pressure_field)
   end subroutine deallocate_density_fields

   subroutine allocate_density_fields()
      call deallocate_density_fields()
      allocate(density_field(n,n,n))
      allocate(pressure_field(n,n,n))
      !!density_field = 0._dl
      !!pressure_field = 0._dl
   end subroutine allocate_density_fields


end module define_fields
