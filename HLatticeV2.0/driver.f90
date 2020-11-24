program GHLattice
  use init
  use evolution
  use io_utils
  implicit none
  run_name = Trim(GetParam(1))
  if(Trim(run_name) .eq. "") run_name = "test"
  call initialize()
  call evolve()
end program GHLattice
