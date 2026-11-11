! Forward 1-D Rayleigh dispersion using exactly the DSurfTomo v1.4 chain:
! Brocher (2005) Vp/rho from Vs, refineGrid2LayerMdl, surfdisp96.
! Built by lant install-external against the unmodified DSurfTomo sources (external/DSurfTomo/src).
!
!   fwd1d model.txt periods.txt igr minthk [iwave]
!
! model.txt   : first line nz, then nz lines "depth(km) vs(km/s)"
! periods.txt : first line kmax, then kmax periods (s)
! igr         : 0 = phase velocity, 1 = group velocity
! minthk      : sublayer thickness divisor used by refineGrid2LayerMdl
! iwave       : optional, 2 = Rayleigh (default), 1 = Love (surfdisp96's codes, as DSurfTomo's wave type)
!
! writes "period velocity" to stdout, one line per period.
program fwd1d
  implicit none
  integer, parameter :: NL = 200, NP = 80   ! NP: the limit of surfdisp96
  integer :: nz, kmax, i, iflsph, iwave, mode, igr, rmax, ios
  real :: minthk
  real, allocatable :: depz(:), vsz(:), vpz(:), rhoz(:)
  real*8, allocatable :: t(:)
  real :: rdep(NL), rvp(NL), rvs(NL), rrho(NL), rthk(NL)
  real*8 :: tRc(NP), cgRc(NP)
  character(len=256) :: mfile, pfile, arg

  call get_command_argument(1, mfile)
  call get_command_argument(2, pfile)
  call get_command_argument(3, arg); read(arg, *) igr
  call get_command_argument(4, arg); read(arg, *) minthk
  iwave = 2
  if (command_argument_count() >= 5) then
     call get_command_argument(5, arg); read(arg, *) iwave
  end if

  open(10, file=trim(mfile), status='old', action='read')
  read(10, *) nz
  allocate(depz(nz), vsz(nz), vpz(nz), rhoz(nz))
  do i = 1, nz
     read(10, *) depz(i), vsz(i)
  end do
  close(10)

  open(11, file=trim(pfile), status='old', action='read')
  read(11, *) kmax
  if (kmax > NP) stop 'too many periods (max 80, the limit of surfdisp96)'
  allocate(t(kmax))
  do i = 1, kmax
     read(11, *) t(i)
  end do
  close(11)
  tRc(1:kmax) = t(1:kmax)

  ! Brocher (2005), identical to CalSurfG.f90 lines 49-53
  do i = 1, nz
     vpz(i) = 0.9409 + 2.0947*vsz(i) - 0.8206*vsz(i)**2 + &
              0.2683*vsz(i)**3 - 0.0251*vsz(i)**4
     rhoz(i) = 1.6612*vpz(i) - 0.4721*vpz(i)**2 + &
               0.0671*vpz(i)**3 - 0.0043*vpz(i)**4 + &
               0.000106*vpz(i)**5
  end do

  iflsph = 1
  mode = 1
  call refineGrid2LayerMdl(minthk, nz, depz, vpz, vsz, rhoz, rmax, rdep, &
                           rvp, rvs, rrho, rthk)
  call surfdisp96(rthk, rvp, rvs, rrho, rmax, iflsph, iwave, mode, igr, &
                  kmax, tRc, cgRc)
  do i = 1, kmax
     write(*, '(F10.4,1X,F12.6)') tRc(i), cgRc(i)
  end do
end program fwd1d
