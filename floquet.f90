module global_variables
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical systems
  complex(8),allocatable :: zpsi(:,:) ! zpsi(2,nk)
  real(8),allocatable :: eps_bk(:,:)
  integer :: nk, nsym
  integer,allocatable :: nk_sym(:)
  real(8) :: a_vec(2,2), a_lattice, b_vec(2,2)
  real(8) :: delta_vec(2,3)
  real(8),allocatable :: kx0(:),ky0(:),kxt(:),kyt(:), kpath(:)
  real(8) :: t0_hop, eps_b, eps_n

! time propagation
  integer :: nt
  real(8) :: dt, Tperiod

! laser fields
  real(8),allocatable :: Act(:,:) 
  real(8) :: E0, omega0, Tpulse0

end module global_variables
!---------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input_variables
  call initialize

  call calc_ground_state

end program main
!---------------------------------------------------------------
subroutine input_variables
  use global_variables
  implicit none

! physical parameters
  eps_b = 3.34d0/27.2114d0
  eps_n = -2.56d0/27.2114d0
  t0_hop = 2.64d0/27.2114d0

! number of grid points
  nsym = 4
  allocate(nk_sym(nsym-1))
  nk_sym = 12
  nk = sum(nk_sym)


! lattice constant
  a_lattice=2.456d0/0.5291772d0  !!2.5AA

! lattice vectors
  a_vec(1,1) = a_lattice*sqrt(3d0)/2d0
  a_vec(2,1) = a_lattice*(0.5d0)

  a_vec(1,2) = a_lattice*sqrt(3d0)/2d0
  a_vec(2,2) = a_lattice*(-0.5d0)

! desplacement
  delta_vec(1,1) = a_lattice/sqrt(3d0)
  delta_vec(2,1) = a_lattice*0d0

  delta_vec(1,2) = a_lattice*(-1d0/(2d0*sqrt(3d0)))
  delta_vec(2,2) = a_lattice*0.5d0

  delta_vec(1,3) = a_lattice*(-1d0/(2d0*sqrt(3d0)))
  delta_vec(2,3) = a_lattice*(-0.5d0)


! time propagation
  omega0   = 0.95d0/27.2114d0  !ev
  Tperiod  = 2d0*pi/omega0
  nt = 32
  dt = Tperiod/nt


end subroutine input_variables
!---------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  real(8) :: volume, kvec_s(2), kvec_e(2), dk, ss, f1, f2
  integer :: ik1, ik2, ik, iks


  allocate(kx0(nk),ky0(nk),kxt(nk),kyt(nk),kpath(nk))
  allocate(zpsi(2,nk), eps_bk(2,nk))


  volume = (a_vec(1,1)*a_vec(2,2)-a_vec(2,1)*a_vec(1,2))*1d0
  
  b_vec(1,1) = 2d0*pi*(a_vec(2,2)*1d0)/volume
  b_vec(2,1) = 2d0*pi*(-a_vec(1,2)*1d0)/volume

  b_vec(1,2) = 2d0*pi*(-1d0*a_vec(2,1))/volume
  b_vec(2,2) = 2d0*pi*(1d0*a_vec(1,1))/volume

  write(*,*)'a-vec'
  write(*,*)a_vec(1,:)
  write(*,*)a_vec(2,:)

  write(*,*)'b-vec'
  write(*,*)b_vec(1,:)
  write(*,*)b_vec(2,:)


  write(*,*)sum(a_vec(:,1)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,1)*b_vec(:,2))/(2d0*pi)
  write(*,*)sum(a_vec(:,2)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,2)*b_vec(:,2))/(2d0*pi)

! K to Gamma
  kvec_s(1) = 2d0/3d0*b_vec(1,1) + 1d0/3d0*b_vec(1,2)
  kvec_s(2) = 2d0/3d0*b_vec(2,1) + 1d0/3d0*b_vec(2,2)
  kvec_e(1) = 0d0*b_vec(1,1) + 0d0*b_vec(1,2)
  kvec_e(2) = 0d0*b_vec(2,1) + 0d0*b_vec(2,2)
  dk = sqrt(sum( (kvec_e(:)-kvec_s(:))**2 ))/nk_sym(1)
  ss = 0d0
  ik = 0
  do iks = 1, nk_sym(1)
    ik = ik + 1
    kpath(ik) = ss + dk
    ss = kpath(ik)
    f1 = dble(iks-1)/nk_sym(1)
    f2 = 1d0 - f1
    kx0(ik) = kvec_s(1)*f2 + kvec_e(1)*f1
    ky0(ik) = kvec_s(2)*f2 + kvec_e(2)*f1
  end do

! Gamma to M
  kvec_s(1) = 0d0*b_vec(1,1) + 0d0*b_vec(1,2)
  kvec_s(2) = 0d0*b_vec(2,1) + 0d0*b_vec(2,2)
  kvec_e(1) = 0.5d0*b_vec(1,1) + 0.5d0*b_vec(1,2)
  kvec_e(2) = 0.5d0*b_vec(2,1) + 0.5d0*b_vec(2,2)
  dk = sqrt(sum( (kvec_e(:)-kvec_s(:))**2 ))/nk_sym(1)
  do iks = 1, nk_sym(2)
    ik = ik + 1
    kpath(ik) = ss + dk
    ss = kpath(ik)
    f1 = dble(iks-1)/nk_sym(2)
    f2 = 1d0 - f1
    kx0(ik) = kvec_s(1)*f2 + kvec_e(1)*f1
    ky0(ik) = kvec_s(2)*f2 + kvec_e(2)*f1
  end do

! M to K'
  kvec_s(1) = 0.5d0*b_vec(1,1) + 0.5d0*b_vec(1,2)
  kvec_s(2) = 0.5d0*b_vec(2,1) + 0.5d0*b_vec(2,2)
  kvec_e(1) = 1d0/3d0*b_vec(1,1) + 2d0/3d0*b_vec(1,2)
  kvec_e(2) = 1d0/3d0*b_vec(2,1) + 2d0/3d0*b_vec(2,2)
  dk = sqrt(sum( (kvec_e(:)-kvec_s(:))**2 ))/nk_sym(1)
  do iks = 1, nk_sym(3)
    ik = ik + 1
    kpath(ik) = ss + dk
    ss = kpath(ik)
    f1 = dble(iks-1)/nk_sym(3)
    f2 = 1d0 - f1
    kx0(ik) = kvec_s(1)*f2 + kvec_e(1)*f1
    ky0(ik) = kvec_s(2)*f2 + kvec_e(2)*f1
  end do


end subroutine initialize
!---------------------------------------------------------------
subroutine calc_ground_state
  use global_variables
  implicit none
  integer :: ik, ik1, ik2
  complex(8) :: zham(2,2), zvec(2,2), zfk
  real(8) :: eps_t(2), eps_bk_l(2,nk)
  real(8) :: kx_t, ky_t

  do ik = 1, nk
    
    kx_t = kx0(ik)
    ky_t = ky0(ik)

    zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

    zham(1,1) = eps_b
    zham(1,2) = t0_hop*zfk
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_n

    call calc_eig_vec_2x2(zham, zvec, eps_t)
    zpsi(:,ik) = zvec(:,1)
    eps_bk(:,ik) = eps_t(:)
    
  end do

  open(20,file='band_map.out')
  do ik = 1, nk
    write(20,"(999e26.16e3)")kpath(ik),kx0(ik),ky0(ik),eps_bk(1:2,ik)
  end do
  close(20)

end subroutine calc_ground_state
!---------------------------------------------------------------  
subroutine calc_eig_vec_2x2(zham, zvec, eps_t)
  implicit none
  complex(8),intent(in) :: zham(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out) :: eps_t(2)
  real(8) :: a,c
  complex(8):: zb, zx, zy

! (a    zb)
! (zb^*  c)

  a = zham(1,1)
  zb = zham(1,2)
  c = zham(2,2)

  eps_t(1) = 0.5d0*( (a+c)-sqrt((a-c)**2+4d0*abs(zb)**2))
  eps_t(2) = 0.5d0*( (a+c)+sqrt((a-c)**2+4d0*abs(zb)**2))

  if(a<c)then
    zy = conjg(zb)/(eps_t(1)-c)
    zvec(1,1) = 1d0/sqrt(1d0+abs(zy)**2)
    zvec(2,1) = zy/sqrt(1d0+abs(zy)**2)

    zx = zb/(eps_t(2)-a)
    zvec(1,2) = zx/sqrt(1d0+abs(zx)**2)
    zvec(2,2) = 1d0/sqrt(1d0+abs(zx)**2)

  else

    zx = zb/(eps_t(1)-a)
    zvec(1,1) = zx/sqrt(1d0+abs(zx)**2)
    zvec(2,1) = 1d0/sqrt(1d0+abs(zx)**2)

    zy = conjg(zb)/(eps_t(2)-c)
    zvec(1,2) = 1d0/sqrt(1d0+abs(zy)**2)
    zvec(2,2) = zy/sqrt(1d0+abs(zy)**2)

  end if



end subroutine calc_eig_vec_2x2
!---------------------------------------------------------------
subroutine init_laser
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, xx


  allocate(Act(2,-1:nt+1))
  Act = 0d0

  do it = 0, nt
    tt = dt*it
    xx = tt -0.5d0*Tpulse0

    if(abs(xx)<= 0.5d0*Tpulse0)then
      Act(1,it) = - E0/omega0*cos(omega0*xx)*cos(pi*xx/Tpulse0)**4
      Act(2,it) = - E0/omega0*sin(omega0*xx)*cos(pi*xx/Tpulse0)**4
    end if

  end do

  open(20,file='laser.out')
  do it = 0, nt
    tt = dt*it
    write(20,"(999e26.16)")tt,Act(:,it)
  end do
  close(20)

end subroutine init_laser
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
