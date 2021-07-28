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

! Floquet calculation
  integer,parameter :: nmax_floquet = 72
  integer,parameter :: ndim_F = 2*(2*nmax_floquet+1)
  complex(8) :: zham_floquet(2,2,-2*nmax_floquet:2*nmax_floquet)

end module global_variables
!---------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input_variables
  call initialize

  call calc_ground_state
  call calc_floquet_state

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
  nk_sym = 64
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
  omega0   = 0.25d0/27.2114d0  !ev
  Tperiod  = 2d0*pi/omega0
  nt = 512
  dt = Tperiod/nt

  E0 = 0.01d0

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
  real(8) :: Epdir(2), ss

! Gamma-K direction
  Epdir(1) = 2d0/3d0*b_vec(1,1) + 1d0/3d0*b_vec(1,2)
  Epdir(2) = 2d0/3d0*b_vec(2,1) + 1d0/3d0*b_vec(2,2)
  ss = sqrt(sum((Epdir)**2))
  Epdir = Epdir/ss

  allocate(Act(2,0:nt-1))
  Act = 0d0

  do it = 0, nt-1
    tt = dt*it
    Act(1,it) = -Epdir(1)*E0/omega0*cos(omega0*tt)
    Act(2,it) = -Epdir(2)*E0/omega0*cos(omega0*tt)
  end do

  open(20,file='laser.out')
  do it = 0, nt-1
    tt = dt*it
    write(20,"(999e26.16)")tt,Act(:,it)
  end do
  close(20)

end subroutine init_laser
!---------------------------------------------------------------
subroutine calc_floquet_state
  use global_variables
  implicit none
  integer :: ik, ii
  real(8) :: eps_F(ndim_F), occ_F(ndim_F)
  call init_laser

  open(40,file='floquet_band.out')
  do ik = 1, nk
    
    call calc_floquet_sub_matrix(ik)
    call diagonalize_floqet_matrix(eps_F, occ_F)

    write(40,"(999e26.16e3)")kpath(ik),(eps_F(ii), occ_F(ii), ii=1,ndim_F)
    
  end do
  close(40)

end subroutine calc_floquet_state
!---------------------------------------------------------------
subroutine calc_floquet_sub_matrix(ik)
  use global_variables
  implicit none
  integer,intent(in) :: ik
  complex(8) :: zham(2,2), zfk
  integer :: ii, it
  real(8) :: kx_t, ky_t

  zham_floquet = 0d0


  do ii = -2*nmax_floquet, 2*nmax_floquet
    
    do it = 0, nt-1
      kx_t = kx0(ik) + Act(1,it)
      ky_t = ky0(ik) + Act(2,it)

      zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
          + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
          + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

      zham(1,1) = eps_b
      zham(1,2) = t0_hop*zfk
      zham(2,1) = conjg(zham(1,2))
      zham(2,2) = eps_n

      zham_floquet(:,:,ii) = zham_floquet(:,:,ii) + zham(:,:)*exp(zi*ii*omega0*dt*it)
    end do
  end do

  zham_floquet = zham_floquet/nt

!! debug
!  zham_floquet = 0d0
!  kx_t = kx0(ik)
!  ky_t = ky0(ik)
!
!  zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
!      + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
!      + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 
!
!  zham(1,1) = eps_b
!  zham(1,2) = t0_hop*zfk
!  zham(2,1) = conjg(zham(1,2))
!  zham(2,2) = eps_n
!  zham_floquet(:,:,0) = zham

end subroutine calc_floquet_sub_matrix
!---------------------------------------------------------------
subroutine diagonalize_floqet_matrix(eps_F,occ_F)
  use global_variables
  implicit none
  real(8),intent(out) :: eps_F(ndim_F), occ_F(ndim_F)
  complex(8),allocatable :: zham_F(:,:)
  integer :: idim, jdim, ii
!==LAPACK
  integer :: nmax
  integer :: lwork
  complex(8),allocatable :: za(:,:),work_lp(:)
  real(8),allocatable :: rwork(:),w(:)
  integer :: info
!==LAPACK

  allocate(zham_F(ndim_F, ndim_F))

!==LAPACK
  nmax = ndim_F
  lwork = 6*nmax**2
  allocate(za(nmax,nmax),work_lp(lwork),rwork(3*nmax-2),w(nmax))
!==LAPACK

  zham_F = 0d0

  do idim = 1, 2*nmax_floquet + 1
    do jdim = 1, 2*nmax_floquet + 1

      ii = idim-jdim
      zham_F(2*(idim-1)+1:2*(idim-1)+2,2*(jdim-1)+1:2*(jdim-1)+2) &
        = zham_floquet(1:2,1:2,ii)
    end do

    zham_F(2*(idim-1)+1,2*(idim-1)+1) &
              = zham_F(2*(idim-1)+1,2*(idim-1)+1) &
              - omega0*(idim-nmax_floquet -1)

    zham_F(2*(idim-1)+2,2*(idim-1)+2) &
              = zham_F(2*(idim-1)+2,2*(idim-1)+2) &
              - omega0*(idim-nmax_floquet -1)
  end do

  za = zham_F
  Call zheev('V', 'U', nmax, za, nmax, w, work_lp, lwork, rwork, info)

  eps_F = w
  idim = ndim_F/2 !2*(nmax_floquet + 1 -1)+1
!  idim = 2*(nmax_floquet + 1 -1)+1
  do ii = 1, ndim_F
    occ_F(ii) = abs(za(idim,ii))**2+abs(za(idim+1,ii))**2
  end do

end subroutine diagonalize_floqet_matrix
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
!---------------------------------------------------------------
