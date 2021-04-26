module global_variables
  implicit none
! mathematical parameters
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical systems
  complex(8),allocatable :: zpsi(:,:) ! zpsi(2,nk)
  real(8),allocatable :: eps_bk(:,:)
  integer :: nk, nk1, nk2
  real(8) :: a_vec(2,2), a_lattice, b_vec(2,2)
  real(8) :: delta_vec(2,3)
  real(8),allocatable :: kx0(:),ky0(:),kxt(:),kyt(:)
  real(8) :: t0_hop, eps_b, eps_n

! time propagation
  integer :: nt
  real(8) :: dt, Tprop

! laser fields
  real(8),allocatable :: Act(:,:) 
  real(8) :: E0, omega0, Tpulse0

! MPI
  include 'mpif.h'
  integer :: Myrank,Nprocs,ierr
  real(8) :: Time_start,Time_now
  integer :: nk_s, nk_e, nk_average, nk_remainder

end module global_variables
!---------------------------------------------------------------
program main
  use global_variables
  implicit none

  call MPI_init(ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)
  Time_start=MPI_WTIME()

  call input_variables
  call initialize

  call calc_ground_state
  call time_propagation


  call MPI_finalize(ierr)
  
end program main
!---------------------------------------------------------------
subroutine input_variables
  use global_variables
  implicit none

! physical parameters
  eps_b = 0.5d0*5.9d0/27.2114d0
  eps_n = -0.5d0*5.9d0/27.2114d0
  t0_hop = -2.64d0/27.2114d0

! number of grid points
  nk1 = 32
  nk2 = 32

! lattice constant
  a_lattice = 2.5d0/0.5291772d0 !! 2.5 AA

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
  Tprop = 10d0/0.024189d0
  dt = 0.02d0
  nt = aint(Tprop/dt) + 1

  E0       = 0.012d0
  omega0   = 0.95d0/27.2114d0  !ev
  Tpulse0  = 10d0/0.024189d0   !fs


end subroutine input_variables
!---------------------------------------------------------------
subroutine initialize
  use global_variables
  implicit none
  real(8) :: volume
  integer :: ik1, ik2, ik

  nk = nk1*nk2


  nk_average = nk/nprocs
  nk_remainder = mod(nk,nprocs)
  if(nk < nprocs)call error_finalize('Error: nk < # of MPI processes')
  if(myrank+1 <= nk_remainder)then
    nk_s = 1 + myrank*(nk_average+1)
    nk_e = nk_s + (nk_average + 1) -1
  else
    nk_s = 1 + nk_remainder*(nk_average+1) + nk_average*(myrank-nk_remainder)
    nk_e = nk_s + nk_average -1
  end if

!  write(*,*)myrank,nk_s,nk_e

  allocate(kx0(nk),ky0(nk),kxt(nk),kyt(nk))
  allocate(zpsi(2,nk_s:nk_e),eps_bk(2,nk))

  volume = (a_vec(1,1)*a_vec(2,2)-a_vec(2,1)*a_vec(1,2))*1d0
  
  b_vec(1,1) = 2d0*pi*(a_vec(2,2)*1d0)/volume
  b_vec(2,1) = 2d0*pi*(-a_vec(1,2)*1d0)/volume

  b_vec(1,2) = 2d0*pi*(-1d0*a_vec(2,1))/volume
  b_vec(2,2) = 2d0*pi*(1d0*a_vec(1,1))/volume

  if(myrank == 0)then
    write(*,*)sum(a_vec(:,1)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,1)*b_vec(:,2))/(2d0*pi)
    write(*,*)sum(a_vec(:,2)*b_vec(:,1))/(2d0*pi),sum(a_vec(:,2)*b_vec(:,2))/(2d0*pi)
  end if

  ik = 0
  do ik1 = 0, nk1-1
    do ik2 = 0, nk2-1
      ik = ik + 1

      kx0(ik) = b_vec(1,1)*ik1/dble(nk1) + b_vec(1,2)*ik2/dble(nk2)
      ky0(ik) = b_vec(2,1)*ik1/dble(nk1) + b_vec(2,2)*ik2/dble(nk2)

    end do
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

  eps_bk_l = 0d0
  do ik = nk_s, nk_e
    
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
    eps_bk_l(:,ik) = eps_t(:)
    
  end do

  call MPI_Allreduce(eps_bk_l, eps_bk, 2*nk, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

  if(myrank == 0)then
    open(20,file='band_map.out')
    ik = 0
    do ik1 = 0, nk1-1
      do ik2 = 0, nk2-1
        ik=ik+1
        write(20,"(999e26.16e3)")kx0(ik),ky0(ik),eps_bk(1:2,ik)
      end do
      write(20,*)
    end do
    close(20)
  end if

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
subroutine time_propagation
  use global_variables
  implicit  none
  real(8) :: Etot, jt(2)
  real(8) :: Act_t(2)
  real(8),allocatable :: pop_dist(:,:) ! pop_dist(nk,2)
  integer :: it
  integer :: ik, ik1, ik2

  allocate(pop_dist(nk,2))

  call init_laser
  
  if(myrank == 0) open(20,file='total_energy.out')
  if(myrank == 0) open(21,file='current.out')

  do it = 0, nt
! |zpsi(t+dt)> = exp(-zi*dt*H(t+dt/2)) |zpsi(t)>
    Act_t(:) = 0.5d0*(Act(:,it+1)+Act(:,it))  ! Act(t+dt/2)=0.5*(Act(t+dt)+Act(t))

    kxt(:) = kx0(:) + Act_t(1)
    kyt(:) = ky0(:) + Act_t(2)

    call dt_evolve 

    Act_t(:) = Act(:,it+1)
    kxt(:) = kx0(:) + Act_t(1)
    kyt(:) = ky0(:) + Act_t(2)
    call calc_current(jt)
    call calc_energy(Etot)
    if(myrank == 0) write(20,"(2e26.16e3)")dt*it,Etot
    if(myrank == 0) write(21,"(999e26.16e3)")dt*it,jt(:)

  end do
  if(myrank == 0) close(20)
  if(myrank == 0) close(21)

  call calc_pop_dist(pop_dist)
  if(myrank == 0)then
    open(22,file='pop_dist_final.out')
    ik = 0
    do ik1 = 0, nk1-1
      do ik2 = 0, nk2-1
        ik=ik+1
        write(22,"(999e26.16e3)")kxt(ik),kyt(ik),pop_dist(ik,1:2)
      end do
      write(22,*)
    end do
    close(22)
  end if


end subroutine time_propagation
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

  if(myrank == 0)then
    open(20,file='laser.out')
    do it = 0, nt
      tt = dt*it
      write(20,"(999e26.16)")tt,Act(:,it)
    end do
    close(20)
  end if

end subroutine init_laser
!---------------------------------------------------------------  
subroutine dt_evolve
  use global_variables
  implicit none
  integer :: ik
  complex(8) :: zham(2,2), zfk, zvec(2), zhvec(2)
  real(8) :: kx_t, ky_t
  complex(8) :: zfactor
  integer :: iexp
  


  do ik = nk_s, nk_e

    kx_t = kxt(ik)
    ky_t = kyt(ik)

    zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

    zham(1,1) = eps_b
    zham(1,2) = t0_hop*zfk
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_n

! |zpsi(t)>    
! |zpsi(t+dt)> = exp(-zi*dt*H(t+dt/2))  |zpsi(t)>    
!== Taylor expansion ==
! exp(-zi*dt*H(t+dt/2)) = 1 -zi*dt*H -0.5d0....

    zvec(:) = zpsi(:,ik)
    zfactor = 1d0
    do iexp = 1, 4
      zfactor = zfactor*(-zi*dt)/iexp
      zhvec = matmul(zham, zvec) 
      zpsi(:,ik) = zpsi(:,ik) + zfactor*zhvec(:)

      zvec = zhvec  
    end do
!== Taylor expansion ==

! |zpsi(t+dt)>
  end do


end subroutine dt_evolve
!---------------------------------------------------------------  
subroutine calc_energy(Etot)
  use global_variables
  implicit none
  real(8),intent(out) :: Etot
  real(8) :: Etot_l
  complex(8) :: zham(2,2), zfk, zvec(2), zhvec(2)
  real(8) :: kx_t, ky_t
  integer :: ik

  Etot_l = 0d0
  do ik = nk_s, nk_e

    kx_t = kxt(ik)
    ky_t = kyt(ik)

    zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

    zham(1,1) = eps_b
    zham(1,2) = t0_hop*zfk
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_n

    zvec(:) = zpsi(:,ik)
    zhvec(:) = matmul(zham,zvec)

!Etot = Etot + <zpsi|H|zpsi>
    Etot_l = Etot_l + sum(conjg(zpsi(1:2,ik))*zhvec(1:2))
!    Etot = Etot + sum(abs(zpsi(1:2,ik))**2)

  end do
  Etot_l = Etot_l/nk
  call MPI_Allreduce(Etot_l, Etot, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)  

end subroutine calc_energy
!---------------------------------------------------------------  
subroutine calc_current(jt)
  use global_variables
  implicit  none
  real(8),intent(out) :: jt(2)
  real(8) :: jt_l(2)
  complex(8) :: zJop_x(2,2),zJop_y(2,2), zfk, zvec(2)
  complex(8) :: zJvec_x(2),zJvec_y(2)
  real(8) :: kx_t, ky_t
  integer :: ik

  jt_l = 0d0
  do ik = nk_s, nk_e

    kx_t = kxt(ik)
    ky_t = kyt(ik)

    zJop_x(1,1) = 0d0
    zJop_x(2,2) = 0d0
    zJop_x(1,2) = -t0_hop*( &
      zi*delta_vec(1,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
     +zi*delta_vec(1,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
     +zi*delta_vec(1,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) )

    zJop_x(2,1) = conjg(zJop_x(1,2))

    zJop_y(1,1) = 0d0
    zJop_y(2,2) = 0d0
    zJop_y(1,2) = -t0_hop*( &
      zi*delta_vec(2,1)*exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
     +zi*delta_vec(2,2)*exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
     +zi*delta_vec(2,3)*exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) )

    zJop_y(2,1) = conjg(zJop_y(1,2))


    zvec(:) = zpsi(:,ik)

    zJvec_x = matmul(zJop_x,zvec)
    zJvec_y = matmul(zJop_y,zvec)

    jt_l(1) = jt_l(1) + sum(conjg(zpsi(:,ik))*zJvec_x(:))
    jt_l(2) = jt_l(2) + sum(conjg(zpsi(:,ik))*zJvec_y(:))

  end do

  jt_l = jt_l/nk
  call MPI_Allreduce(jt_l, jt, 2, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine calc_current
!---------------------------------------------------------------  
subroutine calc_pop_dist(pop_dist)
  use global_variables
  implicit none
  real(8),intent(out):: pop_dist(nk,2)
  real(8) :: pop_dist_l(nk,2)
  complex(8) :: zham(2,2), zfk, zvec(2), zphi(2,2)
  real(8) :: eps_t(2)
  real(8) :: kx_t, ky_t
  integer :: ik

  pop_dist_l = 0d0
  do ik = nk_s, nk_e

    kx_t = kxt(ik)
    ky_t = kyt(ik)

    zfk = exp(zi*(kx_t*delta_vec(1,1)+ky_t*delta_vec(2,1))) &
        + exp(zi*(kx_t*delta_vec(1,2)+ky_t*delta_vec(2,2))) &
        + exp(zi*(kx_t*delta_vec(1,3)+ky_t*delta_vec(2,3))) 

    zham(1,1) = eps_b
    zham(1,2) = t0_hop*zfk
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_n  

    call calc_eig_vec_2x2(zham, zphi, eps_t)

    pop_dist_l(ik,1) = abs(sum(conjg(zphi(:,1))*zpsi(:,ik)) )**2
    pop_dist_l(ik,2) = abs(sum(conjg(zphi(:,2))*zpsi(:,ik)) )**2

  end do

  call MPI_Allreduce(pop_dist_l, pop_dist, 2*nk, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)

end subroutine calc_pop_dist
!-------------------------------------------------------------
subroutine error_finalize(message_t)
  use global_variables
  implicit none
  character(*),intent(in) :: message_t

  if(myrank==0)write(*,"(A)")trim(message_t)
  call MPI_Finalize(ierr)
  stop

end subroutine error_finalize
