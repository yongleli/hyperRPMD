module system
implicit none
real(kind=8),parameter :: pi=3.141592653589793d0
contains
subroutine potential(x,V,dV)
implicit none
real(kind=8),intent(in) :: x
real(kind=8),intent(out) :: V,dV
!!real(kind=8),parameter :: k_force =  0.369127 !! Force constant
real(kind=8),parameter :: k_force =  256.0d0 !! Force constant

V  = 0.5*k_force*(x)**2
dV = k_force*(x)
end subroutine potential

subroutine ana_fft(Y_in,Nbeads)
implicit none
integer,intent(in) :: Nbeads
real(kind=8),intent(inout) :: Y_in(Nbeads)
real(kind=8) :: Y_out(Nbeads)
integer :: k,i, N_half
real(kind=8) :: factor_in

N_half = Nbeads/2
Y_out = 0.0d0

do k=1, Nbeads !! The first term is just the summation of all.
  Y_out(1) = Y_out(1) + Y_in(k)
  do i=2,N_half
    factor_in = 2.0*pi*dble((k-1)*(i-1))/dble(Nbeads)
    Y_out(i) = Y_out(i) + dcos(factor_in)*Y_in(k)
  end do
  Y_out(N_half+1) = Y_out(N_half+1) + (-1)**(k-1)*Y_in(k)
  do i=(N_half+2),Nbeads
    factor_in = 2.0*pi*dble((i-1)*(k-1))/dble(Nbeads)
    Y_out(i) = Y_out(i) + dsin(factor_in)*Y_in(k)
  end do
end do
Y_out = Y_out/dsqrt(dble(Nbeads))

do k = 1, Nbeads
  Y_in(k) = Y_out(k)
end do
end subroutine ana_fft

subroutine ana_ifft(Y_out, Nbeads)
implicit none
integer,intent(in) :: Nbeads
!!real(kind=8),intent(out) :: Y_in(Nbeads)
!!real(kind=8),intent(in):: Y_out(Nbeads)
real(kind=8),intent(inout):: Y_out(Nbeads)
real(kind=8) :: Y_in(Nbeads)
integer :: k,i, N_half
real(kind=8) :: factor_in
real(kind=8) :: Y_temp

N_half = Nbeads/2
Y_in = 0.0d0
!!print *, "N_half+1: ", N_half + 1
do k=1, Nbeads
  Y_in(k) = Y_out(1)
  do i=2,N_half
    factor_in = 2.0*pi*dble((i-1)*(k-1))/dble(Nbeads)
    Y_in(k) = Y_in(k) + 2.0d0*(dcos(factor_in)*Y_out(i)-dsin(factor_in)*Y_out(N_half+i))
  end do

  Y_in(k) = Y_in(k) + (-1.0)**((k-1))*Y_out(N_half+1)
end do

Y_in = Y_in/dsqrt(dble(Nbeads))

do k = 1, Nbeads
 Y_out(k) = Y_in(k)
end do

do k = 3, N_half-1, 2
  Y_temp = Y_out(k)
  Y_out(k) = Y_out(Nbeads-k+2)
  Y_out(Nbeads-k+2) = Y_temp
end do
end subroutine ana_ifft

!!! Cayley, BCOCB. Here only contains COC.
subroutine free_ring_polymer(q,p,Nbeads,mass,dt,beta)
implicit none
real(kind=8),intent(inout) :: q(Nbeads), p(Nbeads)
real(kind=8),intent(in)    :: mass, dt, beta
integer,intent(in)         :: Nbeads
!!real(kind=8)               :: poly(4,Nbeads)
real(kind=8)               :: beta_n, twown, pi_n
!!real(kind=8)               :: q0(Nbeads), p0(Nbeads)
real(kind=8)               :: v(Nbeads)
integer :: N_half, k, j, i, N_half2
real(kind=8)               :: wk, wt, wm, cos_wt, sin_wt, q_new!!, p_new
real(kind=8)               :: wk2, denom, denom_m1, v_new
real(kind=8)               :: c11(Nbeads), c12(Nbeads), c21(Nbeads), c22(Nbeads)
real(kind=8)               :: gamma_in(Nbeads) !! For friction.
real(kind=8)               :: fact_random(Nbeads) !! For friction.
real(kind=8)               :: fact_fric(Nbeads) !! For friction.
real(kind=8)               :: v_xi

        !!! First, transform from p to v^{iso}
        do k = 1, Nbeads
               v(k) = p(k)/(mass) !! v^{iso,(k)}_{i}=p_{i}^{k}/(Pmass_i)
        end do
        ! Transform to normal mode space
        call rfft(q, Nbeads)
        call rfft(v, Nbeads)

!!call ana_fft(q,Nbeads)
!!call ana_fft(p,Nbeads)
!!                call rfft(p(:), Nbeads)
!!                call rfft(q(:), Nbeads)
beta_n = beta / dble(Nbeads)
twown = 2.0d0 / beta_n
pi_n = pi / dble(Nbeads)
!!! Half B.
do k = 1, Nbeads
            if ( mod((k-1),2)==0 ) then
              wk = twown * dsin( dble(k-1) * pi_n *0.5d0 )
            else
              wk = twown * dsin( k * pi_n * 0.5d0 )
            end if

            gamma_in(k) = 2.0d0* wk !! friction parameter.
            fact_fric(k) = dexp(-gamma_in(k)*dt) !! prefactor of O (v(k))
            fact_random(k) = dsqrt((1.0d0-dexp(-2.0d0*gamma_in(k)*dt))/(beta*mass)) !! prefactor of O (\xi)

            wk2 = wk*wk
            denom = 4.0d0 + dt**2 * wk2
            if(denom<1.0E-12 .or. denom>1.0E+12) then
              print *, "Singular!"
              print *, denom
              stop
            end if
            denom_m1 = 1.0d0/denom
            denom_m1 = dsqrt(denom_m1) !! For cay(dt A)^{1/2}!!

            c11(k) =  2.0d0*denom_m1
            c12(k) =  dt*denom_m1
            c21(k) = -dt*denom_m1*wk2
            c22(k) = c11(k)
            !!print '(4ES18.6)', c11, c12, c21, c22

            q_new    = c11(k)*q(k) + c12(k)*v(k)
            v_new    = c21(k)*q(k) + c22(k)*v(k)
            q(k) = q_new
            v(k) = v_new
end do

!!! O. Langevin Thermostat.
do k = 1, Nbeads
   call randomn(v_xi)
   v(k) = fact_fric(k)* v(k) + fact_random(k)*v_xi
end do

!!! Another Half B.
do k = 1, Nbeads
            q_new    = c11(k)*q(k) + c12(k)*v(k)
            v_new    = c21(k)*q(k) + c22(k)*v(k)
            q(k) = q_new
            v(k) = v_new
end do

call irfft(q,Nbeads)
call irfft(v,Nbeads)

        !!! from v^{iso} to p
        do k = 1, Nbeads
               p(k) = v(k)*mass !!
        end do
!!call ana_ifft(q,Nbeads)
!!call ana_ifft(p,Nbeads)
!!call irfft(q(:),Nbeads)
!!call irfft(p(:),Nbeads)
end subroutine free_ring_polymer

subroutine verlet_step(x,p,dt,mass,Nbeads,beta)
implicit none
integer,intent(in)         :: Nbeads
real(kind=8),intent(inout) :: x(Nbeads), p(Nbeads)
real(kind=8),intent(in) :: dt, mass, beta
real(kind=8) :: V, dV
integer :: k

do k=1, Nbeads
  call potential(x(k),V,dV)
  p(k) = p(k) - 0.5d0*dt*dV
end do

if (Nbeads == 1) then
  x = x + p*dt/mass
else
  call free_ring_polymer(x,p,Nbeads,mass,dt,beta)
end if

do k=1, Nbeads
  call potential(x(k),V,dV)
  p(k) = p(k) - 0.5d0*dt*dV
end do
end subroutine verlet_step

subroutine sample_momentum(p, mass, beta, Nbeads,TT)
    implicit none
    integer, intent(in) :: Nbeads
    double precision, intent(in) :: mass
    double precision, intent(in) :: beta
    double precision, intent(in) :: TT
    double precision, intent(out) :: p(Nbeads)

    double precision :: beta_n, dp
    integer :: k

    beta_n = beta / dble(Nbeads)
    dp = sqrt(mass / beta_n)

    do k = 1, Nbeads
        call randomn(p(k))
        p(k) = p(k) * dp
    end do

end subroutine sample_momentum

subroutine MD_run(x_init, p_init, Nstep, N_andersen, dt, mass,Nbeads,beta,i_f1,i_f2,TT,ensemble)
!!!! f1: centroid, f2: bead-bead
implicit none
real(kind=8),intent(in)  :: dt, mass, beta
integer,intent(in)       :: Nstep, N_andersen, Nbeads
integer,intent(in)       :: i_f1, i_f2
real(kind=8),intent(in)  :: x_init(Nbeads), p_init(Nbeads)
character(len=3),intent(in) :: ensemble

real(kind=8)             :: x(Nbeads), p(Nbeads)
real(kind=8)             :: Xc, Pc
real(kind=8)             :: Time
real(kind=8)             :: TT
real(kind=8)             :: ene
integer :: istep, k, N_interval, N_start
!!integer :: i_f1, i_f2
x = x_init
p = p_init
N_interval = 1 !! Output interval.
!!N_start = Nstep/2
N_start = 0
do istep = 1, Nstep

        call verlet_step(x,p,dt,mass,Nbeads,beta)

        Xc = 0.0
        Pc = 0.0
        do k=1, Nbeads
            Xc = Xc + x(k)
            Pc = Pc + p(k)
        end do
        Xc = Xc/dble(Nbeads)
        Pc = Pc/dble(Nbeads)
        Time=(istep*dt)

        !if (istep> N_start .and. mod(istep,N_interval)==0) then
        call calc_ene(x, p, Nbeads, mass, beta, ene)
        write(i_f1,'(ES14.6,3ES18.6)') Time, Xc, Pc, ene/dble(Nbeads)
        do k=1, Nbeads
          write(i_f2,'(ES14.6,2ES18.6)') Time, x(k), p(k)
        end do
        !end if

        !!!!! Andersen thermostat.
        if (ensemble =='NVT') then
          if (mod(istep , N_andersen) == 0) then
              call sample_momentum(p, mass, beta,Nbeads,TT)
          end if
        end if

end do

end subroutine MD_run

subroutine calc_ene(x, p, Nbeads, mass, beta, ene) 
implicit none
integer,intent(in)::Nbeads
real(8),intent(in)::x(Nbeads), p(Nbeads),mass, beta
real(8),intent(out)::ene
integer::ibead
real(8)::V, dV
ene = 0.0d0

do ibead=1,Nbeads-1
  ene = ene + (x(ibead+1) - x(ibead))**2
end do
ene = ene + (x(Nbeads) - x(1))**2
ene = ene * 0.5 * mass * (dble(Nbeads) / beta)**2

do ibead=1,Nbeads
  ene = ene + p(ibead)**2 / mass
  call potential(x(ibead),V,dV)
  ene = ene + V
end do
end subroutine

end module system
