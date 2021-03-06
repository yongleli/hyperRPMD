module system
implicit none
real(kind=8),parameter :: pi=3.141592653589793d0
contains

!!! Cayley.
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
real(kind=8)               :: wk2, denom, denom_m1, c11, c12, c21, c22, v_new

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

do k = 1, Nbeads
            if ( mod((k-1),2)==0 ) then
              wk = twown * dsin( dble(k-1) * pi_n *0.5d0 )
            else
              wk = twown * dsin( k * pi_n * 0.5d0 )
            end if
            wk2 = wk*wk
            denom = 4.0d0 + dt**2 * wk2
            if(denom<1.0E-12 .or. denom>1.0E+12) then
              print *, "Singular!"
              print *, denom
              stop
            end if
            denom_m1 = 1.0d0/denom

            c11 = 8.0d0*denom_m1 - 1.0d0
            c12 = 4.0d0*dt*denom_m1
            c21 = -4.0d0*dt*denom_m1*wk2
            c22 = c11
            !!print '(4ES18.6)', c11, c12, c21, c22

            q_new    = c11*q(k) + c12*v(k)
            v_new    = c21*q(k) + c22*v(k)
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

end module system
