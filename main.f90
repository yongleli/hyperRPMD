program main
use system
implicit none
real(kind=8) :: mass
real(kind=8) :: Temperature
real(kind=8) :: beta
real(kind=8) :: TT
real(kind=8) :: dt
real(kind=8),parameter :: fact = 0.00000316678d0
real(kind=8),parameter :: constants_kB = 1.3806504d-23
integer :: Nbeads
integer :: N_andersen
integer :: Nstep
real(kind=8),allocatable :: x_init(:), p_init(:)
real(kind=8) :: omega, TWO_ZPE, dt_orig
character(len=100) :: filename, filename2
character(len=5)   :: x1, S_Nbeads
character(len=12)  :: S_Temp,x2
character(len=20)  :: fmt
character(len=3)  :: ensemble
integer :: i,k, i_f1, i_f2

!!Nbeads = 
!!Nstep =  40
Temperature = 300.0d0
call getarg(1,S_Nbeads)
read(S_Nbeads,'(I5)') Nbeads
!!call getarg(2,S_Temp)
!!read(S_Temp, '(F8.3)') Temperature
if (Nbeads <= 1) Nbeads = 1
if (Temperature < 0.0d0 ) Temperature =0.0d0

Temperature = 1.0
beta = 1.0
TT = 1.0
dt = 0.01 !! dt is 0.25 fs.
mass = 1.0d0

print *, "Nbeads: ", Nbeads
print *, "Temperature: ", Temperature
!!Nstep = 20000
!!Nstep = 20000
Nstep = 20000

fmt = '(I5.5)' ! an integer of width 5 with zeros at the left
write (x1,fmt) Nbeads ! converting integer to string using a 'internal file'
fmt = '(F8.1)' 
write (x2,fmt) Temperature
filename ='Output'//trim(x1)//'_'//trim(adjustl(x2))//'_centroid.dat'
filename2='Output'//trim(x1)//'_'//trim(adjustl(x2))//'_bead.dat'
allocate(x_init(Nbeads))
allocate(p_init(Nbeads))

!!!!beta = 4.35974417e-18 / (constants_kB * Temperature)
!!!!TT = Temperature * fact 
!!!!dt_orig = 0.1 !! dt is 0.1 fs.
!!!!dt = 41.3414d0*dt_orig !! 1 fs = 41.3414 a.u.
!!!!mass = 1836.0d0*0.5d0
N_andersen = int(dsqrt(dble(Nstep)))

!!!! Initial condition
x_init = 0.0d0
call sample_momentum(p_init, mass, beta, Nbeads,TT)
call sample_momentum(p_init, mass, beta, Nbeads,TT)
call sample_momentum(p_init, mass, beta, Nbeads,TT)

!!!! MD major part
open(100,file=filename)
open(200,file=filename2)
ensemble = 'NVE'
call MD_run(x_init, p_init, Nstep, N_andersen, dt, &
& mass,Nbeads,beta,100,200,TT,ensemble)

close(100)
close(200)

deallocate(x_init)
deallocate(p_init)
end program main
