program ising
integer :: i, j, m, n ! dummy integers
integer, allocatable :: A(:,:) ! matrix containing spins
integer :: nrows, ncols ! number of rows and cols of A
real :: temp, beta ! temperature, inverse temperature
integer :: ConfigType ! starting configuration type
integer :: npass ! number of passes for MC algorithm
integer :: ipass ! the current pass number
integer :: nequil ! number of equilibration steps
integer :: trial_spin ! values of changed spin
real :: high_temp ! starting temp for scan
real :: low_temp ! final temp for scan
real :: temp_interval ! interval between scan points
integer :: nscans ! number of scans (each at diff T)
integer :: iscan ! current scan number
real :: deltaU ! change in energy between 2 configs
real :: log_eta ! log of random number to compare to
real :: magnetization ! magnetization of all spins in lattice
real :: magnetization_ave ! cumulative average magnetization
real :: magnetization2_ave ! cumulative average of mag. squared
real :: energy ! energy of all spins in lattice
real :: energy_ave ! cumulative average of energy
real :: energy2_ave ! cumulative average of energy squared
integer :: output_count ! # times things have been added to averages

nscans=10000

print*, "________________MONTE CARLO 2D ISING MODEL________________"
print*, "Monte Carlo Statistics for 2D Ising Model with"
print*, " periodic boundary conditions."

open(unit=11,file='ising.in',status='old',action='read',form='unformatted')
read(11,*);read(11,*) nrows
read(11,*);read(11,*) ncols
read(11,*);read(11,*) npass
read(11,*);read(11,*) nequil
read(11,*);read(11,*) high_temp
read(11,*);read(11,*) low_temp
read(11,*);read(11,*) temp_interval
read(11,*);read(11,*) ConfigType
close(11)

allocate(A(nrows+2,ncols+2))

open(unit=33,file='magnetization',status='replace',action='write')
write(33,*) "temp ave_magnetization ave_magnetization^2 susceptibility"

open(unit=34,file='energy',status='replace',action='write')
write(34,*) "temp ave_energy ave_energy^2 C_v"

scan_loop: do iscan = 1, nscans
temp = high_temp - temp_interval*(iscan-1)
print*, "Running program for T =", temp

! Initialize variables
beta = 1.0/temp
output_count = 0
energy_ave = 0.0
energy2_ave = 0.0
magnetization_ave = 0.0
magnetization2_ave = 0.0

! Set up the initial spin configuration.
select case(ConfigType)

case(1) ! checkerboard setup
A(1,1) = 1
do i = 1, nrows+1
A(i+1,1) = -A(i,1)
enddo
do j = 1, ncols+1
A(:,j+1) = -A(:,j)
enddo

case(2) ! interface
do i = 1, nrows+2
do j = 1, (ncols+2)/2
A(i,j) = 1
enddo
do j = (ncols+2)/2 + 1, ncols+2
A(i,j) = -1
enddo
enddo
 
case(3) ! unequal interface
do i = 1, nrows+2
do j = 1, (ncols+2)/4
A(i,j) = 1
enddo
do j = (ncols+2)/4 + 1, ncols+2
A(i,j) = -1
enddo
enddo
case default
print*, "Error! Check ConfigType parameter in ising.in"
stop
end select

! Main loop containing Monte Carlo algorithm:
MC_passes: do ipass = 0, npass
! If ipass is greater than nequil (the number of equilibration steps),
! calculate the magnetization and energy:
if (ipass > nequil) then
output_count = output_count + 1
magnetization = sum(A(2:nrows+1,2:nrows+1))/(ncols*nrows*1.0)
magnetization_ave = magnetization_ave + magnetization
magnetization2_ave = magnetization2_ave + magnetization**2

energy = 0.0
do i = 2, nrows + 1
do j = 2, ncols + 1
energy = energy - A(m,n)*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))
enddo
enddo
energy = energy/(ncols*nrows*2.0)
energy_ave = energy_ave + energy
energy2_ave = energy2_ave + energy**2
endif

m = nint((nrows-1)*ranf() + 2) ! choose a random row
n = nint((ncols-1)*ranf() + 2) ! choose a random column
trial_spin = -A(m,n) ! trial spin value

deltaU = -trial_spin*(A(m-1,n)+A(m+1,n)+A(m,n-1)+A(m,n+1))*2
log_eta = dlog(ranf() + 1.0d-10) ! random number 0-1 (+ tiny offset)
if (-beta*deltaU > log_eta) then
A(m,n) = trial_spin
if (m == 2) A(nrows+2,n) = trial_spin
if (m == nrows+1) A(1,n) = trial_spin
if (n == 2) A(m,ncols+2) = trial_spin
if (n == ncols+1) A(m,1) = trial_spin
endif
enddo MC_passes

! Write final spin array to output file
write(33,*) temp, abs(magnetization_ave/output_count), &
magnetization2_ave/output_count, &
beta*(magnetization2_ave/output_count - (magnetization_ave/output_count)**2)
write(34,*) temp, energy_ave/output_count, energy2_ave/output_count, &
(beta**2)*(energy2_ave/output_count - (energy_ave/output_count)**2)

enddo scan_loop
close(33)
close(34)

print*, "Program ising.f90 complete!"

end program ising
