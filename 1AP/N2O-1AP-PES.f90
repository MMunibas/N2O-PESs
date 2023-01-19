module surface
!!!!!!N2O 1A' potential energy surface for spectroscopic study

contains
subroutine calcener(capr,smlr,theta, ener, der)
use RKHS            ! This module needs to be used by your code
implicit none
real*8 :: lambda
real*8, intent(out) :: ener
real*8, intent(in) :: capr, smlr, theta
real*8, dimension(:), intent(out) :: der(3)
real*8,parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
real*8 :: asener, anener, asder, z1, z2, angle
real*8, dimension(:) :: ander(3), x(3)
integer :: kk, ii
type(kernel), save  :: pes1           ! The kernel type is needed to set up and evaluate a RKHS model
logical, save :: stored = .false., kread = .false.
logical, save :: ker1 = .false.

!=====================================================
!
!        3       
!        O      
!  r2   /      
!      /  r3  
!theta/      
!N___/_____N
!1   r1    2 
!            
!=====================================================

if (.not. ker1) then
  inquire(file="pes1.kernel", exist=ker1)   ! file_exists will be true if the file exists and false otherwise
end if

lambda=0.1d-16

if (.not. kread) then
  if (ker1) then
    call pes1%load_from_file("pes1.kernel")
    kread = .true.
  else
    call pes1%read_grid("pes1.csv")
!    print*,"IAMHERE"
    call pes1%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
    call pes1%k1d(2)%init(RECIPROCAL_POWER_N2_M5_KERNEL)   ! choose one-dimensional kernel for dimension 2
    call pes1%k1d(3)%init(RECIPROCAL_POWER_N2_M5_KERNEL) 

!    call pes1%calculate_coefficients_slow(lambda)
    call pes1%calculate_coefficients_fast()

    call pes1%calculate_sums()

    call pes1%save_to_file("pes1.kernel")
!
    kread = .true.
  end if
end if

angle=theta
if (theta>pi/2.0d0) angle=pi-theta
x(1)=(1.0d0-cos(angle))/2.0d0
x(2)=capr
x(3)=smlr

  anener = 0.0d0
  ander=0.0d0
  call pes1%evaluate_fast(x,anener,ander)

    ener = anener
    der(1) = ander(2)
    der(2) = ander(3)
    der(3) = ander(1)*sin(angle)/2.0d0

return

end subroutine calcener

end module
