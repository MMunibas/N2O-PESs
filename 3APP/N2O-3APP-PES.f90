module rep_ker
implicit none
real*8, parameter :: dk26f1 = 1.0d0/14.0d0, dk26f2 = 1.0d0/18.0d0, &
dk24f1 = 2.0d0/15.0d0, dk24f2 = 2.0d0/21.0d0, dk25f1 = 2.0d0/21.0d0,&
dk25f2 = 1.0d0/14.0d0, akf1=2.0d0/3.0d0

contains

function drker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker24, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker24 = dk24f1/xl**5 - dk24f2*xs/xl**6

end function drker24

function ddrker24(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker24, xl, xs

if (x .lt. xi) then
  ddrker24 = -dk24f2/xi**6
else
  ddrker24 = -5.0d0*dk24f1/x**6 + 6.0d0*dk24f2*xi/x**7
end if

end function ddrker24

function drker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker25, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker25 = dk25f1/xl**6 - dk25f2*xs/xl**7

end function drker25

function ddrker25(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker25, xl, xs

if (x .lt. xi) then
  ddrker25 = -dk25f2/xi**7
else
  ddrker25 = -6.0d0*dk25f1/x**7 + 7.0d0*dk25f2*xi/x**8
end if

end function ddrker25

function drker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: drker26, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

drker26 = dk26f1/xl**7 - dk26f2*xs/xl**8

end function drker26

function ddrker26(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: ddrker26, xl, xs

if (x .lt. xi) then
  ddrker26 = -dk26f2/xi**8
else
  ddrker26 = -7.0d0*dk26f1/x**8 + 8.0d0*dk26f2*xi/x**9
end if

end function ddrker26

function atker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: atker23, xl, xs

xl = x
xs = xi
if (x .lt. xi) then
  xl = xi
  xs = x
end if

atker23 = 1.0d0 + xs*xl + 2.0d0*xs**2*xl  - akf1*xs**3

end function atker23

function datker23(x,xi)
implicit none
real*8, intent(in) :: x, xi
real*8 :: datker23, xl, xs

if (x .lt. xi) then
  datker23 = xi + 4.0d0*x*xi  - 3.0d0*akf1*x**2
else
  datker23 = xi + 2.0d0*xi**2
end if

end function datker23


end module rep_ker

module surface2
implicit none
real*8, allocatable, dimension(:,:) :: asy_array1, asy_array2, & 
darray1, darray2
integer :: na1, na2, nda1, nda2

contains

subroutine pes3d(r12,r23,r31,totener,dvdr)
implicit none
real*8,parameter :: m1 = 14.003074d0, m2 = m1, m3 = 15.99491462d0
real*8, intent(in) :: r12, r23, r31
real*8, dimension(:), intent(out) :: dvdr(3)
real*8, intent(out) :: totener
real*8 :: capr, theta
real*8, dimension(:,:) :: derv(3,3), deri(3,3), dw(3,3)
real*8, dimension(:) :: derj(3), w(3), ener(3), r(3), m(3), bigr(3)

!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! N---------N
! 1   r1    2

!=====================================================
!
!PES1                     PES2            PES3
!        3             3                 3
!        O             O                 O
!  r2   /               \               /
!      /  r3             \r1        r1 /
!theta/             theta/\           /\ theta r2
!N___/_____N         r2 /  \         /  \
!1   r1    2           /    N      1N    \
!                     / r3  2             \
!                    N1                r3  N2

    r(1)=r12
    r(2)=r31
    r(3)=r23
    m(1) = m2
    m(2) = m1
    m(3) = m3
    call crdtrf(r,m,capr,theta,derv)
    call calcener(capr,r12,theta, ener(1), derj, 1)
    bigr(1) = capr
        
    deri(1,1) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
    deri(1,3) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
    deri(1,2) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)
  
    r(1)=r23
    r(2)=r31
    r(3)=r12
    m(1) = m2
    m(2) = m3
    m(3) = m1
    call crdtrf(r,m,capr,theta,derv)
    call calcener(capr,r23,theta, ener(2), derj, 2)
    bigr(2) = capr

    deri(2,2) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
    deri(2,3) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
    deri(2,1) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

    r(1)=r31
    r(2)=r23
    r(3)=r12
    m(1) = m1
    m(2) = m3
    m(3) = m2
    call crdtrf(r,m,capr,theta,derv)
    call calcener(capr,r31,theta, ener(3), derj, 3)
    bigr(3) = capr

    deri(3,3) = derj(1)*derv(1,1) + derj(2)*derv(2,1) + derj(3)*derv(3,1)
    deri(3,2) = derj(1)*derv(1,2) + derj(2)*derv(2,2) + derj(3)*derv(3,2)
    deri(3,1) = derj(1)*derv(1,3) + derj(2)*derv(2,3) + derj(3)*derv(3,3)

    call weights(r12, r23, r31, w, dw, bigr)
    totener = sum(ener*w)
   
    dvdr(1) = sum(deri(:,1)*w)+sum(dw(:,1)*ener)
    dvdr(2) = sum(deri(:,2)*w)+sum(dw(:,2)*ener)
    dvdr(3) = sum(deri(:,3)*w)+sum(dw(:,3)*ener)

end subroutine pes3d

subroutine crdtrf(r,m,capr,theta, derv)
implicit none
real*8, dimension(:), intent(in) :: r(3), m(3)
real*8, dimension(:,:), intent(out) :: derv(3,3)
real*8, intent(out) :: capr, theta
real*8, parameter :: pc = sqrt(epsilon(1.0d0))
real*8 :: cmu1, cmu2, rcm1, rcm2

!        3
!       /|\
!      / | \
!   r3/ R|  \r2
!    /   |th \
!   /____|____\ 
!  1    r1    2

cmu1 = m(2) / (m(1)+m(2))
cmu2 = m(1) / (m(1)+m(2))

rcm1 = r(1) * cmu1
rcm2 = r(1) * cmu2

capr = sqrt (r(2)**2/r(1)*rcm1 + r(3)**2/r(1)*rcm2 - rcm1*rcm2 )
if (abs(capr) < pc) capr = pc

theta = (rcm2**2+capr**2-r(2)**2)/2.0d0/rcm2/capr
theta=min(1.0d0,max(-1.0d0,theta))
theta = acos(theta)

derv(1,1) = -cmu1*cmu2*r(1)/capr

derv(1,2) = r(2)*cmu1/capr

derv(1,3) = r(3)*cmu2/capr

derv(2,1) = 1.0d0
derv(2,2) = 0.0d0
derv(2,3) = 0.0d0

derv(3,1) = (derv(1,1)/capr*cos(theta)+cos(theta)/r(1)-(capr*derv(1,1)+rcm2*cmu2)/rcm2/capr)&
            /sqrt(1.0d0-cos(theta)**2)

derv(3,2) = (r(2)/rcm2/capr-derv(1,2)/rcm2+cos(theta)/capr*derv(1,2))/sqrt(1.0d0-cos(theta)**2)

derv(3,3) = (cos(theta)/capr*derv(1,3)-derv(1,3)/rcm2)/sqrt(1.0d0-cos(theta)**2)

return

end subroutine crdtrf

subroutine weights(d1,d2,d3,w,dw, bigr)
implicit none
real*8, intent(in) :: d1, d2, d3
real*8, dimension(:), intent(out) :: w(3)
real*8, dimension(:), intent(in) :: bigr(3)
real*8, dimension(:,:), intent(out) :: dw(3,3)
integer, parameter :: power = 2
real*8, parameter :: wtol = epsilon(1.0d0), dr1 = 1.2d0, dr2 = 1.14d0, dr3 = 1.14d0
real*8 :: wsum, pw1, pw2, pw3, r1, r2, r3

r1=d1
r2=d2
r3=d3

wsum = 0.0d0
do while(wsum < wtol)
  pw1 = exp(-(r1/dr1)**power)
  pw2 = exp(-(r2/dr2)**power)
  pw3 = exp(-(r3/dr3)**power)
  wsum = pw1+pw2+pw3
  if(wsum < wtol) then
    r1=r1-0.1d0
    r2=r2-0.1d0
    r3=r3-0.1d0
  end if
end do

  w(1) = pw1 / wsum
  w(2) = pw2 / wsum
  w(3) = pw3 / wsum

  dw(1,1) =-power*((r1/dr1)**(power-1))*(pw1*pw2+pw1*pw3)/dr1/wsum**2
  dw(1,2) = power*((r2/dr2)**(power-1))*pw1*pw2/dr2/wsum**2
  dw(1,3) = power*((r3/dr3)**(power-1))*pw1*pw3/dr3/wsum**2

  dw(2,1) = power*((r1/dr1)**(power-1))*pw2*pw1/dr1/wsum**2
  dw(2,2) =-power*((r2/dr2)**(power-1))*(pw2*pw1+pw2*pw3)/dr2/wsum**2
  dw(2,3) = power*((r3/dr3)**(power-1))*pw2*pw3/dr3/wsum**2


  dw(3,1) = power*((r1/dr1)**(power-1))*pw3*pw1/dr1/wsum**2
  dw(3,2) = power*((r2/dr2)**(power-1))*pw3*pw2/dr2/wsum**2
  dw(3,3) =-power*((r3/dr3)**(power-1))*(pw3*pw1+pw3*pw2)/dr3/wsum**2

return

end subroutine weights

subroutine calcener(capr,smlr,theta, ener, der, sno)
use rep_ker
use RKHS            ! This module needs to be used by your code
implicit none
real*8 :: lambda
real*8, intent(out) :: ener
real*8, intent(in) :: capr, smlr, theta
integer, intent(in) :: sno
real*8, dimension(:), intent(out) :: der(3)
real*8,parameter :: pi = acos(-1.0d0), piby180 = pi/180.0d0
real*8 :: asener, anener, asder, z1, z2
real*8, dimension(:) :: ander(3), x(3)
integer :: kk, ii
type(kernel), save  :: pes21, pes22!, pes23           ! The kernel type is needed to set up and evaluate a RKHS model
logical, save :: stored = .false., kread = .false.
logical, save :: ker1 = .false., ker2 = .false.!, ker3 = .false.
character (len=80), save :: datapath="./"

if (.not. ker1) then
  inquire(file=trim(datapath)//"pes21.kernel", exist=ker1)   ! file_exists will be true if the file exists and false otherwise
end if

if (.not. ker2) then
  inquire(file=trim(datapath)//"pes22.kernel", exist=ker2)   ! file_exists will be true if the file exists and false otherwise
end if

lambda=0.1d-19

if (.not. stored ) then

open(unit=1001,file=trim(datapath)//"asymp2.dat", status = "old")

read(1001,*)na1
allocate(asy_array1(na1,2))
do ii = 1, na1
  read(1001,*)asy_array1(ii,1), asy_array1(ii,2)
end do

read(1001,*)na2
allocate(asy_array2(na2,2))
do ii = 1, na2
  read(1001,*)asy_array2(ii,1), asy_array2(ii,2)
end do

read(1001,*)nda1
allocate(darray1(nda1,2))
do ii = 1, nda1
  read(1001,*)darray1(ii,1), darray1(ii,2)
end do

read(1001,*)nda2
allocate(darray2(nda2,2))
do ii = 1, nda2
  read(1001,*)darray2(ii,1), darray2(ii,2)
end do

stored = .true.

end if

if (.not. kread) then
!if (ker1 .and. ker2 .and. ker3 ) then
if (ker1 .and. ker2 ) then
  call pes21%load_from_file(trim(datapath)//"pes21.kernel")
  call pes22%load_from_file(trim(datapath)//"pes22.kernel")
!  call pes23%load_from_file(trim(datapath)//"pes23.kernel")
  kread = .true.
else
  call pes21%read_grid(trim(datapath)//"pes21.csv")
!  print*,"IAMHERE"
  call pes21%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
  call pes21%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
  call pes21%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL) 
!
  call pes21%calculate_coefficients_slow(lambda)
!  call pes21%calculate_coefficients_fast()
!
  call pes21%calculate_sums()
!
  call pes21%save_to_file(trim(datapath)//"pes21.kernel")
!
  call pes22%read_grid(trim(datapath)//"pes22.csv")
!  print*,"IAMHERE"
  call pes22%k1d(1)%init(TAYLOR_SPLINE_N2_KERNEL)         ! choose one-dimensional kernel for dimension 1
  call pes22%k1d(2)%init(RECIPROCAL_POWER_N2_M6_KERNEL)   ! choose one-dimensional kernel for dimension 2
  call pes22%k1d(3)%init(RECIPROCAL_POWER_N2_M6_KERNEL)
!
  call pes22%calculate_coefficients_slow(lambda)
!  call pes22%calculate_coefficients_fast()
!
  call pes22%calculate_sums()
!
  call pes22%save_to_file(trim(datapath)//"pes22.kernel")

  kread = .true.
end if
end if

x(1)=(1.0d0-cos(theta))/2.0d0
x(2)=capr
x(3)=smlr

if (sno==1) then
  asener = 0.0d0
  asder=0.0d0
  do kk = 1,na1
    asener = asener + drker26(smlr,asy_array1(kk,1))*asy_array1(kk,2)
    asder = asder + ddrker26(smlr,asy_array1(kk,1))*asy_array1(kk,2)
  end do

  anener = 0.0d0
  ander=0.0d0
  call pes21%evaluate_fast(x,anener,ander)

    ener = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/2.0d0
end if

if (sno==2) then
  asener = 0.0d0
  asder=0.0d0
  do kk = 1,na2
    asener = asener + drker26(smlr,asy_array2(kk,1))*asy_array2(kk,2)
    asder = asder + ddrker26(smlr,asy_array2(kk,1))*asy_array2(kk,2)
  end do

  anener = 0.0d0
  ander=0.0d0
  call pes22%evaluate_fast(x,anener,ander)

    ener = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/2.0d0
end if

if (sno==3) then
  asener = 0.0d0
  asder=0.0d0
  do kk = 1,na2
    asener = asener + drker26(smlr,asy_array2(kk,1))*asy_array2(kk,2)
    asder = asder + ddrker26(smlr,asy_array2(kk,1))*asy_array2(kk,2)
  end do

  anener = 0.0d0
  ander=0.0d0
  call pes22%evaluate_fast(x,anener,ander)

    ener = asener+anener
    der(1) = ander(2)
    der(2) = asder+ander(3)
    der(3) = ander(1)*sin(theta)/2.0d0
end if

return

end subroutine calcener

end module

module callpot2

contains

subroutine n2o3apppes(tmpr, totener, dvdr)
use surface2
implicit none
real*8, dimension (:), intent (in) :: tmpr(3)
real*8, dimension (:), intent (out) :: dvdr(3)
real*8, intent (out) :: totener
real*8, dimension (:) :: xp(3), tmpdvdr(3), r(3)
real*8, parameter :: dx = 0.0005d0
real*8 :: ener
integer :: ii

!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! N---------N
! 1   r1    2
!N+ON 
!NO ON NN
r(1) = tmpr(3)
r(2) = tmpr(1)
r(3) = tmpr(2)

call pes3d(r(1),r(2),r(3),totener,dvdr)
if ((dvdr(1)-1) .eq. dvdr(1) .or. (dvdr(3)-1) .eq. dvdr(3) .or. (dvdr(2)-1) .eq. dvdr(2) ) go to 100
return
!100 print*,"NaN or Infinity"
100 if ((r(1)+r(2)-r(3))<1d-10) then
!==========================================================
!   Forward difference to calculate first derivative      =
!==========================================================
    xp=r
    xp(1)=r(1)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(1)=(ener-totener)/dx

    xp=r
    xp(2)=r(2)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(2)=(ener-totener)/dx

    xp=r
    xp(3)=r(3)-dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(3)=(totener-ener)/dx

    return
end if

if ((r(2)+r(3)-r(1))<1d-10) then
!==========================================================
!   Forward difference to calculate first derivative      =
!==========================================================
    xp=r
    xp(1)=r(1)-dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(1)=(totener-ener)/dx

    xp=r
    xp(2)=r(2)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(2)=(ener-totener)/dx

    xp=r
    xp(3)=r(3)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(3)=(ener-totener)/dx

    return
end if

if ((r(1)+r(3)-r(2))<1d-10) then
!==========================================================
!   Forward difference to calculate first derivative      =
!==========================================================
    xp=r
    xp(1)=r(1)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(1)=(ener-totener)/dx

    xp=r
    xp(2)=r(2)-dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(2)=(totener-ener)/dx

    xp=r
    xp(3)=r(3)+dx
    call pes3d(xp(1),xp(2),xp(3),ener,tmpdvdr)
    dvdr(3)=(ener-totener)/dx

    return
end if

end subroutine n2o3apppes

subroutine diatno(r,ener,der)
use surface2
use rep_ker
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: ener, der
integer :: kk

ener=0.0d0
der=0.0d0
do kk = 1, nda2
    ener = ener + drker26(r,darray2(kk,1))*darray2(kk,2)
    der = der + ddrker26(r,darray2(kk,1))*darray2(kk,2)
end do

return

end subroutine diatno

subroutine diatn2(r,ener,der)
use surface2
use rep_ker
implicit none
real*8, intent(in) :: r
real*8, intent(out) :: ener, der
integer :: kk

ener=0.0d0
der=0.0d0
do kk = 1, nda1
    ener = ener + drker26(r,darray1(kk,1))*darray1(kk,2)
    der = der + ddrker26(r,darray1(kk,1))*darray1(kk,2)
end do

return

end subroutine diatn2

end module
