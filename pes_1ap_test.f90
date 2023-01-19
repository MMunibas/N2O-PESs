program testpes
use surface
!use callpot2
implicit none
real*8 :: capr, smlr, theta, v
real*8, dimension(:) :: dvdr(3)
real*8, parameter :: pi=acos(-1.0d0)

!!!!!!!!!!!!!!!!!!!!!!!!
!        3       
!        O      
!       /      
!      / capr  
!theta/      
!N___/_____N
!1   smlr  2 
!            
!!!!!!!!!!!!!!!!!!!!!!! 

capr=3.5d0   !in bohr
smlr=2.2d0   !in bohr
theta=179.0d0*pi/180.0d0 !in radian

call calcener(capr,smlr,theta,v,dvdr)

write(*,*)"Energy = ", v, "Hartree"

end
