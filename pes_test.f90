program testpes
use callpot1 
!use callpot2
implicit none
real*8 :: rnn, rno, ron, v
real*8, dimension(:) :: dvdr(3), r(3)

!      3
!      O
!     / \
!  r3/   \r2
!   /     \
!  /       \
! N---------N
! 1   r1    2

rnn=2.4d0 !in bohr
rno=4.5d0 !in bohr
ron=2.3d0 !in bohr

r(1)=rnn
r(2)=ron
r(3)=rno

call n2o3appes(r, v, dvdr)
!call n2o3apppes(r, v, dvdr)

write(*,*)"Energy = ", v, "Hartree"
write(*,*)"dV/dr_i (Hartree/bohr) = ", dvdr

end
