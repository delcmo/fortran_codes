module mesh
integer, parameter :: n = 10
end module
program array
use mesh
implicit none
! Declare variables: matrix and vector
real, allocatable, dimension(:,:) :: kernel
real, dimension(n,n) :: mass
real, dimension(n,1) :: U,product
integer i_var,j_var, AllocateStatus
! Loop to set mass matrix entries
do i_var = 1,n,1
   do j_var = 1,n,1
      mass(i_var,j_var)=i_var+j_var
      print *, "Value of the matrix: M(", i_var, ",",j_var,")",mass(i_var,j_var)
   end do
end do
! Loop to set vector
do i_var = 1,n,1
   U(i_var,1) = i_var
   print *, "Value of the vector U(", i_var, ")=",U(i_var,1)
end do
! Product of mass and U
product=matmul(mass,U)
! Output product
do i_var = 1,n,1
   print *, "product=", product(i_var,1)
end do
! Allocate memory
allocate(kernel(n,n), stat=AllocateStatus)
if(AllocateStatus /= 0) then
print *, "error"
end if
! Output kernel
do i_var = 1,n,1
   do j_var = 1,n,1
   print *, "kernel=", kernel(i_var,j_var)
   end do
end do
end program array
