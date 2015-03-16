program test
implicit none
! Interface
interface
! function
function Area_circle_fnct(r)
real Area_circle_fnct
real, intent(in) :: r
end function Area_circle_fnct
! subroutine
subroutine compute_area(r2, area)
real, intent(in) :: r2
real, intent(out) :: area
end subroutine compute_area
end interface
! Declarations
real answer_real, x, y, radius, area
integer answer_int, a, b
! Compute the sum of two reals
print *, "Enter two reals"
read *, x
read *, y
answer_real=x+y
print *, "The answer is ", answer_real
! Compute the sum of two integers
print *, "Enter two integers"
read *, a
read *, b
answer_int=a+b
print *, "The answer is ", answer_int
! Compute the radius of a circle using a function
print *, "Enter a radius"
read *, radius
print *, "The area of the circle of radius ", radius, " is ",  Area_circle_fnct(radius)
! Compute the area of a circle using a subroutine
print *, "Enter a radius"
read *, radius
call compute_area(radius, area)
print *, "The area of the circle of radius ", radius, " is ", area
end program test
!
!----------------------------------------------------------------
! End of program
!----------------------------------------------------------------
!
! Function to compute the area of a circle
!
function Area_circle_fnct(r)
implicit none
real :: Area_circle_fnct
real, intent(in) :: r
real, parameter :: Pi = 3.1415927
Area_circle_fnct=Pi*r*r
end function Area_circle_fnct
!
! Subroutine to compute the area of a circle
!
subroutine compute_area(r2, area)
implicit none
real, intent(in) :: r2
real, intent(out) :: area
real, parameter :: Pi = 3.1415927
area=Pi*r2*r2
end subroutine compute_area
