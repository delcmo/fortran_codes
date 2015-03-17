!----------!
!- Module -!
!----------!
module mesh
! final time
double precision, parameter :: end_time = 1.E-2
! number of cells
integer, parameter :: nb_cell = 100
! number of degree of freedom
integer, parameter :: dof = nb_cell+1
! initial left and right values
double precision, parameter :: u_left = 2.
double precision, parameter :: u_right = 1.
! bc values
double precision, parameter :: u_left_bc = 2.
double precision, parameter :: u_right_bc = 1.
! initial position of the step
double precision, parameter :: x_step = 0.5
! length and discretization parameters
double precision, parameter :: L = 1.
double precision, parameter :: delta_x = L / nb_cell
integer, parameter :: p = 1
double precision, parameter :: Jwx = 0.5*delta_x
double precision, parameter :: cfl = 0.5
end module
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
!!-----------------!!
!!- Start program -!!
!!-----------------!!
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
program array
use mesh
implicit none
!!!!!!!!!!!!!!!!!!!!!
!! Start interface !!
!!!!!!!!!!!!!!!!!!!!!
interface
!
!!!! List of functions !!!!
!
! test function
!
function phi(x,i)
double precision phi
double precision, intent(in) :: x
integer, intent(in) :: i
end function phi
!
! grad test function
!
function grad_phi(x,i)
double precision grad_phi
double precision, intent(in) :: x
integer, intent(in) :: i
end function grad_phi
!
! flux function
!
function flux(u,flux_type)
double precision flux
double precision, intent(in) :: u
character(len=6), intent(in) :: flux_type
end function flux
!
! dissipative function
!
function viscous_flux(grad_u,kappa)
double precision viscous_flux
double precision, intent(in) :: grad_u,kappa
end function viscous_flux
!
! Compute viscous coefficients
!
function viscous_coeff(uq,delta_x,dof)
integer, intent(in) :: dof
double precision viscous_coeff
double precision, intent(in) :: uq,delta_x
double precision viscous_flux
end function viscous_coeff
!
!!!! End list of functions !!!!
!
!!!! List of subroutines !!!!
!
! Weight and quadrature points
!
subroutine quad_pts_and_weight(p, xq, wq)
integer, intent(in) :: p
double precision, dimension (p+1,1), intent(out) :: xq, wq
end subroutine
!
! Compute mass matrix
!
subroutine mass_matrix(p,nb_cell,xq,wq,Jwx,M)
integer, intent(in) :: p,nb_cell
double precision, intent(in) :: Jwx
double precision, dimension (p+1,1), intent(in) :: xq,wq
double precision, dimension (p*nb_cell+1,p*nb_cell*1), intent(out) :: M
end subroutine mass_matrix
!
! Steady-state residual
!
subroutine stt_residual(p,nb_cell,xq,wq,Jwx,U,stt_R)
integer, intent(in) :: p,nb_cell
double precision, intent(in) :: Jwx
double precision, dimension(p+1,1), intent(in) :: xq,wq
double precision, dimension(p*nb_cell+1,1), intent(in) :: U
double precision, dimension(p*nb_cell+1,1), intent(out) :: stt_R
end subroutine stt_residual
!
! Apply boundry conditions
!
subroutine apply_bc(dof,u_left_bc,u_right_bc,bc_type,U,stt_R)
integer, intent(in) :: dof
double precision, intent(in) :: u_left_bc,u_right_bc
character(len=7), intent(in) :: bc_type
double precision, dimension(dof,1), intent(in) :: U
double precision, dimension(dof,1), intent(inout) :: stt_R
end subroutine apply_bc
!
! Compute the numerical solution at each time step
!
subroutine compute_next_time_step(U,U_old, invMass, stt_residual, delta_t,dof)
integer, intent(in) :: dof
double precision, dimension(dof,1), intent(in) :: U_old
double precision, dimension(dof,1), intent(out) :: U
double precision, dimension(dof,1), intent(in) :: stt_residual
double precision, dimension(dof,dof), intent(in) :: invMass
double precision, intent(in) :: delta_t
end subroutine compute_next_time_step
!
! Compute time step based on cfl number
!
subroutine delta_t_from_cfl(delta_t,delta_x,U,cfl,dof)
integer, intent(in) :: dof
double precision, dimension(dof,1), intent(in) :: U
double precision, intent(in) :: cfl,delta_x
double precision, intent(out) :: delta_t
end subroutine delta_t_from_cfl
!
!!!! End list of subroutines !!!!
!
end interface
!
!!!!!!!!!!!!!!!!!!!
!! End interface !!
!!!!!!!!!!!!!!!!!!!
!
! Declare variables: matrix and vectors
!
double precision, dimension(dof,dof) :: Mass, invMass
double precision, dimension(dof,1) :: U,Uold,Uolder,stt_R,X
integer i_var,j_var, q_var, AllocateStatus, n, info, time_step_nb
double precision time,delta_t
double precision, dimension(dof,1) :: work
integer, dimension(dof,1) :: ipiv
!
! Declare vectors storing quadrature potins and weights
!
double precision, dimension(p+1,1) :: xq, wq
!
! Output parameters set by the user
!
print *, "---------------------"
print *, "Parameters used in the simulation:"
print *, "Length of the 1-D domain:", L, "m"
print *, "Number of elements:",nb_cell
print *, "CFL number:",cfl
!
! Loop to set the initial values
!
print *, "---------------------"
print *, "Set the initial conditions"
do i_var = 1,dof,1
   X(i_var,1)=(i_var-1)*delta_x
   if (X(i_var,1)<x_step) then 
      U(i_var,1)=u_left
   else
      U(i_var,1)=u_right
   endif
   X(i_var,1)=(i_var-1)*delta_x
end do
Uold=U
Uolder=U
print *, "Done setting the initial conditions"
!
! Call the subroutine to get the quadrature points and weights
!
call quad_pts_and_weight(p,xq,wq)
!
! Compute mass matrix
!
print *, "---------------------"
print *, "Compute the mass matrix and its invert"
Mass(:,:)=0.
call mass_matrix(p,nb_cell,xq,wq,Jwx,Mass)
!
! Invert mass matrix with lapack
!
invMass = Mass
n = size(Mass,1)
call DGETRF(n, n, invMass, n, ipiv, info)
if (info /= 0) then
     stop 'Matrix is singular!'
end if
call DGETRI(n, invMass, n, ipiv, work, n, info)
if (info /= 0) then
     stop 'Matrix inversion failed!'
end if
print *, "Done with computing mass matrix and its invert"
!
! Loop over time to get the numerical solution at each time step
!
print *, "---------------------"
print *,""
print *, "Start loop over time:"
time=0.d0
time_step_nb=1
do while (time<end_time)
   print *, "time step number=",time_step_nb
   print *, "    time =",time,"seconds"
   call stt_residual(p,nb_cell,xq,wq,Jwx,U,stt_R)
   call apply_bc(dof,u_left_bc,u_right_bc,"bc_type",U,stt_R)
   call delta_t_from_cfl(delta_t,delta_x,Uold,cfl,dof)
   print *, "    time step =",delta_t,"seconds"
   print *, ""
   call compute_next_time_step(U,Uold, invMass, stt_R, delta_t,dof)
   Uolder=Uold
   Uold=U
   time=time+delta_t
   time_step_nb=time_step_nb+1
!  do i_var=1,dof,1
!      print *, "U(",i_var,",1)=",Uold(i_var,1)
!   end do
end do
print *, "End of loop over time."
print *, "---------------------"
end program array
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
!!-----------------!!
!!-- End program --!!
!!-----------------!!
!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Start declaring functions !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! test function
!
function phi(x,i)
implicit none
double precision phi
double precision, intent(in) :: x
integer, intent(in) :: i
if ( i == 1 ) then
phi=0.5*(1-x)
else if (i==2) then
phi=0.5*(1+x)
else
print *, "Wrong integer"
stop
end if
end function phi
!
! grad test function
!
function grad_phi(x,i)
implicit none
double precision grad_phi
double precision, intent(in) :: x
integer, intent(in) :: i
if (i==1) then
   grad_phi=-0.5
else if (i==2) then
   grad_phi=0.5
else
   print *, "Wrong integer"
end if
end function grad_phi
!
! flux function
!
function flux(u,flux_type)
implicit none
double precision flux
double precision, intent(in) :: u
character(len=6), intent(in) :: flux_type
select case (flux_type)
case ("burger")
flux=0.5*u*u
case ("linear")
flux=u
case default
print *, "Error: the flux type provided in function 'flux' is not implemented"
stop
end select
end function flux
!
! dissipative function
!
function viscous_flux(grad_u,kappa)
implicit none
double precision viscous_flux
double precision, intent(in) :: grad_u,kappa
viscous_flux=kappa*grad_u
end function viscous_flux
!
! Compute maximum eigenvalue
!
function viscous_coeff(uq,delta_x,dof)
implicit none
integer, intent(in) :: dof
double precision viscous_coeff
double precision, intent(in) :: uq,delta_x
double precision viscous_flux
viscous_flux=0.5*delta_x*abs(uq)
end function viscous_coeff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Done declaring functions !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Start declaring subroutines !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! weight and quadrature points
!
subroutine quad_pts_and_weight(p,xq,wq)
implicit none
integer, intent(in) :: p
double precision, dimension(p+1,1), intent(out) :: xq, wq
integer i_var
do i_var=1,p+1,1
   wq(i_var,1)=1
end do
xq(1,1)=-1./sqrt(3.)
xq(2,1)=1./sqrt(3.)
end subroutine
!
! mass matrix
!
subroutine mass_matrix(p,nb_cell,xq,wq,Jwx,M)
implicit none
integer, intent(in) :: p,nb_cell
double precision, intent(in) :: Jwx
double precision, dimension(p+1,1), intent(in) :: xq,wq
double precision, dimension(p*nb_cell+1,p*nb_cell+1), intent(out) :: M
integer row_var,i_var,j_var,q_var,dof
double precision sum,phi
dof=p*nb_cell+1
do row_var=1,nb_cell,p
   do j_var=1,p+1,1
      do i_var=1,p+1,1
         sum=0.
         do q_var=1,p+1,1
            sum=sum + wq(q_var,1)*phi(xq(q_var,1),i_var)*phi(xq(q_var,1),j_var)*Jwx
         end do
         M(row_var+i_var-1,row_var+j_var-1)=M(row_var+i_var-1,row_var+j_var-1)+sum
      end do
   end do
end do
end subroutine mass_matrix
!
! steady-state residual
!
subroutine stt_residual(p,nb_cell,xq,wq,Jwx,U,stt_R)
implicit none
integer, intent(in) :: p,nb_cell
double precision, intent(in) :: Jwx
double precision, dimension (p+1,1), intent(in) :: xq,wq
double precision, dimension (p*nb_cell+1,1), intent(in) :: U
double precision, dimension (p*nb_cell+1,1), intent(out) :: stt_R
integer row_var,i_var,j_var,q_var,dof
double precision uq,grad_uq,sum,phi,grad_phi,flux,viscous_flux,visc_coeff,viscous_coeff
dof=p*nb_cell+1
stt_R(:,:)=0.
do row_var=1,nb_cell,1
   do j_var=1,p+1,1
      sum=0.
      do q_var=1,p+1,1
         uq=0.
         grad_uq=0.
         do i_var=1,p+1,1
            uq=uq+phi(xq(q_var,1),i_var)*U(p*row_var+i_var-1,1)
            grad_uq=grad_uq+grad_phi(xq(q_var,1),i_var)*U(p*row_var+i_var-1,1)
         end do
         visc_coeff=viscous_coeff(uq,2*Jwx,dof)
         sum=sum+grad_phi(xq(q_var,1),j_var)*flux(uq,"burger")-viscous_flux(grad_uq,0d0)*grad_phi(xq(q_var,1),j_var)/Jwx
      end do
      stt_R(p*row_var+j_var-1,1)=stt_R(p*row_var+j_var-1,1)+sum
   end do
end do
end subroutine stt_residual
!
! Bounday conditions
!
subroutine apply_bc(dof,u_left_bc,u_right_bc,bc_type,U,stt_R)
implicit none
integer, intent(in) :: dof
double precision, intent(in) :: u_left_bc,u_right_bc
character(len=7), intent(in) :: bc_type 
double precision, dimension(dof,1), intent(in) :: U
double precision, dimension(dof,1), intent(inout) :: stt_R
double precision phi, flux
stt_R(1,1)=stt_R(1,1)+flux(u_left_bc,"burger")*phi(-1.d0,1) ! left bc
stt_R(dof,1)=stt_R(dof,1)-flux(u_right_bc,"burger")*phi(1.d0,2) ! right bc
end subroutine apply_bc
!
! Compute the numerical solution at each time step
!
subroutine compute_next_time_step(U,U_old, invMass, stt_residual, delta_t,dof)
implicit none
integer, intent(in) :: dof
double precision, dimension(dof,1), intent(in) :: U_old
double precision, dimension(dof,1), intent(out) :: U
double precision, dimension(dof,1), intent(in) :: stt_residual
double precision, dimension(dof,dof), intent(in) :: invMass
double precision, intent(in) :: delta_t
double precision, dimension(dof,1) :: A
U=matmul(invMass,stt_residual)
U=U*delta_t
U=U+U_old
end subroutine compute_next_time_step
!
! Compute time step based on cfl number
!
subroutine delta_t_from_cfl(delta_t,delta_x,U,cfl,dof)
implicit none
integer, intent(in) :: dof
double precision, dimension(dof,1), intent(in) :: U
double precision, intent(in) :: cfl,delta_x
double precision, intent(out) :: delta_t
delta_t=cfl*delta_x/maxval(U)
end subroutine delta_t_from_cfl
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Done declaring subroutines !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
