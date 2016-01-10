MODULE global_data
implicit real*8(a-h,l-z)
!     implicit integer*8(i,j,k)
!     implicit none
!
!   This global_data MODULE is for the DEM_1D code which follows
!
!     im  = total number of mesh points
!     ieq = number of equations to solve in each fluid:
!              1- alpha, 2- rho, 3- rho*u, 4- rho*E)
!     ifl = number of fluids (restricted to 2 in the present version)
      integer(kind=4), parameter :: im=102, ieq=4, ifl=2
!
      integer(kind=8), parameter :: neq = ifl*ieq*(im-2)
!
!      integer(kind=4), parameter :: ONE=1, TWO=2
!
!     Cell centers array: distance coordiate for cell centers, xx(im)
      double precision :: xx(im)
!
!     Geometry arrays:
!       vol        = cell volumes
!       surfent    = left cross section of a cell
!       surfsort   = right cross section of a cell
!       surffluide = cross section available for the flow at cell boundaries
      double precision :: vol(im), surfent(im), surfsort(im), surffluide(im-1)
!
!     Conservative variables array:
!        {alpha, alpha*rho, alpha*rho*u, alpha*rho*E}
      double precision :: uc(ieq,im,ifl)
!
!     Cell edge flux arrays
!       flag     =  Lagrangian fluxes
!       fkeuleur =  Eulerian fluxes for phase k
      double precision :: flag(ieq,im,ifl,ifl), fkeuler(ieq,im,ifl,ifl)

!      double precision :: rhs(ieq,im,ifl)
!
!     Stiffened gas EOS coefficient arrays
      double precision :: gammaf(ifl), pinff(ifl)
      double precision :: qf(ifl), qpf(ifl)
!
!     contact : contact surfaces between fluids at cell boundaries
!     xstar   : phase function at cell boundaries
      double precision :: contact(im,ifl,ifl), xstar(ifl,im,ifl,ifl)
!
!     sautplus : jump in characteristic function for cell on right side of face
!     sautmoins : jump in characteristic function for cell on left side of face
      double precision :: sautplus(ifl,im,ifl,ifl), sautmoins(ifl,im,ifl,ifl)
!
!     Relaxation arrays:
!       zk     :   acoustic impedance for phase k
!       pk     :   pressure of phase k
!       alphak :   volume fraction of phase k
!       uk     :   velocity of phase k
      double precision :: zk(ifl),pk(ifl),alphak(ifl),uk(ifl)
!
!     Source terms
      double precision :: source(ifl,ieq)
!
!     Initial data in primitive variables:
!       prim0l = initial on the left
!       prim0r = initial on the right
      double precision :: prim0l(ieq,ifl), prim0r(ieq,ifl)
!
!     Primitive variables array:
!       {alpha, rho, u, p}
      double precision :: prim(ieq,im,ifl)
!
!     Arrays used locally by the solver
      double precision :: uur(ieq), uul(ieq), fr(ieq), fl(ieq), fsol(ieq)
!
!     Mixture variables (used for post processing only)
      double precision :: romel(im),umel(im),pmel(im)
!
!~~~~~~~~~~~~~~~~~~~~ INPUT DATA  ~~~~~~~~~~~~~~~~~~~~~
!     kmax   : number of time steps
!     kprint : printing frequency
      integer(kind=8), parameter:: kmax=200000
      integer(kind=8) :: kprint=10000
!
!     initial time step size:
!                    subsequently calculated in subroutine Use_BruteForce_Method
!     double precision :: dt=0.04
      double precision :: dt=1.d-7
!
!     CFL used for the computations
!       Normally cfl is set to 0.5d0, but may need to
!       be smaller if relaxation is turned on (irelax=1).
      double precision :: cfl=0.5d0
!
!     irelax=0 --> no relaxation is used
!     irelax=1 --> finite rate relaxation is used
      integer(kind=4), parameter:: irelax=0
!
!     Water-vapor stiffened gas EOS parameters; double precision
      double precision :: gamma1=2.35d0,pinf1=1.d9,cv1=1816.d0  ! water
      double precision :: gamma2=1.43d0,pinf2=0.d0,cv2=1040.d0  ! vapor
      double precision :: q1 = -1167.0d3, qp1 =   0.0d0         ! water
      double precision :: q2 =  2030.0d3, qp2 = -23.0d3         ! vapor
!
!     Initialization of a few constants:
!       temps = time; updated each time step in subroutine Use_BruteForce_Method
!       cmax  = maximum wave speed in the entire domain (used for CFL);
!               updated each time step in subroutine Derv
      double precision :: cmax=0.d0,temps=0.d0
!
!~~~~~~~~~~~~~~~~~~ GEOMETRY ~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     tube length
      double precision :: ltube=2.d0
!
!     dx= constant cell spacing, set in subroutine Init
      double precision :: dx
!
!     Indexes for the geometric junctions: 
      integer(kind=4), parameter:: idisc1=26 ! index at which the first segment ends
      integer(kind=4), parameter:: idisc2=51 ! index at which the second segment ends
      integer(kind=4), parameter:: idisc3=76 ! index at which the third segment ends
      integer(kind=4), parameter:: idisc4=102 ! index at which the fourth segment ends
!
!     Inlet and outlet cross sections for the each segment, 1 to 4:
!       sent = inlet cross-sectional area
!       ssort= outlet cross-sectional area
      double precision :: sent1=0.2d0,  ssort1=0.2d0
      double precision :: sent2=0.2d0,  ssort2=0.2d0
      double precision :: sent3=0.2d0,  ssort3=0.2d0
      double precision :: sent4=0.2d0,  ssort4=0.2d0

END MODULE global_data

!=======================================================================
!======================= DEM_1D ========================================
!=======================================================================
!
PROGRAM DEM_1D
use global_data
implicit real*8(a-h,l-z)

!
!  Code for solving 1-D two phase flows in ducts of variable
!  cross-section, including mass transfer, using the Discrete
!  Equations Method (DEM).
!
!  Reference: Ray A. Berry, Richard Saurel, Olivier LeMetayer,
!	      "The discrete equation method (DEM) for fully compressible,
!             two-phase flows in ducts of spatially varying cross-section,"
!             Nuclear Engineering and Design 240 (2010) 3797-3818.	
!
!  Coded by:  Ray Berry, Idaho National Laboratory
!             Richard Saurel, Polytech'Marseille
!             Tammi Grimmett, Idaho National Laboratory
!	
!
!    x=0.                                                      x=xmax
!    /|<----------------- problem domain ------------------------>|/
!    /|                                                           |/
!    /|      |        |      |      |      |        |      |      |/
! 1  /|   2  |  ::::  | i-1  |  i   | i+1  |  ::::  | im-2 | im-1 |/ im
!    /|      |        |      |      |      |        |      |      |/
!     1      2       i-2    i-1     i     i+1            im-2   im-1
!
!             cell indexing
!          cell-face indexing
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Initialization
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call init
	k = 0
	call Print_Results(k)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Explicit method
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    call Use_BruteForce_method

stop
END PROGRAM DEM_1d
!
!=======================================================================
!======================= Use_BruteForce_method =========================
!=======================================================================

SUBROUTINE Use_BruteForce_method
use global_data
implicit real*8(a-h,l-z)

   double precision :: rhs(ieq,im,ifl)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Time loop begins
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
   do k=1,kmax
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Initialize rhs array
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO ifld = 1,ifl
         DO icell = 1,im
            DO ivar = 1,ieq
               rhs(ivar,icell,ifld) = 0.0d0
            END DO
         END DO
      END DO

      call Derv(temps, rhs)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Do time update
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=2,im-1
         do iu=1,ieq
            do kf=1,ifl
               uc(iu,i,kf) = uc(iu,i,kf) + dt*rhs(iu,i,kf)
            enddo
         enddo
      enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  Correction of volume fraction round off errors
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   do i=1,im
      sum=0.d0
      do kf=1,ifl  ! sum of volume fractions
         sum = sum+uc(1,i,kf)
      enddo

      do kf=1,ifl
         uc(1,i,kf) = uc(1,i,kf)/sum   ! TO CHECK.  Normalize
      enddo
   enddo

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Computation of primitive variables
!    (important for varying cross sections)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     call UpdatePrimFromUC(ieq,im,ifl,gammaf,pinff,qf,uc,prim)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Time step calculation
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     dt=cfl*dx/cmax
     temps=temps+dt

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!     Print variables
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (mod(k,kprint).eq.0)then
         print*,'iter',k, 'time',temps, 'dt',dt, 'cmax', cmax
         call Print_Results(k)
      endif
!	  if (k.gt.142800)then
!		 print*,'iter',k, 'time',temps, 'dt',dt, 'cmax', cmax
!      endif 
!	  if (k.eq.142903)then
!		  call Print_Results(k)
!	  endif

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!  End of the time loop
   end do
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

   return
END SUBROUTINE Use_BruteForce_method
!
!=======================================================================
!======================= Derv ==========================================
!=======================================================================
!
SUBROUTINE Derv(t,rhs)
use global_data
implicit real*8(a-h,l-z)

   double precision :: t
   double precision :: rhs(ieq,im,ifl)
   double precision :: gammaL, gammaR
   double precision :: pinfL, pinfR
   double precision :: qL, qR
!
!   BOUNDARY CONDITIONS
!
!    Absorption BC
!        do kf=1,ifl
!        do iu=1,ieq
!        uc(iu,1,kf)=uc(iu,2,kf)
!        uc(iu,im,kf)=uc(iu,im-1,kf)
!        enddo
!        enddo
!
!**********************************************************
!
!   Closed left boundary, ustar=0, at left edge of cell 2
!
      goto 10
      kfL = 1
      kfR = 1
      call closed_LB(ieq,im,ifl,kfL,kfR,t,prim,gamma1,pinf1,  &
                        q1,pstar,fkeuler)
!
      kfL = 2
      kfR = 2
      call closed_LB(ieq,im,ifl,kfL,kfR,t,prim,gamma2,pinf2,  &
                        q2,pstar,fkeuler)
!
      kfL = 1
      kfR = 2
      call closed_LB(ieq,im,ifl,kfL,kfR,t,prim,gamma2,pinf2,  &
                        q2,pstar,fkeuler)
      flag(1,1,1,2) = 0.d0
      flag(2,1,1,2) = 0.d0
      flag(3,1,1,2) = pstar
      flag(4,1,1,2) = 0.d0
!
      kfL = 2
      kfR = 1
      call closed_LB(ieq,im,ifl,kfL,kfR,t,prim,gamma1,pinf1,  &
                        q1,pstar,fkeuler)
      flag(1,1,2,1) = 0.d0
      flag(2,1,2,1) = 0.d0
      flag(3,1,2,1) = pstar
      flag(4,1,2,1) = 0.d0
!
!     set volume fractions
      uc(1,1,1) = uc(1,2,1)
      uc(1,1,2) = uc(1,2,2)
   10 continue
!
!   Prescribed pressure at the INLET
!      goto 20
         kfL = 1
         kfR = 1
         gammaL= gamma1
         gammaR= gamma1
         pinfL = pinf1
         pinfR = pinf1
         qL    = q1
         qR    = q1
!        rho0=1361.d0  ! tank density
         rho0=1000.d0  ! tank density
         p0  =10.d5    ! tank pressure
!      get ustar, pstar, Eulerian flux at inlet face with 1,1 arrangement
         call inlet_BC(ieq,im,ifl,kfL,kfR,rho0,p0,t,prim,  &
                          gammaL,gammaR,pinfL,pinfR,qL,qR, &
                          pr,zr,ur,ustar,pstar,fkeuler)
!        print*,'ustar LIQ 11',ustar, pstar, t

!*****************************
         kfL = 2
         kfR = 2
         gammaL= gamma2
         gammaR= gamma2
         pinfL = pinf2
         pinfR = pinf2
         qL    = q2
         qR    = q2
!         rho0=1.d0    ! tank density
         rho0=5.4766d0    ! tank density
         p0  =10.d5   ! tank pressure
!     get ustar, pstar, Eulerian flux at inlet face with 2,2 arrangement
         call inlet_BC(ieq,im,ifl,kfL,kfR,rho0,p0,t,prim,  &
                          gammaL,gammaR,pinfL,pinfR,qL,qR, &
                          pr,zr,ur,ustar,pstar,fkeuler)
!        print*,'ustar GAS 22',ustar, pstar

!*****************************
!
!       CONTACT  12
!   Prescribed pressure at the INLET
         kfL = 1
         kfR = 2
         gammaL= gamma1
         gammaR= gamma2
         pinfL = pinf1
         pinfR = pinf2
         qL    = q1
         qR    = q2
!        rho0=1361.d0  ! tank density
         rho0=1000.d0  ! tank density
         p0  =10.d5    ! tank pressure
!     get ustar, pstar, Eulerian flux at inlet face with 1,2 arrangement
         call inlet_BC(ieq,im,ifl,kfL,kfR,rho0,p0,t,prim,  &
                          gammaL,gammaR,pinfL,pinfR,qL,qR, &
                          pr,zr,ur,ustar,pstar,fkeuler)
!        print*,'ustar LIQ 12',ustar
!
!     get Lagrangian flux at inlet face with 1,2 arrangement
         flag(1,1,1,2) = -ustar
         flag(2,1,1,2) = 0.d0
         flag(3,1,1,2) = pstar
         flag(4,1,1,2) = pstar*ustar

!*****************************
! 
!       CONTACT  21
!   Prescribed pressure at the INLET
         kfL = 2
         kfR = 1
         gammaL= gamma2
         gammaR= gamma1
         pinfL = pinf2
         pinfR = pinf1
         qL    = q2
         qR    = q1
!         rho0=1.d0     ! tank density
         rho0=5.4766d0     ! tank density
         p0  =10.d5    ! tank pressure
!     get ustar, pstar, Eulerian flux at inlet face with 2,1 arrangement
         call inlet_BC(ieq,im,ifl,kfL,kfR,rho0,p0,t,prim,  &
                          gammaL,gammaR,pinfL,pinfR,qL,qR, &
                          pr,zr,ur,ustar,pstar,fkeuler)
!        print*,'ustar GAS 21',ustar
!
!     get Lagrangian flux at inlet face with 2,1 arrangement
         flag(1,1,2,1) = -ustar
         flag(2,1,2,1) = 0.d0
         flag(3,1,2,1) = pstar
         flag(4,1,2,1) = pstar*ustar

!*****************************
!    Set inlet boundary volume fractions
          uc(1,1,1) = 0.5d0
          uc(1,1,2) = 0.5d0
          prim(1,1,1) = 0.5d0
          prim(2,1,2) = 0.5d0
!
   20 continue
!**********************************************************
!
!   Closed right boundary, ustar=0, at right edge of cell im-1
!
      goto 30
      kfL = 1
      kfR = 1
      call closed_RB(ieq,im,ifl,kfL,kfR,t,prim,gamma1,pinf1,  &
                        q1,pstar,fkeuler)
!
      kfL = 2
      kfR = 2
      call closed_RB(ieq,im,ifl,kfL,kfR,t,prim,gamma2,pinf2,  &
                        q2,pstar,fkeuler)
!
      kfL = 1
      kfR = 2
	  j = im-1
      call closed_RB(ieq,im,ifl,kfL,kfR,t,prim,gamma1,pinf1,  &
                        q1,pstar,fkeuler)
      flag(1,j,1,2) = 0.d0
      flag(2,j,1,2) = 0.d0
      flag(3,j,1,2) = pstar
      flag(4,j,1,2) = 0.d0
!
      kfL = 2
      kfR = 1
	  j = im-1
      call closed_RB(ieq,im,ifl,kfL,kfR,t,prim,gamma2,pinf2,  &
                        q2,pstar,fkeuler)
      flag(1,j,2,1) = 0.d0
      flag(2,j,2,1) = 0.d0
      flag(3,j,2,1) = pstar
      flag(4,j,2,1) = 0.d0
!
!     set volume fractions
      uc(1,im,1) = uc(1,im-1,1)
      uc(1,im,2) = uc(1,im-1,2)
   30 continue
!
!   Prescribed pressure at the OUTLET
!      goto 40
!       liquid
         kfL = 1
         kfR = 1
         pstar = 5.d5
!         pstar = 9.7d5
!   get rhostar, ustar, estar, Eulerian flux at outlet face with 1,1 arrangement
         call outlet_BC(ieq,im,ifl,kfL,kfR,prim,gamma1,pinf1,q1,pstar, &
                        rhostar,ustar,estar,fkeuler)

!----------
!       gas
         kfL = 2
         kfR = 2
         pstar = 5.d5
!         pstar = 9.7d5
!   get rhostar, ustar, estar, Eulerian flux at outlet face with 2,2 arrangement
         call outlet_BC(ieq,im,ifl,kfL,kfR,prim,gamma2,pinf2,q2,pstar, &
                        rhostar,ustar,estar,fkeuler)
!
!----------
!       Contact 12
!   Prescribed pressure at the OUTLET
         kfL = 1
         kfR = 2
         j = im-1
         pstar = 5.d5
!         pstar = 9.7d5
!   get rhostar, ustar, estar, Eulerian flux at outlet face with 1,2 arrangement
         call outlet_BC(ieq,im,ifl,kfL,kfR,prim,gamma1,pinf1,q1,pstar, &
                        rhostar,ustar,estar,fkeuler)
!        print*,'ustar LIQ 12',ustar
!
!   get Lagrangian flux at outlet face with 1,2 arrangement
         flag(1,j,1,2) = -ustar
         flag(2,j,1,2) = 0.d0
         flag(3,j,1,2) = pstar
         flag(4,j,1,2) = pstar*ustar
! 
!----------
!       Contact 21
         kfL = 2
         kfR = 1
         j = im-1
         pstar = 5.d5
!         pstar = 9.7d5
!   get rhostar, ustar, estar, Eulerian flux at outlet face with 2,1 arrangement
         call outlet_BC(ieq,im,ifl,kfL,kfR,prim,gamma2,pinf2,q2,pstar, &
                        rhostar,ustar,estar,fkeuler)
!
!   get Lagrangian flux at outlet face with 2,1 arrangement
         flag(1,j,2,1) = -ustar
         flag(2,j,2,1) = 0.d0
         flag(3,j,2,1) = pstar
         flag(4,j,2,1) = pstar*ustar
!
!----------
!    Set outlet boundary volume fractions (not used)
          uc(1,im,1) = uc(1,im-1,1)
          uc(1,im,2) = uc(1,im-1,2)
          prim(1,im,1) = prim(1,im-1,1)
          prim(1,im,2) = prim(1,im-1,2)
!
   40 continue
!*******************************************************
!
!  Maximum wave speed reset to zero
        cmax=0.d0
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Computation of Eulerian and Lagrangian fluxes
!      --- d means right
!      --- g means left
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! Loop over each internal cell face and compute fluxes
       do i=2,im-2
!       Loop over each contact pair
        do kfd=1,ifl ! kfd= index for the fluid on right
         do kfg=1,ifl ! kfg= index for the fluid on left
!               Store right and left conservative variables and fluxes
!               of the Euler equations on each side of the cell boundary
                do iu=1,ieq
                 uur(iu) = uc(iu,i+1,kfd)/uc(1,i+1,kfd)
                enddo
                gammaR= gammaf(kfd)
                pinfR = pinff(kfd)
                qR    = qf(kfd)
                call calcflux(uur,fr,cr,gammaR,pinfR,qR)
!
                do iu=1,ieq
                 uul(iu) = uc(iu,i,kfg)/uc(1,i,kfg)
                enddo
                gammaL= gammaf(kfg)
                pinfL = pinff(kfg)
                qL    = qf(kfg)
                call calcflux(uul,fl,cl,gammaL,pinfL,qL)
!               Euler flux
!
!               velocities on right and left
                ur = uc(3,i+1,kfd)/uc(2,i+1,kfd)
                ul = uc(3,i,kfg)/uc(2,i,kfg)
!
!             Lagrangian fluxes computation
              call hllclag(uur,uul,fr,fl,ur,cr,ul,cl,fsol, &
                           gammaR,gammaL,pinfR,pinfL,qR,qL)
!
              do iu=1,ieq
               flag(iu,i,kfg,kfd) = fsol(iu)
              enddo
!
!             Eulerian fluxes computation
              call hllc(uur,uul,fr,fl,ur,cr,ul,cl,fsol, &
                        gammaR,gammaL,pinfR,pinfL,qR,qL)
!
              do iu=1,ieq
               fkeuler(iu,i,kfg,kfd) = fsol(iu)
              enddo
!
!          Update of the maximum wave speed (CFL)
!           if(t.ge.0.74547148556480625)then
!      print*, 'i',i, 'kfg',kfg, 'kfd',kfd
!	  print*, 'alphal',uc(1,i,kfg), 'alphar',uc(1,i+1,kfd)
!		   endif
           cmax=dmax1(cmax,dabs(ur)+cr,dabs(ul)+cl)
         enddo
        enddo
       enddo
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       Computation of the indicator functions at cell boundaries (faces)
!       and corresponding jumps
!          saut  - jump in xstar
!          plus  - plus
!          moins - minus
!       So  sautplus  - jump for cell on plus side of face
!           sautmoins - jump for cell on minus side of face
!      DO NOT MODIFY THIS BLOCK OF CODE!!!	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       do i=1,im-1 ! loop over cell faces
           do kf=1,ifl ! kf = index for a given fluid
              do kfd=1,ifl ! loop for the 4 possible pairs
              do kfg=1,ifl
                if(kf.ne.kfg.and.kf.ne.kfd)then
                 xstar(kf,i,kfg,kfd)    = 0.d0
                 sautplus(kf,i,kfg,kfd) = 0.d0
                 sautmoins(kf,i,kfg,kfd)= 0.d0
                else
                    if(kfg.eq.kfd)then
                        xstar(kf,i,kf,kf)    = 1.d0
                        sautplus(kf,i,kf,kf) = 0.d0
                        sautmoins(kf,i,kf,kf)= 0.d0
                    elseif(kf.eq.kfg)then
                      if((-flag(1,i,kfg,kfd)).gt.0.d0)then
                        xstar(kf,i,kfg,kfd)    = 1.d0
                        sautplus(kf,i,kfg,kfd) =-1.d0
                        sautmoins(kf,i,kfg,kfd)= 0.d0
                      else
                        xstar(kf,i,kfg,kfd)    = 0.d0
                        sautplus(kf,i,kfg,kfd) = 0.d0
                        sautmoins(kf,i,kfg,kfd)=-1.d0
                      endif
                    elseif(kf.eq.kfd)then
                      if((-flag(1,i,kfg,kfd)).lt.0.d0)then
                       xstar(kf,i,kfg,kfd)    = 1.d0
                       sautplus(kf,i,kfg,kfd) = 0.d0
                       sautmoins(kf,i,kfg,kfd)= 1.d0
                      else
                       xstar(kf,i,kfg,kfd)    = 0.d0
                       sautplus(kf,i,kfg,kfd) = 1.d0
                       sautmoins(kf,i,kfg,kfd)= 0.d0
                      endif
                    endif
                  endif
                enddo
              enddo
            enddo
       enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!       Computation of normalized contact lengths
!          --- d means right
!          --- g means left	
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
         do i=1,im-1 ! loop over cell faces
           alpha1g = uc(1,i,1)
           alpha1d = uc(1,i+1,1)
!
           alpha2g = uc(1,i,2)
           alpha2d = uc(1,i+1,2)
!
           contact(i,1,1) = dmin1(alpha1g,alpha1d)
           contact(i,2,2) = dmin1(alpha2g,alpha2d)
!
           alpha1g = alpha1g-contact(i,1,1)
           alpha2g = alpha2g-contact(i,2,2)
!
           alpha1d = alpha1d-contact(i,1,1)
           alpha2d = alpha2d-contact(i,2,2)
!
           contact(i,1,2) = dmin1(alpha1g,alpha2d)
           contact(i,2,1) = dmin1(alpha2g,alpha1d)
         enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Computation of the two phase flux and RHS (time) update
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do i=2,im-1 ! loop over all cell centers
       do kf=1,ifl ! fluid index
        do iu=1,ieq ! equation index
!
! signalplus - Eulerian (or conservative) two phase flux for right cell boundary
! signalmoins - Eulerian (or conservative) two phase flux for left cell boundary
          signalplus = 0.d0
          signalmoins= 0.d0
!
!Ray's convention, for convention modification below:
! signalplusnc - Lagrangian (nonconservative) two phase flux for right cell boundary
! signalmoinsnc - Lagrangian (nonconservative) two phase flux for left cell boundary
!
          signalplusnc = 0.d0
          signalmoinsnc= 0.d0
        do kfd=1,ifl
        do kfg=1,ifl
         signalplus = signalplus+xstar(kf,i,kfg,kfd) &
          *surffluide(i)*contact(i,kfg,kfd)*fkeuler(iu,i,kfg,kfd)
         signalmoins= signalmoins+xstar(kf,i-1,kfg,kfd) &
          *surffluide(i-1)*contact(i-1,kfg,kfd)*fkeuler(iu,i-1,kfg,kfd)
!
!  Original (wrong!)
!        signalplusnc = signalplusnc+sautplus(kf,i-1,kfg,kfd) &
!         *surffluide(i)*contact(i-1,kfg,kfd)*flag(iu,i-1,kfg,kfd)
!        signalmoinsnc= signalmoinsnc+sautmoins(kf,i,kfg,kfd) &
!         *surffluide(i-1)*contact(i,kfg,kfd)*flag(iu,i,kfg,kfd)
!  This is wrong because the surffluide indices need to match the indices
!    of the other terms in each expression.  This will make it correct
!    because the signalplusnc and signalmoinsnc are summed, but the 
!    convention is confusing because now signalplusnc is the Lagrangian flux
!    at the left face of the cell and signalmoinsnc is the Lagrangian flux
!    at the rigt face of the cell.
!
!  Ray's modification (correct but inconsistent with Richard's naming convention)
        signalplusnc = signalplusnc+sautmoins(kf,i,kfg,kfd) &
         *surffluide(i)*contact(i,kfg,kfd)*flag(iu,i,kfg,kfd)
        signalmoinsnc= signalmoinsnc+sautplus(kf,i-1,kfg,kfd) &
         *surffluide(i-1)*contact(i-1,kfg,kfd)*flag(iu,i-1,kfg,kfd)
!
        enddo
        enddo
!!$            uc(iu,i,kf) = uc(iu,i,kf) &
!!$                        -dt/vol(i)*(signalplus-signalmoins) &
!!$                        +dt/vol(i)*(signalplusnc+signalmoinsnc)
!
            rhs(iu,i,kf) = -1.0d0/vol(i)*(signalplus-signalmoins) &
                           +1.0d0/vol(i)*(signalplusnc+signalmoinsnc)
!
!            uc(iu,i,kf) = uc(iu,i,kf)+dt*rhs(iu,i,kf)
!
         alpha= prim(1,i,kf)
         p    = prim(4,i,kf) 
         if(iu.eq.3)then
!!$          uc(iu,i,kf) = uc(iu,i,kf) &
!!$                  +dt/vol(i)*p*alpha*(surffluide(i)-surffluide(i-1))
          rhs(iu,i,kf) = rhs(iu,i,kf) &
                  +1.0d0/vol(i)*p*alpha*(surffluide(i)-surffluide(i-1))
         endif
        enddo
       enddo
      enddo
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Relaxation terms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
       if(irelax.eq.1)then
         do i=2,im-1
            call Relaxation_Terms(i)
!
            do iu=1,ieq
             do kf=1,ifl
!               uc(iu,i,kf) = uc(iu,i,kf)+dt*source(kf,iu)
               rhs(iu,i,kf) = rhs(iu,i,kf)+source(kf,iu)
!               uc(iu,i,kf) = uc(iu,i,kf)+dt*rhs(iu,i,kf)
             enddo
            enddo

         enddo
       endif


      return
END SUBROUTINE Derv

!======================================================================
!======================= Relaxation_Terms =============================
!======================================================================

SUBROUTINE Relaxation_Terms(i)
use global_data
implicit real*8(a-h,l-z)

!  BE CAREFUL NOT TO USE i FOR A LOOP INDEX IN THIS SUBROUTINE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Relaxation terms
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      do kf=1,ifl
         alpha= uc(1,i,kf)
         rho  = uc(2,i,kf)/alpha
         u    = uc(3,i,kf)/uc(2,i,kf)
         e    = uc(4,i,kf)/uc(2,i,kf)-u**2/2.d0
!
         gamma= gammaf(kf)
         pinf = pinff(kf)
         q    = qf(kf)
         call eospression(rho,e,gamma,pinf,q,p,c)
!
         zk(kf)    = rho*c
         pk(kf)    = p
         alphak(kf)= alpha
         uk(kf)    = u
         do iu=1,ieq
           source(kf,iu) = 0.d0
         enddo
      enddo
!
!   Compute interface pressure and velocity
      pinter = (zk(1)*pk(2)+zk(2)*pk(1))/(zk(1)+zk(2))
      vinter = (zk(1)*uk(1)+zk(2)*uk(2))/(zk(1)+zk(2))
!
!   Compute interfacial specific area (area/volume):	
!    ssvmax is approximately the maximum interfacial area concentration,
!    which is used to scale the following form function:
!      6.75*(alpha**2)*(1-alpha)
!   This function, concave down, with peak of approximately 1 at alpha=2/3
!   which goes to zero at alpha=0 and at alpha=1.
!      ssvmax = 1700.d0
!      ssvmax = 50.d0
      ssvmax = 0.d0
      vert_offset = 0.d0
!      ssv    = ssvmax*6.75d0*alphak(2)*alphak(2)*alphak(1)
	  ssv    = ssvmax*4.d0*alphak(1)*alphak(2) - vert_offset
	  if (ssv.lt.0.d0) ssv=0.d0
!
!   Specify or compute pressure and velocity relaxation coefficients	
      mu     = 1.d0/(zk(1)+zk(2))
      lambda = zk(1)*zk(2)*mu
!   Note: if lambda =0, no velocity relaxation: stratified uniform flow
!
!   Compute phasic temperatures
      kf = 1
      alpha1= uc(1,i,kf)
      rho1  = uc(2,i,kf)/alpha1
      u1    = uc(3,i,kf)/uc(2,i,kf)
      e1    = uc(4,i,kf)/uc(2,i,kf)-u1**2/2.d0
!      gamma1= gammaf(kf)
!      pinf1 = pinff(kf)
!      q1    = qf(kf)
      call eospression(rho1,e1,gamma1,pinf1,q1,p1,c1)
      t1   = (p1+pinf1)/(cv1*(gamma1-1.d0)*rho1)
!
      kf = 2
      alpha2= uc(1,i,kf)
      rho2  = uc(2,i,kf)/alpha2
      u2    = uc(3,i,kf)/uc(2,i,kf)
      e2    = uc(4,i,kf)/uc(2,i,kf)-u2**2/2.d0
!      gamma2= gammaf(kf)
!      pinf2 = pinff(kf)
!      q2    = qf(kf)
      call eospression(rho2,e2,gamma2,pinf2,q2,p2,c2)
      t2   = (p2+pinf2)/(cv2*(gamma2-1.d0)*rho2)
!
!      print*,'t1',t1,'t2',t2
!
!   For relaxation only, no interphase heat or mass transfer
      hconv1 = 0.d0
      hconv2 = 0.d0
      mdot   = 0.d0
      ti = (t1+t2)/2.d0
      rhoi = rho1
      goto 50
!      
!   INTERPHASE HEAT and MASS TRANSFER MODEL
!    Mass transfer model is based on heat transfer at interface	
!
!    Assume the interface pressure is at saturation condition,
!    i.e., Psat=pinter.
!	
!    Compute the interface saturation temperature (ti=Tsat(pinter),
!    phasic enthalpies at Tsat, and heat of vaporization at Tsat	
      call tsat(pinter,ti,h1,h2,lvap)
!      
!    Assume interface density is saturated liquid density, rhol(pinter)      
      rhoi = (pinter+pinf1)/((gamma1-1.d0)*cv1*ti)
!
!    Compute total saturated phasic enthalpies
      ht1 = h1+0.5d0*vinter*vinter
      ht2 = h2+0.5d0*vinter*vinter
!      
!    Set the convective heat transfer coefficients      
      conduct1 = 0.5d0
      rdrop    = 3.d0*alpha1/ssv
      hconv1   = 5.d0*conduct1/rdrop
      conduct2 = 0.026
      visc2    = 134.4d-7
      cp2      = gamma2*cv2
      ddrop    = 2.d0*rdrop
      Re       = rho2*ddrop*dabs(u2-u1)/visc2
      Pr       = visc2*cp2/conduct2
      Nu       = 2.d0 + 0.6d0*(Re**0.5d0)*(Pr**0.33d0)
      hconv2   = conduct2*Nu/ddrop
!
!    If no interfacial heat and mass transfer are desired, set
!    interfacial heat transfer coefficients to zero here, i.e.
!      hconv1 = 0.d0
!      hconv2 = 0.d0
!	
!    Compute the interphase mass flow rate	
      mdot = (hconv1*(t1-ti) + hconv2*(t2-ti))/lvap
!
   50 continue
   
!    Set gravity, friction, and hydraulic diameter
!      grav = 0.d0
      grav = -9.81d1
	  fric = 2.d3
	  pi   = 4.0D0*DATAN(1.0D0)
	  area_ave = 0.5d0*(surffluide(i)+surffluide(i-1))
	  dh = dsqrt(4.d0*area_ave/pi)
!	   lambda = 1.d-1 * lambda
!
!    Compute source terms:	
      source(1,1) =   mu*ssv*(pk(1)-pk(2)) - mdot*ssv/rhoi
      source(1,2) = - mdot*ssv
      source(1,3) =   lambda*ssv*(uk(2)-uk(1)) - mdot*ssv*vinter &
                    + alpha1*rho1*grav &
					- 0.5d0*fric*alpha1*rho1*u1*dabs(u1)/dh  	  
      source(1,4) = - pinter*source(1,1) &
                    + vinter*lambda*ssv*(uk(2)-uk(1)) &
                    + mdot*ssv*((pinter/rhoi)-ht1) &
                    + ssv*hconv1*(ti-t1) &
					+ alpha1*rho1*grav*u1
!
	  source(2,1) = - mu*ssv*(pk(1)-pk(2)) + mdot*ssv/rhoi
	  source(2,2) = + mdot*ssv
	  source(2,3) = - lambda*ssv*(uk(2)-uk(1)) + mdot*ssv*vinter &
	                + alpha2*rho2*grav &
					- 0.5d0*fric*alpha2*rho2*u2*dabs(u2)/dh
      source(2,4) =   pinter*source(1,1) &
                    - vinter*lambda*ssv*(uk(2)-uk(1)) &
                    - mdot*ssv*((pinter/rhoi)-ht2) &
                    + ssv*hconv2*(ti-t2) &
					+ alpha2*rho2*grav*u2
!
return
END SUBROUTINE relaxation_terms

!=======================================================================
!======================= tsat ==========================================
!=======================================================================

SUBROUTINE tsat(p,temp,h1,h2,lvap)
!-----------------------------
!  Computes the saturation temperature for water/vapor for a specified
!    pressure using the stiffened gas equation of state for each phase.
!  Also computes the heat of vaporization at this saturation temperature.
!  Other quantities can be computed and returned as commented at the end.
!-----------------------------
!   input:
!     pressure (of liquid if vaporization, of vapor if condensation), p
!   output:
!     saturation temperature at input pressure, temp
!     heat of vaporization, lvap
use global_data
implicit real*8(a-h,l-z)
!
      cp1   = gamma1*cv1
      cp2   = gamma2*cv2
!
!  Form constants
      E = cp2-cv2
      A = (cp1-cp2+qp2-qp1)/E
      B = (q1-q2)/E
      C = (cp2-cp1)/E
      D = (cp1-cv1)/E
      if ((p+pinf1).le.0.d0.or.(p+pinf2).le.0.d0) then
        print*, 'negative p+pinf?, phase ? '
      endif
!      R = A + dlog( (p+pinf1)**D/(p+pinf2) )
      R= A + D*dlog(p+pinf1) - dlog(p+pinf2)
!
!  Initial guess of temp, Tsat(p), for Newton iteration
      temp = 400.d0
!
      do i=1,100
        ff = B/temp + C*dlog(temp) + R
        df = C/temp - B/temp/temp
        if (dabs(ff).lt. 1.d-6) goto 25
        temp = temp - ff/df
        if (i.eq.100) print*, 'i= ', i
!       print*, 'i= ', i, 'temp =', temp
        if (temp .lt. 0.d0) then
	  print*, 'saturation temperature negative'
	  print*, 'i= ', i, 'temp =', temp
        endif
      enddo
   25 continue
!
!  Compute the heat of vaporization at the saturation temperature
      h1 = gamma1*cv1*temp + q1
      h2 = gamma2*cv2*temp + q2
      lvap = h2 - h1
!  Compute other thermodynamic quantities at saturation temperature
!      v1  = (gamma1-1.d0)*cv1*temp/(p+pinf1)
!      v2  = (gamma2-1.d0)*cv2*temp/(p+pinf2)
!      e1  = (p+gamma1*pinf1)*v1/(gamma1-1.d0) + q1
!      e2  = (p+gamma2*pinf2)*v2/(gamma2-1.d0) + q2
!      s1  = qp1 + cv1*dlog(temp**gamma1/(p+pinf1)**(gamma1-1.d0))
!      s2  = qp2 + cv2*dlog(temp**gamma2/(p+pinf2)**(gamma2-1.d0))
!      g1  = (gamma1*cv1-qp1)*temp + q1 &
!            -cv1*temp*dlog(temp**gamma1/(p+pinf1)**(gamma1-1.d0))
!      g2  = (gamma2*cv2-qp2)*temp + q2 &
!            -cv2*temp*dlog(temp**gamma2/(p+pinf2)**(gamma2-1.d0))
!
return
END SUBROUTINE tsat

!======================================================================
!======================= UpdatePrimFromUC =============================
!======================================================================

SUBROUTINE UpdatePrimFromUC(ieq,im,ifl,gammaf,pinff,qf,uc,prim)
!--------------------------------------------
!    Computation of primitive variables
!    (important for varying cross sections)
!--------------------------------------------
implicit none

!... Dummy variables
      integer :: ieq, im, ifl
      double precision :: uc(ieq,im,ifl)
      double precision :: prim(ieq,im,ifl)
      double precision :: gammaf(ifl)
      double precision :: pinff(ifl)
      double precision :: qf(ifl)

!..  Local variables
      integer :: kf, i
      double precision :: alpha, rho, u, e, p, c
      double precision :: gamma, pinf, q
!
      DO kf=1,ifl
         DO i=2,im-1
            alpha= uc(1,i,kf)
            rho  = uc(2,i,kf)/alpha
            u    = uc(3,i,kf)/uc(2,i,kf)
            e    = uc(4,i,kf)/uc(2,i,kf)-u**2/2.d0
            gamma= gammaf(kf)
            pinf = pinff(kf)
            q    = qf(kf)
            CALL eospression(rho,e,gamma,pinf,q,p,c)
            prim(1,i,kf) = alpha
            prim(2,i,kf) = rho
            prim(3,i,kf) = u
            prim(4,i,kf) = p
         END DO
      END DO

return
END SUBROUTINE UpdatePrimFromUC

!=======================================================================
!======================= calcflux ======================================
!=======================================================================

SUBROUTINE calcflux(uu,f,c,gamma,pinf,q)
implicit real*8(a-h,l-z)
       parameter (ieq=4)
       dimension f(ieq),uu(ieq)
             rho =  uu(2)
             u   =uu(3)/uu(2)
             if(dabs(u).lt.1.d-8)then
               u    = 0.d0
               uu(3)= 0.d0
             endif
             e   = uu(4)/uu(2)-u**2/2.d0
             call eospression(rho,e,gamma,pinf,q,p,c)
!            computation of the pressure by the EOS
!
          f(1) = 0.d0
          f(2) = rho*u
          f(3) = rho*u*u+p
          f(4) = rho*u*(e+u*u/2.d0)+p*u
return
END SUBROUTINE calcflux
!
!=======================================================================
!======================= eospe =========================================
!=======================================================================

SUBROUTINE eospe(rho,p,gamma,pinf,q,e)
implicit real*8(a-h,l-z)
       e = (p+gamma*pinf)/rho/(gamma-1.d0) + q
return
END SUBROUTINE eospe
!
!=======================================================================
!======================= eospression ===================================
!=======================================================================

SUBROUTINE eospression(rho,e,gamma,pinf,q,p,c)
implicit real*8(a-h,l-z)
       p = (gamma-1.d0)*rho*(e-q)-gamma*pinf
       c = dsqrt(gamma*(p+pinf)/rho)
return
END SUBROUTINE eospression
!
!=======================================================================
!======================= hllc ==========================================
!=======================================================================

      SUBROUTINE hllc(uur,uul,fr,fl,ur,cr,ul,cl,fsol, &
                      gammaR,gammaL,pinfR,pinfL,qR,qL)
       implicit real*8(a-h,o-z)
       parameter (ieq=4)
       dimension uur(ieq),uul(ieq),fr(ieq),fl(ieq)
       dimension qq(ieq),uhll(ieq),fsol(ieq)
          sr1 = ur+cr
          sr2 = ur-cr
!
          sl1 = ul+cl
          sl2 = ul-cl
!   Davis
          sl  = dmin1(sl2,sr2)
          sr  = dmax1(sr1,sl1)
!
	  xmachL = dabs(ul)/cl
	  xmachR = dabs(ur)/cr
!
      if(xmachL.ge.0.3d0.or.xmachR.ge.0.3d0) then
!       Use the HLLC method
        if(sl.gt.0.d0)then
          do iu=1,ieq
           fsol(iu) = fl(iu)
          enddo
        elseif(sr.lt.0.d0)then
             do iu=1,ieq
              fsol(iu) = fr(iu)
             enddo
        else
!            Calculate speed Sm
             do iu=1,ieq
              uhll(iu)= (sr*uur(iu)-sl*uul(iu)-(fr(iu)-fl(iu)))/(sr-sl)
             enddo 
             sm = uhll(3)/uhll(2)
!            Calculate the corresponding star states
              if(sm.gt.0.d0)then
                  do iu=1,ieq
                   qq(iu) = fl(iu)-sl*uul(iu)
                  enddo
                  s = sl
              else
                  do iu=1,ieq
                   qq(iu) = fr(iu)-sr*uur(iu)
                  enddo
                  s = sr
              endif
               alphastar= 0.d0
               rostar   = qq(2)/(sm-s)
               pstar    = qq(3)-sm*qq(2)
               estar    = (qq(4)-sm*pstar)/qq(2)
!
               uhll(1) = 1.d0
               uhll(2) = rostar
               uhll(3) = rostar*sm
               uhll(4) = rostar*estar
!
              do iu=1,ieq
               fsol(iu) = qq(iu)+s*uhll(iu)
              enddo
        endif
      else
!       Use low-Mach modified flux method
	    rhoL = uul(2)
	    rhoR = uur(2)
	    eL   = uul(4)/uul(2) - 0.5d0*ul*ul
	    eR   = uur(4)/uur(2) - 0.5d0*ur*ur
	    pL   = (gammaL-1.d0)*rhoL*(eL-qL) - gammaL*pinfL
	    pR   = (gammaR-1.d0)*rhoR*(eR-qR) - gammaR*pinfR
	    pstar= (pL+pR)/2.d0
	    ustar= (rhoL*cl*ul+rhoR*cr*ur+pL-pR)/(rhoL*cl+rhoR*cr)
	    if(ustar.ge.0.d0) then
	      rhostar= rhoL + (pstar-pL)/(cl**2)
!          estar  = uul(4)
          estar  = qL+(pstar+gammaL*pinfL)/rhostar/(gammaL-1.d0)
	    else
	      rhostar= rhoR + (pstar-pR)/(cr**2)
!          estar  = uur(4)
          estar  = qR+(pstar+gammaR*pinfR)/rhostar/(gammaR-1.d0)
	    endif
        fsol(1) = 0.d0
        fsol(2) = rhostar*ustar
        fsol(3) = rhostar*ustar*ustar + pstar
!        fsol(4) = estar*ustar + pstar*ustar
        fsol(4) = rhostar*(estar+0.5d0*ustar**2)*ustar + pstar*ustar
      endif
      return
      END SUBROUTINE hllc

!=======================================================================
!======================= hllclag =======================================
!=======================================================================

      SUBROUTINE hllclag(uur,uul,fr,fl,ur,cr,ul,cl,fsol, &
                         gammaR,gammaL,pinfR,pinfL,qR,qL)
       implicit real*8(a-h,o-z)
       parameter (ieq=4)
       dimension uur(ieq),uul(ieq),fr(ieq),fl(ieq)
       dimension qq(ieq),uhll(ieq),fsol(ieq)
          sr1 = ur+cr
          sr2 = ur-cr
!
          sl1 = ul+cl
          sl2 = ul-cl
!   Davis
          sl  = dmin1(sl2,sr2)
          sr  = dmax1(sr1,sl1)
!
	  xmachL = dabs(ul)/cl
	  xmachR = dabs(ur)/cr
!
      if(xmachL.ge.0.3d0.or.xmachR.ge.0.3d0) then
!       Calculate the speed Sm
          do iu=1,ieq
           uhll(iu) = (sr*uur(iu)-sl*uul(iu)-(fr(iu)-fl(iu)))/(sr-sl)
          enddo 

          sm = uhll(3)/uhll(2)
           do iu=1,ieq
            qq(iu) = fr(iu)-sr*uur(iu)
           enddo
                 pstar   =qq(3)-sm*qq(2)
!
                 fsol(1) =-sm
                 fsol(2) = 0.d0
                 fsol(3) = pstar
                 fsol(4) = pstar*sm
      else
!       Use low-Mach modified flux method
        rhoL = uul(2)
        rhoR = uur(2)
        eL   = uul(4)/uul(2) - 0.5d0*ul*ul
        eR   = uur(4)/uur(2) - 0.5d0*ur*ur
        pL   = (gammaL-1.d0)*rhoL*(eL-qL) - gammaL*pinfL
        pR   = (gammaR-1.d0)*rhoR*(eR-qR) - gammaR*pinfR
        pstar= (pL+pR)/2.d0
        ustar= (rhoL*cl*ul+rhoR*cr*ur+pL-pR)/(rhoL*cl+rhoR*cr)
        fsol(1) = -ustar
        fsol(2) = 0.d0
        fsol(3) = pstar
        fsol(4) = pstar*ustar
      endif
      return
      END SUBROUTINE hllclag

!=======================================================================
!======================= init ==========================================
!=======================================================================
!
SUBROUTINE init
use global_data

   integer :: i, kf
   double precision :: alpha,rho,u,p,e,q
   double precision :: gamma
   double precision :: pinf
   double precision :: dl1, dl2, dl3, dl4
   double precision :: dsdx1, dsdx2, dsdx3, dsdx4
   double precision :: xent, xsort
   double precision :: area0, pi
!  double precision :: 
!  double precision :: prim(ieq,im,ifl)
!  double precision :: uc(ieq,im,ifl)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!          Initialization
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! compute (constant) cell size
!   dx = ltube/dfloat(im-1)
   dx = ltube/dfloat(im-2)

!  Thermodynamic data
   gammaf(1)= gamma1
   pinff(1) = pinf1
!
   gammaf(2)= gamma2
   pinff(2) = pinf2
!
   qf(1)    = q1
   qf(2)    = q2
!
! ordering of primitive variables array (prim0) for initial data:
!      1= alpha, 2= rho, 3= u , 4= p
!
!  INITIAL CONDITIONS
!
!   first index  = variable
!   second index = fluid

!    LEFT STATE (first segment)
!    fluid 1
       prim0l(1,1) = 0.5d0
       prim0l(2,1) = 1000.d0
       prim0l(3,1) = -0.d0
       prim0l(4,1) = 10.d5
!       prim0l(4,1) = 9.7d5
!    fluid 2 
       prim0l(1,2) = 0.5d0
!       prim0l(2,2) = 1.d0
       prim0l(2,2) = 0.54816d0
       prim0l(3,2) = -0.d0
       prim0l(4,2) = 10.d5
!       prim0l(4,2) = 9.7d5
!    RIGHT STATE (second segment)
!    fluid 1
       prim0r(1,1) = 0.5d0
       prim0r(2,1) = 1000.d0
       prim0r(3,1) = -0.d0
       prim0r(4,1) = 5.d5
!       prim0r(4,1) = 9.7d5
!    fluid 2
       prim0r(1,2) = 0.5d0
!       prim0r(2,2) = 1.d0
       prim0r(2,2) = 0.54816d0
       prim0r(3,2) = -0.d0
       prim0r(4,2) = 5.d5
!       prim0r(4,2) = 9.7d5
!
!   Conservative variable vectors
!
        do i=1,idisc2
!         cell center position
!          xx(i) = dx*(i-1)+dx/2.d0
          xx(i) = dx*(i-2)+dx/2.d0
          do kf=1,ifl
             alpha= prim0l(1,kf)
             rho  = prim0l(2,kf)
             u    = prim0l(3,kf)
             p    = prim0l(4,kf)
!
             gamma= gammaf(kf)
             pinf = pinff(kf)
             q    = qf(kf)
             call eospe(rho,p,gamma,pinf,q,e)
!
             prim(1,i,kf) = alpha
             prim(2,i,kf) = rho
             prim(3,i,kf) = u
             prim(4,i,kf) = p
!
             uc(1,i,kf) = alpha
             uc(2,i,kf) = alpha*rho
             uc(3,i,kf) = alpha*rho*u
             uc(4,i,kf) = alpha*rho*(e+u**2/2.d0)

          enddo
        enddo
!
        do i=idisc2+1,idisc4
!         cell center position
!          xx(i) = dx*(i-1)+dx/2.d0
          xx(i) = dx*(i-2)+dx/2.d0
          do kf=1,ifl
             alpha= prim0r(1,kf)
             rho  = prim0r(2,kf)
             u    = prim0r(3,kf)
             p    = prim0r(4,kf)
!
             gamma= gammaf(kf)
             pinf = pinff(kf)
             q    = qf(kf)
             call eospe(rho,p,gamma,pinf,q,e)
!
             prim(1,i,kf) = alpha
             prim(2,i,kf) = rho
             prim(3,i,kf) = u
             prim(4,i,kf) = p
!
             uc(1,i,kf) = alpha
             uc(2,i,kf) = alpha*rho
             uc(3,i,kf) = alpha*rho*u
             uc(4,i,kf) = alpha*rho*(e+u**2/2.d0)

          enddo
        enddo
!
!    Geometrical calculations 
!
! PIECEWISE LINEAR NOZZLE	
!    Segment lengths
       dl1 = xx(idisc1)+dx/2.d0
       dl2 = xx(idisc2)+dx/2.d0-dl1
       dl3 = xx(idisc3)+dx/2.d0-dl1-dl2
       dl4 = xx(idisc4)+dx/2.d0-dl1-dl2-dl3
!    Slope for each segment
!       dsdx1=(ssort1-sent1)/dl1
       dsdx1 = (ssort1-sent1)/(dl1-dx)
       dsdx2 = (ssort2-sent2)/dl2
       dsdx3 = (ssort3-sent3)/dl3
!       dsdx4 = (ssort4-sent4)/dl4
       dsdx4 = (ssort4-sent4)/(dl4-dx)
!
!     For each cell we compute the left cross section (surfent)
!     as well as the right one (surfsort)
        do i=1,idisc1   
          xent = xx(i)-dx/2.d0 ! left cell boundary position
          xsort= xx(i)+dx/2.d0! right cell boundary position
          surfent(i) = dsdx1*xent+sent1
          surfsort(i)= dsdx1*xsort+sent1
!      Linear relation for surfaces = rectangular cross sections
!      (constant depth)
        enddo
        do i=idisc1+1,idisc2
          xent = xx(i)-dx/2.d0-dl1
          xsort= xx(i)+dx/2.d0-dl1
          surfent(i) = dsdx2*xent+sent2
          surfsort(i)= dsdx2*xsort+sent2
        enddo
        do i=idisc2+1,idisc3
          xent = xx(i)-dx/2.d0-dl1-dl2
          xsort= xx(i)+dx/2.d0-dl1-dl2
          surfent(i) = dsdx3*xent+sent3
          surfsort(i)= dsdx3*xsort+sent3
        enddo
        do i=idisc3+1,idisc4
          xent = xx(i)-dx/2.d0-dl1-dl2-dl3
          xsort= xx(i)+dx/2.d0-dl1-dl2-dl3
          surfent(i) = dsdx4*xent+sent4
          surfsort(i)= dsdx4*xsort+sent4
        enddo
!      goto 100
! SMOOTH COSINE NOZZLE	
!... nozzle area multiplier coefficient, and set PI
        area0 = 1.d0
        PI    = 4.0D0*DATAN(1.0D0)

!... Grid Spacing, node positions, and cell-edge areas
      do i=1,im
        xent = xx(i)-dx/2.d0
        xsort= xx(i)+dx/2.d0
        surfent(i) = area0*(1.D0+0.5D0*DCOS(xent*2.d0*PI))
        surfsort(i)= area0*(1.D0+0.5D0*DCOS(xsort*2.d0*PI))
!        surfent(i) = area0*(1.D0+0.5D0*DCOS(xent*1.d0*PI))
!        surfsort(i)= area0*(1.D0+0.5D0*DCOS(xsort*1.d0*PI))
      enddo
  100 continue
      surfent(1)  = surfent(2)
      surfsort(1) = surfent(2)
      surfent(im) = surfsort(im-1)
      surfsort(im)= surfsort(im-1)
!
!     Cell volume computation
        do i=1,im
          vol(i) = 0.5d0*dx*(surfent(i)+surfsort(i))
        enddo
!
!     Cross section available for the fluids flow
        do i=1,im-1
          surffluide(i) = dmin1(surfsort(i),surfent(i+1))
        enddo

!!$       print *," "
!!$       print *,"In init"
!!$       do j = 1,im
!!$          do k = 1,ifl
!!$             do i = 1,ieq
!!$                write(6,'(i4,i4,i4,f17.6)')i,j,k,uc(i,j,k)
!!$             end do
!!$          end do
!!$       end do
!!$       print *," "

!
!... Initialize temporal derivatives
!     CALL Derv(ieq, tv, y, ydot)
!
return
END SUBROUTINE init

!=======================================================================
!======================= inlet_BC ======================================
!=======================================================================
!
SUBROUTINE inlet_BC(ieq,im,ifl,kfL,kfR,rho0,p0,t,prim, &
                    gammaL,gammaR,pinfL,pinfR,qL,qR, &
                    pr,zr,ur,ustar,pstar,f)
implicit none

!... This subroutine computes the value of ustar at the inlet
!    boundary for the DEM method.

!... Dummy arguments
      integer :: kfL,kfR                         ! input
      integer :: ieq, im, ifl                    ! input
      double precision :: rho0, p0               ! input
      double precision :: t                      ! input
      double precision :: prim(ieq,im,ifl)       ! input
      double precision :: gammaL, pinfL          ! input
      double precision :: gammaR, pinfR          ! input
      double precision :: qL, qR                 ! input
      double precision :: pr,zr,ur               ! output
      double precision :: ustar, pstar           ! output
      double precision :: f(ieq,im,ifl,ifl)      ! output
!     double precision :: 

!... Local variables
      integer :: i
      double precision :: gm
      double precision :: rhor, er, cr
      double precision :: dk, dh
      double precision :: ff, df
      double precision :: pib,pfb,trise,pi,pmean,pamp
      double precision :: rhostar, estar
      double precision :: alphar
      double precision :: e0, c0, z0


!  Tank inlet pressure

!      pib  = 2.0d5
!      pfb  = 10.d5
!      trise= 1.0d0
!      pi   = 4.d0*datan(1.d0)
!      pmean= (pib+pfb)/2.d0
!      pamp = pmean-pib

!     if(t.le.trise) p0=pmean+pamp*dsin(pi*t/trise - pi/2.d0)
!      p0 = pmean+pamp*dsin(pi*t/trise - pi/2.d0)
!      p0 = 2.d5+8.d5*(1.d0-dexp(-t/trise))
!      p0 = 10.d5

      alphar= prim(1,2,kfR)
      rhor  = prim(2,2,kfR)
      ur    = prim(3,2,kfR)
      pr    = prim(4,2,kfR)
      call eospe(rhor,pr,gammaR,pinfR,qR,er)
      call eospression(rhor,er,gammaR,pinfR,qR,pr,cr)
      zr    = rhor*cr
!      zr    = rhor*(cr-ur)
!
      call eospe(rho0,p0,gammaL,pinfL,qL,e0)
      call eospression(rho0,e0,gammaL,pinfL,qL,p0,c0)
      z0    = rho0*c0

      gm    = gammaL-1.d0
      dk    = (p0+pinfL)/rho0**gammaL
      dh    = gammaL*(p0+pinfL)/rho0/gm
!
!  initial guess of ustar for newton iteration
      ustar = (p0-pr+zr*ur)/(zr+z0)
!      ustar = 10.d0
!      ustar=300.d0
!      if(t.gt.0.d0) ustar=ur
!
      do i=1,100
!         print*, 'gammaL= ', gammaL, 'pinfL= ', pinfL
!	  print*, 'p0= ', p0, 'rho0= ', rho0
!         print*, 'dh= ', dh
!         print*, 'ustar= ', ustar, 'cr= ', cr
!	  print*, 'kfL= ', kfL, 'kfR= ', kfR
!         if((dh-0.5d0*ustar**2).lt.0.d0) print*, 'yes'
!
         ff = dk*((dh-0.5d0*ustar**2)/dk*gm/gammaL)**(gammaL/gm) &
             -pinfL-pr+zr*(ur-ustar)
         df =-(dk*(gm/gammaL/dk)**(gammaL/gm)*gammaL/gm*ustar* &
                (dh-0.5d0*ustar**2)**(1.d0/gm) + zr)
         if (dabs(ff) .LT. 1.d-6) goto 25
         ustar = ustar - ff/df
!
         if (ustar.Lt.0.d0) then
            write(6,'(a)')'Inlet fluid: convergence failure'
!            print*, 'kfL= ', kfL, 'kfR= ', kfR
!            print*, 'ustar= ', ustar, 'ff= ', ff
!            print*, 'pr= ', pr, 'ur= ', ur, 'rhor= ', rhor
!	     print*, 'dk= ', dk, 'dh= ', dh
         endif
!
      enddo
25    CONTINUE

      pstar   = pr-zr*(ur-ustar)
      rhostar = rho0*((pstar+pinfL)/(p0+pinfL))**(1.d0/gammaL)
      prim(2,1,kfr) = 2.d0*rhostar - prim(2,2,kfr)
      prim(3,1,kfr) = 2.d0*ustar - prim(3,2,kfr)
      prim(4,1,kfr) = 2.d0*pstar - prim(4,2,kfr)
      f(1,1,kfL,kfR) = 0.d0
      f(2,1,kfL,kfR) = rhostar*ustar
      f(3,1,kfL,kfR) = rhostar*ustar**2+pstar
      call eospe(rhostar,pstar,gammaL,pinfL,qL,estar)
      f(4,1,kfL,kfR) = rhostar*ustar*(estar+ustar**2/2.d0) &
                      +pstar*ustar

return
END SUBROUTINE inlet_BC

!=======================================================================
!======================= outlet_BC =====================================
!=======================================================================
!
SUBROUTINE outlet_BC(ieq,im,ifl,kf1,kf2,prim,gamma1,pinf1,q1,pstar, &
                            rhostar,ustar,el,f)
implicit none

!... This subroutine computes the value of ustar at the outlet
!    boundary for the DEM method.

!... Dummy arguments
      integer :: kf1,kf2                          ! input
      integer :: ieq, im, ifl                     ! input
      double precision :: prim(ieq,im,ifl)        ! input
      double precision :: gamma1, pinf1, q1       ! input
!     double precision :: pinf2                   ! input
      double precision :: pstar                   ! input
      double precision :: rhostar                 ! output  
      double precision :: ustar                   ! output
      double precision :: el                      ! output
      double precision :: f(ieq,im,ifl,ifl)       ! output

!... Local variables
      integer :: j
      double precision :: rhol, pl, ul, zl, cl, dml
      double precision :: estar


      j = im-1
!     alpha1l=prim(1,j,kf1)
      rhol= prim(2,j,kf1)
      ul  = prim(3,j,kf1)
      pl  = prim(4,j,kf1)
      call eospe(rhol,pl,gamma1,pinf1,q1,el)
      call eospression(rhol,el,gamma1,pinf1,q1,pl,cl)
      zl  = rhol*cl
      dml = dabs(ul)/cl

      if(dml.lt.1.d0)then    !subsonic
!         pstar  = 1.d5
!         pstar  = 9.7d5
         rhostar= rhol+(pstar-pl)/cl**2
         ustar  = (pl-pstar)/zl+ul
         call eospe(rhostar,pstar,gamma1,pinf1,q1,estar)
         prim(2,im,kf1) = 2.d0*rhostar - prim(2,im-1,kf1)
         prim(3,im,kf1) = 2.d0*ustar - prim(3,im-1,kf1)
         prim(4,im,kf1) = 2.d0*pstar - prim(4,im-1,kf1)
         f(1,j,kf1,kf2) = 0.d0
         f(2,j,kf1,kf2) = rhostar*ustar
         f(3,j,kf1,kf2) = rhostar*ustar**2+pstar
         f(4,j,kf1,kf2) = rhostar*ustar*(estar+ustar**2/2.d0) &
                         +pstar*ustar
         el = estar
      else                   !supersonic
         f(1,j,kf1,kf2) = 0.d0
         f(2,j,kf1,kf2) = rhol*ul
         f(3,j,kf1,kf2) = rhol*ul**2+pl
!        el=(pl+gamma1*pinf1)/rhol/(gamma1-1.d0)
!        call eospe(rhol,pl,gamma1,pinf1,q1,el)
         f(4,j,kf1,kf2) = rhol*ul*(el+ul**2/2.d0) &
                          +pl*ul
      endif

return
END SUBROUTINE outlet_BC

!=======================================================================
!======================= close_LB ======================================
!=======================================================================
!
SUBROUTINE closed_LB(ieq,im,ifl,kfL,kfR,t,prim,gamma,pinf,q,pstar,f)
implicit real*8(a-h,l-z)
!
!... This subroutine computes the left edge fluxes for a
!    closed left boundary
!
      double precision :: prim(ieq,im,ifl)
      double precision :: f(ieq,im,ifl,ifl)
      j=2
      prim(1,1,kfR)= prim(1,j,kfR)
      rho          = prim(2,j,kfR)
      u            = prim(3,j,kfR)
      p            = prim(4,j,kfR)
      call eospe(rho,p,gamma,pinf,q,e)
      call eospression(rho,e,gamma,pinf,q,p,c)
      pstar = p - rho*c*u
      f(1,1,kfL,kfR) = 0.d0
      f(2,1,kfL,kfR) = 0.d0
      f(3,1,kfL,kfR) = pstar
      f(4,1,kfL,kfR) = 0.d0
      
return
END SUBROUTINE closed_LB

!=======================================================================
!======================= close_RB ======================================
!=======================================================================
!
SUBROUTINE closed_RB(ieq,im,ifl,kfL,kfR,t,prim,gamma,pinf,q,pstar,f)
implicit real*8(a-h,l-z)
!
!... This subroutine computes the right edge fluxes for a
!    closed right boundary
!
      double precision :: prim(ieq,im,ifl)
      double precision :: f(ieq,im,ifl,ifl)
      j=im-1
      prim(1,im,kfL)= prim(1,j,kfL)
      rho           = prim(2,j,kfL)
      u             = prim(3,j,kfL)
      p             = prim(4,j,kfL)
      call eospe(rho,p,gamma,pinf,q,e)
      call eospression(rho,e,gamma,pinf,q,p,c)
      pstar = p + rho*c*u
      f(1,j,kfL,kfR) = 0.d0
      f(2,j,kfL,kfR) = 0.d0
      f(3,j,kfL,kfR) = pstar
      f(4,j,kfL,kfR) = 0.d0
      
return
END SUBROUTINE closed_RB

!=======================================================================
!======================= Print_Results =================================
!=======================================================================
!
SUBROUTINE Print_Results(k)
use global_data
implicit real*8(a-h,l-z)

!... Print results contained in uc and prim

!   double precision :: prim(ieq,im,ifl)
!   double precision :: uc(ieq,im,ifl)

   integer :: kf, k, i
   double precision :: alpha, rho, u, e, p, c, t
   double precision :: gamma, cv, pinf

!  Open the files for the calculated data output:
   open (unit=11,file='res1.dat',position='append')
   open (unit=12,file='res2.dat',position='append')
   open (unit=21,file='mel.dat',position='append')


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!   Print variables (to Tecplot x-y plot format files)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!  Fluid variables
!
      do kf=1,ifl
         write(10+kf,'(a,f10.5,a)') 'TITLE = "T= ',temps,' "'
!         write(10+kf,'(a,a)') 'VARIABLES = "XX","alpha","density","velocity' &
!               ,'","pressure","Mach","Temperature"'
         write(10+kf,'(a,a)') 'VARIABLES = "XX","alpha","density","velocity' &
               ,'","pressure","m-flow","Temperature"'
!        write(10+kf,'(a36,a31)') 'VARIABLES = "XX","density","velocity' &
!               ,'","pressure","vol","surffluide"'
         write(10+kf,'(a,i4,a)') 'ZONE T= "', k/kprint ,'"  I=  100'
      do i=2,im-1
              alpha= uc(1,i,kf)
              rho  = uc(2,i,kf)/alpha
              u    = uc(3,i,kf)/uc(2,i,kf)
              e    = uc(4,i,kf)/uc(2,i,kf)-u**2/2.d0
              gamma= gammaf(kf)
              pinf = pinff(kf)
              q    = qf(kf)
               call eospression(rho,e,gamma,pinf,q,p,c)
               if(kf.eq.1) then
                  cv = cv1
                else
                  cv = cv2
               endif
              t    = (p+pinf)/(cv*(gamma-1)*rho)
              area = 0.5d0*(surffluide(i-1)+surffluide(i))
!       write(10+kf,140) xx(i),prim(1,i,kf),prim(2,i,kf),prim(3,i,kf), &
!             prim(4,i,kf)/1.d6,u/c,t
       write(10+kf,140) xx(i),prim(1,i,kf),prim(2,i,kf),prim(3,i,kf), &
             prim(4,i,kf)/1.d6,uc(3,i,kf)*area,t
!       write(10+kf,140) xx(i),prim(1,i,kf),prim(2,i,kf),prim(3,i,kf), &
!             prim(4,i,kf)/1.d6,contact(i,kf,kf)*fkeuler(2,i,kf,kf)*surffluide(i),t
      enddo
       write(10+kf,*)
      enddo
!
!  Mixture variables
!

      WRITE(21,'(a,f11.8,a)') 'TITLE = "T= ',temps,' "'
!      WRITE(21,'(a36,a31)') 'VARIABLES = "XX","density","velocity', &
!           '","pressure","vol","surffluide"'
      WRITE(21,'(a36,a31)') 'VARIABLES = "XX","density","velocity', &
           '","pressure","vol","m-flow"'
      WRITE(21,'(a,i4,a)') 'ZONE T= "', k/kprint ,'"  I=  100'
      DO i=2,im-1
         romel(i) = 0.d0
         pmel(i)  = 0.d0
         umel(i)  = 0.d0
      DO kf=1,ifl
           alpha= uc(1,i,kf)
           rho  = uc(2,i,kf)/alpha
           u    = uc(3,i,kf)/uc(2,i,kf)
           e    = uc(4,i,kf)/uc(2,i,kf)-u**2/2.d0
           gamma= gammaf(kf)
           pinf = pinff(kf)
           q    = qf(kf)
           CALL eospression(rho,e,gamma,pinf,q,p,c)
           romel(i) = romel(i)+uc(2,i,kf)
           pmel(i)  = pmel(i)+alpha*p
           umel(i)  = umel(i)+uc(3,i,kf)
       END DO
         area   = 0.5d0*(surffluide(i-1)+surffluide(i))
         umel(i)= umel(i)/romel(i)
!         WRITE(21,140) xx(i),romel(i),umel(i),pmel(i),vol(i), &
!                       surffluide(i)
         WRITE(21,150) xx(i),romel(i),umel(i),pmel(i),vol(i), &
                       romel(i)*umel(i)*area
       END DO
       WRITE(21,*)
140    FORMAT(7(e16.7,1x))
150    FORMAT(6(e16.7,1x))

       CLOSE(11)
       CLOSE(12)
       CLOSE(21)

return
END SUBROUTINE Print_Results

!=======================================================================
!======================= END ===========================================
!=======================================================================

