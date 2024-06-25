!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

SUBROUTINE tsfun1DLIN(x,ph,phd)

!*******************************************************************************
!*** DESCRIPTION: LINEAR BASIS FUNCTIONS EVALUATION (1D DOMAIN)
!*******************************************************************************

IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION,INTENT(OUT)::ph(2)
DOUBLE PRECISION,INTENT(OUT)::phd(2)

!*** BASIS FUNCTIONS  
ph(1)=1.d0-x
ph(2)=x
phd(1)=-1.d0
phd(2)=1.d0
     
END SUBROUTINE

!*******************************************************************************

SUBROUTINE  tsfun2DLIN(x,y,phi,phic,phie)

!*******************************************************************************
!*** DESCRIPTION: LINEAR BASIS FUNCTIONS EVALUATION FOR TRIAGULAR ELEMENTS
!*** (2D DOMAIN)
!*******************************************************************************

IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION,INTENT(IN)::y
DOUBLE PRECISION,INTENT(OUT)::phi(6)
DOUBLE PRECISION,INTENT(OUT)::phic(6)
DOUBLE PRECISION,INTENT(OUT)::phie(6)

!*** BASIS FUNCTIONS
phi(1)=1.d0-x-y
phi(2)=x
phi(3)=y
!*** BASIS FUNCTIONS DERIVATIVES
phic(1)=-1.d0
phic(2)=1.d0
phic(3)=0.d0
phie(1)=-1.d0
phie(2)=0.d0
phie(3)=1.d0

END SUBROUTINE  

!*******************************************************************************

SUBROUTINE  tsfun2DQUAD(x,y,phi,phic,phie)

!*******************************************************************************
!*** DESCRIPTION: QUADRATIC BASIS FUNCTIONS EVALUATION FOR TRIAGULAR ELEMENTS
!*** (2D DOMAIN)
!*******************************************************************************

IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION,INTENT(IN)::y
DOUBLE PRECISION,INTENT(OUT)::phi(6)
DOUBLE PRECISION,INTENT(OUT)::phic(6)
DOUBLE PRECISION,INTENT(OUT)::phie(6)

!*** BASIS FUNCTIONS 
phi(1)=(1.d0-x-y)*(1.d0-2.d0*x-2.d0*y)
phi(2)=x*(2.d0*x-1.d0)
phi(3)=y*(2.d0*y-1.d0)
phi(4)=4.d0*x*(1.d0-x-y)
phi(5)=4.d0*x*y
phi(6)=4.d0*y*(1.d0-x-y)
!*** BASIS FUNCTIONS DERIVATIVES
phic(1)=4.d0*x+4.d0*y-3.d0
phic(2)=4.d0*x-1.d0
phic(3)=0.d0
phic(4)=-8.d0*x-4.d0*y+4.d0
phic(5)=4.d0*y
phic(6)=-4.d0*y
phie(1)=4.d0*x+4.d0*y-3.d0
phie(2)=0.d0
phie(3)=4.d0*y-1.d0
phie(4)=-4.d0*x
phie(5)=4.d0*x
phie(6)=-4.d0*x-8.d0*y+4.d0

END SUBROUTINE 

!*******************************************************************************

SUBROUTINE tsfun2DQUADSTR(x,y,phi,phic,phie)

!*******************************************************************************
!*** DESCRIPTION: QUADRATIC BASIS FUNCTIONS EVALUATION FOR RECTANGULAR ELEMENTS
!*** (2D DOMAIN)
!*******************************************************************************

IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION,INTENT(IN)::y
DOUBLE PRECISION,INTENT(OUT)::phi(9)
DOUBLE PRECISION,INTENT(OUT)::phic(9)
DOUBLE PRECISION,INTENT(OUT)::phie(9)
DOUBLE PRECISION::l1
DOUBLE PRECISION::l2
DOUBLE PRECISION::l3
DOUBLE PRECISION::dl1
DOUBLE PRECISION::dl2
DOUBLE PRECISION::dl3
DOUBLE PRECISION::c

l1(c)=2.d0*c**2-3.d0*c+1.d0
l2(c)=-4.d0*c**2+4.d0*c
l3(c)=2.d0*c**2-c
dl1(c)=4.d0*c-3.d0
dl2(c)=-8.d0*c+4.d0
dl3(c)=4.d0*c-1.d0
!*** BASIS FUNCTIONS 
phi(1)=l1(x)*l1(y)
phi(2)=l1(x)*l2(y)
phi(3)=l1(x)*l3(y)
phi(4)=l2(x)*l1(y)
phi(5)=l2(x)*l2(y)
phi(6)=l2(x)*l3(y)
phi(7)=l3(x)*l1(y)
phi(8)=l3(x)*l2(y)
phi(9)=l3(x)*l3(y)
!*** BASIS FUNCTIONS DERIVATIVES
phic(1)=l1(y)*dl1(x)
phic(2)=l2(y)*dl1(x)
phic(3)=l3(y)*dl1(x)
phic(4)=l1(y)*dl2(x)
phic(5)=l2(y)*dl2(x)
phic(6)=l3(y)*dl2(x)
phic(7)=l1(y)*dl3(x)
phic(8)=l2(y)*dl3(x)
phic(9)=l3(y)*dl3(x)
phie(1)=l1(x)*dl1(y)
phie(2)=l1(x)*dl2(y)
phie(3)=l1(x)*dl3(y)
phie(4)=l2(x)*dl1(y)
phie(5)=l2(x)*dl2(y)
phie(6)=l2(x)*dl3(y)
phie(7)=l3(x)*dl1(y)
phie(8)=l3(x)*dl2(y)
phie(9)=l3(x)*dl3(y)

END SUBROUTINE 

!*******************************************************************************

SUBROUTINE gauPointsWeights(eltype,nelnodes,gp,wgp,Numgps)

!*******************************************************************************
!*** DESCRIPTION: GAUSS POINTS AND WEIGHTS EVALUATION
!*******************************************************************************

IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::eltype
INTEGER(KIND=8),INTENT(OUT)::nelnodes
INTEGER(KIND=8),INTENT(OUT)::Numgps
DOUBLE PRECISION,INTENT(OUT)::gp(3)
DOUBLE PRECISION,INTENT(OUT)::wgp(3)

!*** ELEMENT TYPE 
!*** (1=LINEAR TRIANGULAR, 2=QUADRATIC TRIANGULAR, 3=QUADRATIC RECTANGULAR)
SELECT CASE (eltype)
  CASE(1)
    nelnodes=3 !*** NUMBER OF NODES ON ONE ELEMENT
    gp(1)=1.d0/3.d0 !*** GAUSS POINTS AND WEIGHTS
    wgp(1)=1.d0
    Numgps=1 !*** NUMBER OF GAUSS POINTS
  CASE(2)
    nelnodes=6
    gp(1)=2.d0/3.d0
    gp(2)=1.d0/6.d0
    gp(3)=1.d0/6.d0
    wgp(1)=1.d0/3.d0
    wgp(2)=1.d0/3.d0
    wgp(3)=1.d0/3.d0
    Numgps=3 !*** NUMBER OF GAUSS POINTS
  CASE(3)
    nelnodes=9
    gp(1)=(1.d0-SQRT(3.d0/5.d0))/2.d0
    gp(2)=0.5d0
    gp(3)=(1.d0+SQRT(3.d0/5.d0))/2.d0
    wgp(1)=5.d0/18.d0
    wgp(2)=8.d0/18.d0
    wgp(3)=5.d0/18.d0
    Numgps=3
END SELECT

END SUBROUTINE 

!*******************************************************************************
