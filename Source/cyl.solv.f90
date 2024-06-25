!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

SUBROUTINE CYLSolver

!*******************************************************************************
!*** DESCRIPTION: CONVENTIONAL YOUNG-LAPLACE EQUATION SOLVER
!*******************************************************************************

USE CommonVars
USE initsol
USE nodnum
USE defpar
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::GENcnt
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::AMcntyl
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::am2yl
DOUBLE PRECISION::normsum

!*** SOLVER INITIALIZATION
CALL CHDIR(ADJUSTL(Dir))
CALL defparamCYL !*** READ CYL PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))

!*** ALLOCATE MEMORY
ALLOCATE (uYL(noukYL),uoYL(noukYL),am1(50*noukYL),am2(50*noukYL),xvec(noukYL), &
res(noukYL),thptYL(nnodesYL),dthptYL(nnodesYL),nopYL(nellYL,3))

!*** GENERATE (OR READ) INITIAL SOLUTION FOR THE CONVENTIONAL YOUNG-LAPLACE
CALL mesh1DLIN(0.d0,Pi/2,nellYL,thptYL,dthptYL) !*** EVENLY DISTRIBUTED MESH
DO i=1,nnodesYL
  uoYL(i)=1.25
ENDDO
uoYL(nnodesYL+1)=1.d0

!*** NODAL NUMBERING
CALL nodnumbCYL

!*** COUNTERS INITIALIZATION
GENcnt=1
NRcnt=1
inr=0

!*** MATRICES INITIALIZATION
DO i=1,noukYL
  xvec(i)=0.d0 !*** INITIALIZE X VECTOR
  uYL(i)=uoYL(i) !*** INITIALIZE SOLUTION
ENDDO
AMcnt=0
DO i=1,SIZE(am2)
  am2(i)=0.d0 !*** INITIALIZE JACOBIAN MATRIX PATTERN
ENDDO

SELECT CASE (GENcnt)
  CASE(1)
!*** FIND SPARSE JACOBIAN MATRIX PATTERN FOR THE CONVENTIONAL YOUNG-LAPLACE
    CALL neigElem(noukYL,nellYL,nelnodesYL,nopYL,kindYL,neighYL,noukYL)
    CALL JacobPtrn(noukYL,nelnodesYL,nopYL,kindYL,neighYL)
!*** SAVE JACOBIAN MATRIX PATTERN
    ALLOCATE(am2yl(AMcnt))
    AMcntyl=AMcnt
    DO i=1,AMcnt
      am2yl(i)=am2(i)
    ENDDO
  CASE DEFAULT
!*** LOAD JACOBIAN MATRIX PATTERN
    AMcnt=AMcntyl
    ALLOCATE(am1(AMcnt),am2(AMcnt),xvec(noukYL),res(noukYL))
    DO i=1,AMcnt
      am2(i)=am2yl(i)
    ENDDO
    DO i=1,noukYL
      xvec(i)=0
      uYL(i)=uoYL(i)
    ENDDO
END SELECT

!*** SOLVE AUGMENTED YOUNG-LAPLACE EQUATION
DO WHILE (inr==0)!*** START NEWTON-RAPHSON ITERATIONS
  DO i=1,noukYL
    uYL(i)=uYL(i)+xvec(i)
  ENDDO
  NRcnt=NRcnt+1
  CALL axbCYL !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING
  CALL ErrorNormCalc(noukYL,xvec,normsum) !*** CALCULATE ERROR NORM
  IF (normsum<=PPres) THEN
    inr=1
  ENDIF
  WRITE(*,'(a,E10.3)')'AYL ERROR NORM:',normsum
  WRITE(*,'(a,i4,/)')'NEWTON ITERATIONS:',NRcnt-1
  IF (normsum>1.d3) THEN !*** CHECK FOR DIVERGENCE
    WRITE(*,'(a)')'CYL ERROR NORM>1.e3 - KILLING PROCESS...'
    STOP
  ENDIF
  IF (normsum<=PPres) THEN
    inr=1
  ENDIF
ENDDO

!*** PRINT RESULTS
CALL printres(GENcnt)

!*** CLEAR MEMORY
DEALLOCATE(am1,am2,xvec,res,uYL,uoYL,)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE axbCYL

!*******************************************************************************
!*** DESCRIPTION: BOUNDARY CONDITION IMPLEMENTATION AND SYSTEM SOLVING 
!*** FOR THE COVENTIONAL YOUNG-LAPLACE EQUATION
!*******************************************************************************

USE CommonVars
USE qsort_module
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::k
INTEGER(KIND=8)::icnt1
INTEGER(KIND=8)::icnt2
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::am2temp
DOUBLE PRECISION::jacob
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::am1temp

!*** MATRICES INITIALIZATION
DO i=1,noukYL
  res(i)=0.d0
ENDDO
DO i=1,SIZE(am1)
  am1(i)=0.d0
ENDDO

!*** ASSEMBLY RESIDUALS (RES) AND JACOBIAN (AM) MATRICES
DO nell=1,nellYL
  CALL abfindCYL
ENDDO

!*** VOLUME CONSERVATION
res(nnodesYL+1)=res(nnodesYL+1)+2.d0

!*** CONTACT ANGLE BOUNDARY CONDITION
res(nnodesYL)=res(nnodesYL)+uYL(nnodesYL)*COS(thYoungYL)
icnt1=(nnodesYL-1)*noukYL+nnodesYL
jacob=-COS(thYoungYL)
CALL buildjacob_nor(icnt1,jacob)

!*** SORT JACOBIAN MATRIX
IF (NRcnt==1) THEN
  ALLOCATE(am1temp(AMcnt),am2temp(AMcnt))
  DO i=1,AMcnt
    am1temp(i)=am1(i)
    am2temp(i)=am2(i)
  ENDDO
  DEALLOCATE(am1)
  DEALLOCATE(am2)
  ALLOCATE(am1(AMcnt))
  ALLOCATE(am2(AMcnt))
  DO i=1,AMcnt
    am1(i)=am1temp(i)
    am2(i)=am2temp(i)
  ENDDO
  DEALLOCATE(am1temp)
  DEALLOCATE(am2temp)
!  CALL qsort(am2,am1) !*** THERE IS NO NEED TO SORT HERE
ENDIF

!*** SOLVE THE LINEAR SYSTEM USING MUMPS DIRECT SOLVER
noukLS=noukYL
CALL linsolver

END SUBROUTINE

!*******************************************************************************

SUBROUTINE abfindCYL

!*******************************************************************************
!*** DESCRIPTION: RESIDUALS & JACOBIAN MATRICES ASSEMBLY FOR THE
!*** COVENTIONAL YOUNG-LAPLACE EQUATION
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::n
INTEGER(KIND=8)::m
INTEGER(KIND=8)::m1
INTEGER(KIND=8)::n1
INTEGER(KIND=8)::ngl(2)
INTEGER(KIND=8)::icnt1
DOUBLE PRECISION::gpAYL
DOUBLE PRECISION::wgpAYL
DOUBLE PRECISION::r
DOUBLE PRECISION::rd
DOUBLE PRECISION::th
DOUBLE PRECISION::th1
DOUBLE PRECISION::ph(2)
DOUBLE PRECISION::phd(2)
DOUBLE PRECISION::phth(2)
DOUBLE PRECISION::jacob

DO i=1,2
  ngl(i)=nopYL(nell,i)
ENDDO

!*** GAUSS POINTS AND WEIGHTS
gpAYL=0.5d0
wgpAYL=1.d0

CALL tsfun1DLIN(gpAYL,ph,phd) !*** CALCULATE BASIS FUNCTIONS

!*** INITIALIZE VARIABLES
r=0.d0; rd=0.d0; th=0.d0; th1=0.d0;

!*** ISOPARAMETRIC MAPPING
DO n=1,2
  th=th+thptYL(ngl(n))*ph(n) !*** TH=POLAR COORDINATE AT THE GAUSS POINT
  th1=th1+thptYL(ngl(n))*phd(n) !*** TH1=THE C DERIVATIVE OF TH
  r=r+uYL(ngl(n))*ph(n) !*** R=RADIAL COORDINATE AT THE GAUSS POINT
ENDDO

DO n=1,2
  phth(n)=phd(n)/th1 !*** PHTH=TH DERIVATIVE OF THE BASIS FUNCTION
ENDDO

DO n=1,2
  rd=rd+uYL(ngl(n))*phth(n) !*** RD=TH DERIVATIVE OF R
ENDDO

!*** ASSEMBLY RESIDUALS

!*** VOLUME CONSERVATION TERM
res(nnodesYL+1)=res(nnodesYL+1)-wgpAYL*th1*(r**3)*SIN(th)

DO m=1,2
  m1=ngl(m)

!*** GRAVITY TERM
  res(m1)=res(m1)-wgpAYL*th1*ph(m)*((gYL*denYL*(roYL**2))/stYL)                &
  *(r**3)*COS(th)*SIN(th)

!*** CURVATURE TERMS
  res(m1)=res(m1)-wgpAYL*th1*(((2.d0*(r**2)*((SIN(th))**2)+(rd**2)*            &
  ((SIN(th))**2))*ph(m))/(((SIN(th)**2))*(r**2+rd**2))**0.5)
  res(m1)=res(m1)-wgpAYL*th1*((r*rd*(SIN(th)**2)*phth(m))/(((SIN(th)**2))      &
  *(r**2+rd**2))**0.5)

!*** PRESSURE TERMS
  res(m1)= res(m1)+wgpAYL*th1*ph(m)*uYL(nnodesYL+1)*(r**2)*SIN(th)

!*** D(RES)/D(K) TERMS
  icnt1=(m1-1)*noukYL+(nnodesYL+1)
  jacob=-wgpAYL*th1*ph(m)*(r**2)*SIN(th)
  CALL buildjacob_nor(icnt1,jacob)

!*** D(VOLUME CONSERVATION TERM)/D(R) TERMS
  icnt1=((nnodesYL+1)-1)*noukYL+m1
  jacob=wgpAYL*th1*3.d0*(r**2)*ph(m)*SIN(th)
  CALL buildjacob_nor(icnt1,jacob)

  DO n=1,2
    n1=ngl(n)

!*** RES=CURVATURE TERMS
    icnt1=(m1-1)*noukYL+n1
    jacob=+(-(2.d0*r**2*SIN(th)**2+rd**2*SIN(th)**2)*ph(m)*r*wgpAYL*th1*       &
    SIN(th)**2/((r**2+rd**2)*SIN(th)**2)**1.5)*ph(n)+(4.d0*ph(m)*r*wgpAYL*th1* &
    SIN(th)**2/((r**2+rd**2)*SIN(th)**2)**0.5)*ph(n)+(-(2.d0*r**2*SIN(th)**2+  &
    rd**2*SIN(th)**2)*ph(m)*rd*wgpAYL*th1*SIN(th)**2/((r**2+rd**2)*SIN(th)**2) &
    **1.5)*phth(n)+(2.d0*ph(m)*rd*wgpAYL*th1*SIN(th)**2/((r**2+rd**2)*SIN(th)  &
    **2)**0.5)*phth(n)+(-phth(m)*r**2*rd*wgpAYL*th1*SIN(th)**4/((r**2+rd**2)*  &
    SIN(th)**2)**1.5)*ph(n)+(phth(m)*rd*wgpAYL*th1*SIN(th)**2/((r**2+rd**2)*   &
    SIN(th)**2)**0.5)*ph(n)+(-phth(m)*r*rd**2*wgpAYL*th1*SIN(th)**4/((r**2+    &
    rd**2)*SIN(th)**2)**1.5)*phth(n)+(phth(m)*r*wgpAYL*th1*SIN(th)**2/((r**2+  &
    rd**2)*SIN(th)**2)**0.5)*phth(n)

!*** RES=PRESSURE TERMS
    jacob=jacob-wgpAYL*th1*ph(m)*uYL(nnodesYL+1)*2.d0*r*ph(n)*SIN(th)

!*** RES=GRAVITY TERMS
    jacob=jacob+wgpAYL*th1*ph(m)*((gYL*denYL*(roYL**2))/stYL)*3.d0*(r**2)*     &
    ph(n)*COS(th)*SIN(th) 
    
    CALL buildjacob_nor(icnt1,jacob)

  ENDDO

ENDDO

END SUBROUTINE

!*******************************************************************************
