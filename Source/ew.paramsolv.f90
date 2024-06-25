!********************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!********************************************************************************

SUBROUTINE EWParamSolver

!********************************************************************************
!*** DESCRIPTION: SOLVE FORCE BALANCE (AUGMENTED YOUNG-LAPLACE) AND 
!*** ELECTRIC FIELD DISTRIBUTION EQUATIONS
!*** PSEUDO ARC-LENGTH PARAMETRIC CONTINUATION
!********************************************************************************

USE CommonVars
USE initsol
USE nodnum
USE defpar
IMPLICIT NONE
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::AMcntyl
INTEGER(KIND=8)::dcouplflg
INTEGER(KIND=8)::GENcnt
INTEGER(KIND=8)::IntCallCnt
INTEGER(KIND=8)::i
INTEGER(KIND=8)::icnt1
INTEGER(KIND=8)::alloc
INTEGER(KIND=8)::solnum
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::am2yl
DOUBLE PRECISION::normsum
DOUBLE PRECISION::dcp_normsum
DOUBLE PRECISION::MinDropletHeightInit
DOUBLE PRECISION::sdiff
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::urYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uthYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uorYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uothYL
CHARACTER(LEN=1)::ElemOrder

!*** SOLVER INITIALIZATION
dsc=0.01d0 !*** DS STEP FOR PREUDO ARC-LENGTH CONTINUATION
OPEN(22,FILE='ResCon.txt') !*** OPEN CONTINUATION RESULTS FILE
WRITE(22,'(a)')'DropS RESULTS'
WRITE(22,'(a)')'VOLTAGE, MAXIMUM DROPLET HEIGHT'
CALL CHDIR(ADJUSTL(Dir))
CALL defparamYL !*** READ AYL PARAMETERS
nelnodesYL=7 !*** NUMBER OF UNKNOWNS FOR EACH ELEMENT
noukYL=noukYL+1 !*** AUGMENT THE UNKNOWNS WITH THE CONTINUATION PARAMETER
noukuoYL=noukuoYL+1
CALL defparamEI !*** READ EI PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))

!*** ALLOCATE MEMORY
ALLOCATE (uYL(noukYL),uoYL(noukYL),sptYL(nnodesYL),dsptYL(nnodesYL),xvec(noukYL),nopYL(nellYL,7),           &
res(noukYL),am1(50*noukYL),am2(50*noukYL),xptYLGP(nellYL),yptYLGP(nellYL),efieldYLGP(nellYL),               &
sptuoYL(nnodesYL),urYL(nnodesYL),uthYL(nnodesYL),uorYL(nnodesYL),uothYL(nnodesYL),sptYLGP(nellYL),          &
YLnodes(nnodesYL),uYLFirst(noukYL),uYLSecond(noukYL),dUdS(noukYL))

 !*** NODAL NUMBERING
CALL nodnumbAYL

!*** READ AUGMENTED YOUNG-LAPLACE FIRST INITIAL SOLUTION
solnum=10
OPEN(solnum,FILE='data_in_uminus1.txt',ACTION='READ')
CALL Sol1DRead(solnum,noukuoYL,nnodesuoYL,uoYL,sptuoYL)
uYLFirst(2*nnodesYL+1)=uoYL(2*nnodesuoYL+1)
uYLFirst(2*nnodesYL+2)=uoYL(2*nnodesuoYL+2)
uYLFirst(2*nnodesYL+3)=uoYL(2*nnodesuoYL+3)

!*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE FIRST INITIAL SOLUTION
CALL mesh1DLIN(0.d0,uYLFirst(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** EVENLY DISTRIBUTED MESH
DO i=1,nnodesuoYL
  uorYL(i)=uoYL(i)
ENDDO
icnt1=1
DO i=nnodesuoYL+1,2*nnodesuoYL
  uothYL(icnt1)=uoYL(i)
  icnt1=icnt1+1
ENDDO
CALL interpsol1D(uorYL,sptuoYL,nnodesuoYL,urYL,sptYL,nnodesYL)
CALL interpsol1D(uothYL,sptuoYL,nnodesuoYL,uthYL,sptYL,nnodesYL)
DO i=1,nnodesYL
  uYLFirst(i)=urYL(i)
ENDDO
icnt1=1
DO i=nnodesYL+1,2*nnodesYL
  uYLFirst(i)=uthYL(icnt1)
  icnt1=icnt1+1
ENDDO
uYLFirst(nnodesYL)=uoYL(nnodesuoYL)
uYLFirst(2*nnodesYL)=uoYL(2*nnodesuoYL)

!*** READ AUGMENTED YOUNG-LAPLACE SECOND INITIAL SOLUTION
solnum=10
OPEN(solnum,FILE='data_in_uo.txt',ACTION='READ')
CALL Sol1DRead(solnum,noukuoYL,nnodesuoYL,uoYL,sptuoYL)
uYLSecond(2*nnodesYL+1)=uoYL(2*nnodesuoYL+1)
uYLSecond(2*nnodesYL+2)=uoYL(2*nnodesuoYL+2)
uYLSecond(2*nnodesYL+3)=uoYL(2*nnodesuoYL+3)

!*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE SECOND INITIAL SOLUTION
CALL mesh1DLIN(0.d0,uYLSecond(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** EVENLY DISTRIBUTED MESH
DO i=1,nnodesuoYL
  uorYL(i)=uoYL(i)
ENDDO
icnt1=1
DO i=nnodesuoYL+1,2*nnodesuoYL
  uothYL(icnt1)=uoYL(i)
  icnt1=icnt1+1
ENDDO
CALL interpsol1D(uorYL,sptuoYL,nnodesuoYL,urYL,sptYL,nnodesYL)
CALL interpsol1D(uothYL,sptuoYL,nnodesuoYL,uthYL,sptYL,nnodesYL)
DO i=1,nnodesYL
  uYLSecond(i)=urYL(i)
ENDDO
icnt1=1
DO i=nnodesYL+1,2*nnodesYL
  uYLSecond(i)=uthYL(icnt1)
  icnt1=icnt1+1
ENDDO
uYLSecond(nnodesYL)=uoYL(nnodesuoYL)
uYLSecond(2*nnodesYL)=uoYL(2*nnodesuoYL)

!*** GENERATE AUGMENTED YOUNG-LAPLACE INITIAL SOLUTION
sdiff=0.d0 !*** SDIFF=NORM(uYLSecond-uYLFirst)
DO i=1,noukYL
  sdiff=sdiff+(uYLSecond(i)-uYLFirst(i))**2
ENDDO
sdiff=SQRT(sdiff)
DO i=1,noukYL
  dUdS(i)=(uYLSecond(i)-uYLFirst(i))/sdiff
  uYL(i)=uYLSecond(i)+dUdS(i)*dsc
ENDDO
MinDropletHeightInit=uYL(1)*COS(uYL(nnodesYL+1))

!*** MATRICES INITIALIZATION
DO i=1,noukYL
  xvec(i)=0.d0 !*** INITIALIZE X VECTOR
ENDDO
AMcnt=0
DO i=1,SIZE(am2)
  am2(i)=0.d0 !*** INITIALIZE JACOBIAN MATRIX PATTERN
ENDDO

!*** EIKONAL SOLUTION
solnum=10
OPEN(solnum,FILE='EikonalSolution.txt',ACTION='READ') !*** OPEN EIKONAL SOLUTION FILE
CALL Sol2DRead(solnum,xptEI,yptEI,uEI,nopEI,nnodesEI,nellEI,eltypeEI) !*** READ EIKONAL SOLUTION
alloc=100
CALL neigElem(nnodesEI,nellEI,nelnodesEI,nopEI,kindEI,neighEI,alloc) !*** FIND NEIGHBORING ELEMENTS AT EIKONAL'S DOMAIN

!*** ELECTRIC FIELD DISTRIBUTION INITIALIZATION
CALL CHDIR(ADJUSTL(Dir))
CALL defparamEL !*** READ ELECTRIC FIELD DISTIBUTION PROBLEM PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))
DO i=1,nellYL
  efieldYLGP(i)=0.d0
ENDDO

!*** COUNTERS & MATRICES INITIALIZATION
DCPcnt=0
GENcnt=0
IntCallCnt=0
dcouplflg=0
DEALLOCATE(uoYL)
ALLOCATE(uoYL(noukYL))
DO i=1,noukYL
  uoYL(i)=uYL(i)
ENDDO

!*** START PSEUDO ARC-LENGTH CONTINUATION
DO WHILE (uYL(1)*COS(uYL(nnodesYL+1))>MinDropletHeight)

!*** FLAGS & COUNTERS INITIALIZATION
  GENcnt=GENcnt+1
  IntCallCnt=IntCallCnt+1
  inr=0
  NRcnt=1 
  lastelEI=1

  SELECT CASE (GENcnt)
    CASE(1)
      SELECT CASE (JacobMatP)
        CASE (0)
!*** FIND SPARSE JACOBIAN MATRIX PATTERN FOR THE AUGMENTED YOUNG-LAPLACE
          CALL neigElem(noukYL,nellYL,nelnodesYL,nopYL,kindYL,neighYL,noukYL)
          CALL JacobPtrn(noukYL,nelnodesYL,nopYL,kindYL,neighYL)
!*** PRINT JACOBIAN MATRIX PATTERN
!         OPEN(200,FILE='JacobMat.txt')
!         WRITE(200,*)AMcntyl
!         DO i=1,AMcnt
!           WRITE(200,*)am2yl(i)
!         ENDDO
!         CLOSE(200)
        CASE (1)
!*** READ JACOBIAN MATRIX PATTERN
          OPEN(200,FILE='JacobMat.txt')
          READ(200,*)AMcnt
          DO i=1,AMcnt
            READ(200,*)am2(i)
          ENDDO
          CLOSE(200)
      END SELECT
!*** SAVE JACOBIAN MATRIX PATTERN
      ALLOCATE(am2yl(AMcnt))
      AMcntyl=AMcnt
      DO i=1,AMcnt
        am2yl(i)=am2(i)
      ENDDO

    CASE DEFAULT
!*** LOAD JACOBIAN MATRIX PATTERN
      SELECT CASE (JacobMatP)
        CASE (0)
          AMcnt=AMcntyl
          ALLOCATE(am1(AMcnt),am2(AMcnt),xvec(noukYL),res(noukYL))
          DO i=1,AMcnt
            am2(i)=am2yl(i)
          ENDDO
        CASE (1)
          OPEN(200,FILE='JacobMat.txt')
          READ(200,*)AMcnt
          ALLOCATE(am1(AMcnt),am2(AMcnt),xvec(noukYL),res(noukYL))
          DO i=1,AMcnt
            READ(200,*)am2(i)
          ENDDO
          CLOSE(200)
      END SELECT
      DO i=1,noukYL
        xvec(i)=0
        uYL(i)=uoYL(i)
      ENDDO
  END SELECT

!*** SOLVE AUGMENTED YOUNG-LAPLACE EQUATION
  DO WHILE (inr==0)!*** STAR NEWTON-RAPHSON ITERATIONS
    DO i=1,noukYL
      uYL(i)=uYL(i)+xvec(i)
    ENDDO
    NRcnt=NRcnt+1
    CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE - EVENLY DISTRIBUTED MESH
    CALL axbVparamYL !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING
    CALL ErrorNormCalc(noukYL,xvec,normsum) !*** CALCULATE ERROR NORM
    IF (normsum<=PPres.OR.GENcnt==1) THEN
      inr=1
    ENDIF
    WRITE(*,'(a,E10.3)')'AYL ERROR NORM:',normsum
    IF (normsum>1.d3) THEN
      WRITE(*,'(a)')'AYL ERROR NORM > 1.e3 - KILLING PROCESS'
      STOP
    ENDIF
  ENDDO

!*** ERROR NORM COMPUTATION FOR DECOUPLED PROBLEM
  DO i=1,noukYL
    xvec(i)=uoYL(i)-uYL(i)
  ENDDO
  CALL ErrorNormCalc(noukYL,xvec,dcp_normsum)
  IF (dcp_normsum<=PPres*10.) THEN
    dcouplflg=1
  ENDIF
  WRITE(*,'(a,E10.3)')'DECOUPLED ERROR:',dcp_normsum

!*** UPDATE THE INITIAL GUESS FOR AUGMENTED YOUNG-LAPLACE EQUATION
  DO i=1,noukYL
    uoYL(i)=uYL(i)
  ENDDO

!*** CLEAR MEMORY
  DEALLOCATE(am1,am2,xvec,res)
  
!*** EXPORT DROPLET GEOMETRY
  CALL geogenEL
  
!*** 2D MESHING
  SELECT CASE (meshGen)
    CASE(1) !*** GENERATE UNSTRUCTURED MESH USING TRIANGLE
      SELECT CASE (eltypeEL)
        CASE(1)
          CALL SYSTEM("triangle -A -a -q25 -Q dropletgeom.poly")
        CASE(2)
          CALL SYSTEM("triangle -A -a -q25 -Q -o2 dropletgeom.poly")
      END SELECT
    CASE(2)  !*** GENERATE UNSTRUCTURED MESH USING GMSH
      WRITE(ElemOrder,'(i1)')eltypeEL
      CALL SYSTEM("gmsh dropletgeom.geo -2 -v 2 -algo del2d -rand 1.e-9 -clmax 1.e-1 -order " //ElemOrder)
  END SELECT
  
!*** READ 2D MESH FILE
  CALL readmeshEL

!*** ALLOCATE MEMORY
  ALLOCATE(am1(50*noukEL),am2(50*noukEL),xvec(noukEL),uEL(noukEL),res(noukEL))

!*** FLAGS AND MATRICES INITIALIZATION
  inr=0
  NRcnt=1
  AMcnt=0
  DO i=1,SIZE(am2)
    am2(i)=0.d0 !*** INITIALIZE JACOBIAN MATRIX PATTERN
  ENDDO
  DO i=1,noukEL
    uEL(i)=bcEL(i) !*** INITIALIZE SOLUTION
    xvec(i)=0.d0 !*** INITIALIZE X VECTOR
  ENDDO

!*** FIND SPARSE JACOBIAN MATRIX PATTERN FOR ELECTRIC FIELD DISTRIBUTION PROBLEM
  alloc=100
  CALL neigElem(nnodesEL,nellEL,nelnodesEL,nopEL,kindEL,neighEL,alloc)
  CALL JacobPtrn(noukEL,nelnodesEL,nopEL,kindEL,neighEL)

!*** SOLVE ELECTRIC FIELD DISTRIBUTION EQUATIONS
  DO WHILE (inr==0)!*** START NEWTON-RAPHSON ITERATIONS
    DO i=1,noukEL
      uEL(i)=uEL(i)+xvec(i)
    ENDDO
    NRcnt=NRcnt+1
    CALL axbEL !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING
    CALL ErrorNormCalc(noukEL,xvec,normsum) !*** CALCULATE ERROR NORM
    IF (normsum<=PPres) THEN
      inr=1
    ENDIF
    WRITE(*,'(a,E10.3)')'EL ERROR NORM:',normsum
    IF (normsum>1.d5) THEN
      WRITE(*,'(a)')'EL ERROR NORM > 1.e5 - KILLING PROCESS...'
      EXIT
    ENDIF
  ENDDO

!*** ELECTRID FIELD EXTRUSION
  CALL efextr

!*** PRINT RESULTS
  IF (dcouplflg==1) THEN
    DCPcnt=DCPcnt+1
    CALL elfieldElms !*** ELECTRIC FIELD CALCULATION
    CALL printres(DCPcnt)
    WRITE(22,*)VltEL,uYL(1)
    WRITE(*,'(a,i4,/)')'DECOUPLED PROBLEM ITERATIONS',IntCallCnt
    WRITE(*,'(a,1x,f6.2,1x,a,/)')'PARAMETRIC SOLVER PROGRESS:',(1.d0-ABS(MinDropletHeight-(uYL(1)*          &
    COS(uYL(nnodesYL+1))))/ABS(MinDropletHeight-MinDropletHeightInit))*100.d0,'%'

!*** PRODUCE AUGMENTED YOUNG-LAPLACE INITIAL SOLUTION
    DO i=1,noukYL
      uYLFirst(i)=uYLSecond(i)
      uYLSecond(i)=uYL(i)
    ENDDO
    sdiff=0.d0
    DO i=1,noukYL
      sdiff=sdiff+(uYLSecond(i)-uYLFirst(i))**2
    ENDDO
    sdiff=SQRT(sdiff)
!*** MODIFY DS STEP FOR PREUDO ARC-LENGTH CONTINUATION (DSC)
!    SELECT CASE(IntCallCnt)
!      CASE(1)
!        dsc=1.5d0*dsc
!      CASE(2)
!        dsc=1.1d0*dsc
!      CASE(3)
!        dsc=1.05d0*dsc
!      CASE(5:6)
!        dsc=dsc/1.5d0
!      CASE(7:)
!        dsc=dsc/2.d0
!    END SELECT
    DO i=1,noukYL
      dUdS(i)=(uYLSecond(i)-uYLFirst(i))/sdiff
      uYL(i)=uYLSecond(i)+dUdS(i)*dsc
    ENDDO
    DEALLOCATE(efieldELElm)
    dcouplflg=0
    IntCallCnt=0
  ENDIF

!*** CLEAR MEMORY
  DEALLOCATE(am1,am2,xvec,res,xptEL,yptEL,nopEL,ncodEL,bcEL,nellidEL,uEL,kindEL,neighEL)

ENDDO

!*** CLOSE FILES
CLOSE(22)

!*** CLEAR MEMORY
DEALLOCATE(uYL,uoYL,sptYL,dsptYL,nopYL,xptYLGP,yptYLGP,efieldYLGP,am2yl,xptEI,yptEI,uEI,nopEI,kindEI,       &
neighEI,sptYLGP)

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE axbVparamYL

!************************************************************************************************************
!*** DESCRIPTION: BOUNDARY CONDITION IMPLEMENTATION AND SYSTEM SOLVING FOR AUGNENTED YOUNG-LAPLACE EQUATION
!************************************************************************************************************

USE CommonVars
USE qsort_module
IMPLICIT NONE
INTEGER(KIND=8)::i
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
  CALL abfindVparamYL
ENDDO

!*** VOLUME CONSERVATION
res(2*nnodesYL+1)=res(2*nnodesYL+1)+Pi

!*** EXTRA EQUATION FOR THE SMAX UNKNOWN
res(2*nnodesYL+2)=-uYL(2*nnodesYL)
icnt1=((2*nnodesYL+2)-1)*noukYL+2*nnodesYL
jacob=1.d0
CALL buildjacob_nor(icnt1,jacob)

!*** PSEUDO ARC-LENGTH CONTINUATION EQUATION
DO i=1,noukYL
  res(2*nnodesYL+3)=res(2*nnodesYL+3)-(+dUdS(i))*(uYL(i)-uYLSecond(i))
ENDDO
res(2*nnodesYL+3)=res(2*nnodesYL+3)-(-dsc)
DO icnt2=1,noukYL
  jacob=dUdS(icnt2)
  icnt1=((2*nnodesYL+3)-1)*noukYL+icnt2
  CALL buildjacob_nor(icnt1,jacob)
ENDDO

!*** IMPOSE NEUMANN BOUNDARY CONDITIONS
res(1)=res(1)-(-Pi/2.d0)
res(nnodesYL)=res(nnodesYL)-(+3.d0*Pi/2.d0)

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
  CALL qsort(am2,am1)
ENDIF

!*** IMPOSE DIRICHLET BOUNDARY CONDITIONS
res(nnodesYL+1)=0.d0
DO icnt2=1,noukYL
  jacob=0.d0
  icnt1=((nnodesYL+1)-1)*noukYL+icnt2
  CALL buildjacob_dir(icnt1,jacob)
ENDDO
jacob=1.d0
icnt1=((nnodesYL+1)-1)*noukYL+nnodesYL+1
CALL buildjacob_dir(icnt1,jacob)

!*** SOLVE THE LINEAR SYSTEM USING MUMPS DIRECT SOLVER
noukLS=noukYL
CALL linsolver

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE abfindVparamYL

!************************************************************************************************************
!*** DESCRIPTION: RESIDUALS & JACOBIAN MATRICES ASSEMBLY FOR AUGMENTED YOUNG-LAPLACE EQUATION
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::n
INTEGER(KIND=8)::m
INTEGER(KIND=8)::m1
INTEGER(KIND=8)::m2
INTEGER(KIND=8)::n1
INTEGER(KIND=8)::n2
INTEGER(KIND=8)::ngl(6)
INTEGER(KIND=8)::icnt1
DOUBLE PRECISION::gpAYL
DOUBLE PRECISION::wgpAYL
DOUBLE PRECISION::r
DOUBLE PRECISION::rd
DOUBLE PRECISION::drd
DOUBLE PRECISION::th
DOUBLE PRECISION::thd
DOUBLE PRECISION::dthd
DOUBLE PRECISION::s
DOUBLE PRECISION::ds
DOUBLE PRECISION::s1
DOUBLE PRECISION::ds1
DOUBLE PRECISION::ph(2)
DOUBLE PRECISION::phd(2)
DOUBLE PRECISION::phs(2)
DOUBLE PRECISION::dphs(2)
DOUBLE PRECISION::dudrei(2)
DOUBLE PRECISION::dudthei(2)
DOUBLE PRECISION::null1
DOUBLE PRECISION::null2
DOUBLE PRECISION::ueigp
DOUBLE PRECISION::fregion
DOUBLE PRECISION::at
DOUBLE PRECISION::jacob
DOUBLE PRECISION::bne
DOUBLE PRECISION::LSact=1.d0

DO i=1,4
  ngl(i)=nopYL(nell,i)
ENDDO

!*** GAUSS POINTS AND WEIGHTS
gpAYL=0.5d0
wgpAYL=1.d0
 
CALL tsfun1DLIN(gpAYL,ph,phd) !*** CALCULATE BASIS FUNCTIONS

!*** INITIALIZE VARIABLES
r=0.d0; rd=0.d0; drd=0.d0; th=0.d0; thd=0.d0; dthd=0.d0; s=0.d0; ds=0.d0; s1=0.d0; ds1=0.d0

!***ISOPARAMETRIC MAPPING
DO n=1,2
  s=s+sptYL(ngl(n))*ph(n)
  s1=s1+sptYL(ngl(n))*phd(n) !*** S1=THE C DERIVATIVE OF S
  r=r+uYL(ngl(n))*ph(n) !*** R=RADIAL COORDINATE AT THE GAUSS POINT
  th=th+uYL(ngl(n+2))*ph(n) !*** TH=ANGULAR COORDINATE AT THE GAUSS POINT
  ds=ds+dsptYL(ngl(n))*ph(n) !*** DS=SMAX DERIVATIVE OF S
  ds1=ds1+dsptYL(ngl(n))*phd(n) !*** DS1=SMAX DERIVATIVE OF S1
ENDDO

DO n=1,2
  phs(n)=phd(n)/s1 !*** PHS=S DERIVATIVE OF PHI THE  BASIS FUNCTION
  dphs(n)=-(phd(n)/(s1**2))*ds1 !*** DPHS=SMAX DERIVATIVE OF PHS
ENDDO

DO n=1,2
  rd=rd+uYL(ngl(n))*phs(n) !*** RD=S DERIVATIVE OF R
  drd=drd+uYL(ngl(n))*dphs(n) !*** DRD=SMAX DERIVATIVE OF RD   
  thd=thd+uYL(ngl(n+2))*phs(n) !*** THD=S DERIVATIVE OF TH
  dthd=dthd+uYL(ngl(n+2))*dphs(n) !*** DTHD=SMAX DERIVATIVE OF THD
ENDDO

!*** SAVE GAUSS POINTS COORDINATES
sptYLGP(nell)=s
xptYLGP(nell)=r*SIN(th)
yptYLGP(nell)=r*COS(th)

IF (flatSolidUp==0) THEN
!*** CALCULATE EIKONAL DERIVATIVES AT YOUNG-LAPLACE NODES
  IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
!*** UNSTRUCTURED EIKONAL MESH
!  CALL eikonval(uYL(ngl(1)),uYL(ngl(3)),dudrei(1),dudthei(1),null1)
!  CALL eikonval(uYL(ngl(2)),uYL(ngl(4)),dudrei(2),dudthei(2),null1)
!*** STRUCTURED EIKONAL MESH
    CALL eikonvalstr(uYL(ngl(1)),uYL(ngl(3)),dudrei(1),dudthei(1),null1)
    CALL eikonvalstr(uYL(ngl(2)),uYL(ngl(4)),dudrei(2),dudthei(2),null1)
  ENDIF
ENDIF

!*** ELECTRIC BOND NUMBER
VltEL=uYL(2*nnodesYL+3)
bne=(eoEL*VltEL**2)/(stYL*roYL)

!*** ASSEMBLY RESIDUALS

!*** VOLUME CONSERVATION TERM
res(2*nnodesYL+1)=res(2*nnodesYL+1)-wgpAYL*s1*(r**2)*thd

DO m=1,2
  m1=ngl(m)
  m2=ngl(m+2)

!*** ARC-LENGTH EQUATION TERMS
  res(m2)=res(m2)-(+wgpAYL*s1*ph(m)*(rd**2))
  res(m2)=res(m2)-(+wgpAYL*s1*ph(m)*(r**2)*(thd**2))
  res(m2)=res(m2)-(-wgpAYL*s1*ph(m))

!*** YOUNG-LAPLACE EQUATION TERMS

!*** GRAVITY TERM
!  res(m1)=res(m1)-(wgpAYL*s1*ph(m)*GcYL*r*COS(th)*SQRT(r**2*thd**2+rd**2)) 

!*** SOLID-LIQUID INTERACTION TERMS
  IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
    SELECT CASE (flatSolidUp)
      CASE(0)
!        CALL eikonval(r,th,null1,null2,ueigp) !*** UNSTRUCTURED EIKONAL MESH
        CALL eikonvalstr(r,th,null1,null2,ueigp) !*** STRUCTURED EIKONAL MESH
      CASE(1)
        ueigp=r*COS(th) !*** ONLY FOR FLAT SOLID SURFACE
    END SELECT
!    fregion=4.d0*WParamYL*((sigmaYL/(ueigp+SmallNumbYL))**8-(sigmaYL/((ueigp+SmallNumbYL)))**6) !*** C1=8, C2=6
   fregion=4.d0*WParamYL*((sigmaYL/(ueigp+SmallNumbYL))**12-(sigmaYL/((ueigp+SmallNumbYL)))**10) !*** C1=12, C2=10
    res(m1)=res(m1)-(+LSact*wgpAYL*s1*ph(m)*(fregion*SQRT(r**2*thd**2+rd**2)))
  ENDIF

!*** CURVATURE TERMS
  IF (ATAN2(r*thd,rd).LE.0.) THEN
    at=ATAN2(r*thd,rd)+2.d0*Pi !*** CHANGE ATAN2 SIGN
  ELSE
    at=ATAN2(r*thd,rd)
  ENDIF
  res(m1)=res(m1)-(wgpAYL*s1*(ph(m)*thd-phs(m)*at))

!*** PRESSURE TERMS
  res(m1)=res(m1)-(-wgpAYL*s1*ph(m)*uYL(2*nnodesYL+1)*SQRT(r**2*thd**2+rd**2))

!*** ELECTRIC TERMS
  res(m1)=res(m1)-(-wgpAYL*s1*ph(m)*bne*efieldYLGP(nell)*SQRT(r**2*thd**2+rd**2)/2.)

!*** ASSEMBLY JACOBIAN MATRIX

!*** D(RES)/D(K) TERMS
  icnt1=(m1-1)*noukYL+(2*nnodesYL+1)
  jacob=+(-wgpAYL*s1*ph(m)*SQRT(r**2*thd**2+rd**2))
  CALL buildjacob_nor(icnt1,jacob)

!*** D(VOLUME CONSERVATION TERM)/D(R) TERMS
  icnt1=((2*nnodesYL+1)-1)*noukYL+m1
  jacob=+2.d0*wgpAYL*s1*r*thd*ph(m)
  CALL buildjacob_nor(icnt1,jacob)

!*** D(VOLUME CONSERVATION TERM)/D(TH) TERMS
  icnt1=((2*nnodesYL+1)-1)*noukYL+m2
  jacob=wgpAYL*s1*(r**2)*phs(m)
  CALL buildjacob_nor(icnt1,jacob)

!*** D(RES)/D(SMAX) TERMS
  icnt1=(m1-1)*noukYL+(2*nnodesYL+2)
!*** RES=CURVATURE TERMS
  jacob=(wgpAYL*(ph(m)*thd-phs(m)*at))*ds1                                                                  &
  +(-(phs(m)*r*rd/(r**2*thd**2+rd**2)-ph(m))*s1*wgpAYL)*dthd                                                &
  +(phs(m)*r*s1*thd*wgpAYL/(r**2*thd**2+rd**2))*drd                                                         &
  +(-s1*wgpAYL*at)*dphs(m)
!*** RES=PRESSURE TERMS
  jacob=jacob+(-wgpAYL*ph(m)*uYL(2*nnodesYL+1)*SQRT(r**2*thd**2+rd**2))*ds1                                 &
  +(-uYL(2*nnodesYL+1)*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*dthd                               &
  +(-uYL(2*nnodesYL+1)*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*drd
!*** RES=GRAVITY TERMS
!  jacob=jacob+(+wgpAYL*ph(m)*GcYL*r*COS(th)*SQRT(r**2*thd**2+rd**2))*ds1                                    &
!  +(GcYL*ph(m)*r**3*s1*thd*wgpAYL*COS(th)/SQRT(r**2*thd**2+rd**2))*dthd                                     &
!  +(GcYL*ph(m)*r*rd*s1*wgpAYL*COS(th)/SQRT(r**2*thd**2+rd**2))*drd
!*** RES=ELECTRIC TERMS
  jacob=jacob+(-0.5d0*SQRT(r**2*thd**2+rd**2)*bne*efieldYLGP(nell)*ph(m)*wgpAYL)*ds1                        &
  +(-0.5d0*bne*efieldYLGP(nell)*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*dthd                      &
  +(-0.5d0*bne*efieldYLGP(nell)*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*drd
!*** RES=SOLID-LIQUID INTERACTION TERMS
!*** C1=8, C2=6
  IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
!    jacob=jacob+((LSact*wgpAYL*ph(m)*(fregion*SQRT(r**2*thd**2+rd**2))))*ds1                                &
!    +(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*LSact*WParamYL*           &
!    ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*drd                                                         &
!    +(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*LSact*WParamYL*           &
!    ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*dthd
!*** C1=12, C2=10
    jacob=jacob+((LSact*wgpAYL*ph(m)*(fregion*SQRT(r**2*thd**2+rd**2))))*ds1                                &	
    +(-4.d0*LSact*WParamYL*ph(m)*rd*s1*wgpAYL*(sigmaYL**10/(SmallNumbYL+ueigp)**10-                         &
    sigmaYL**12/(SmallNumbYL+ueigp)**12)/SQRT(r**2*thd**2+rd**2))*drd                                       &
    +(-4.d0*LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL*(sigmaYL**10/(SmallNumbYL+ueigp)**10-                   &
    sigmaYL**12/(SmallNumbYL+ueigp)**12)/SQRT(r**2*thd**2+rd**2))*dthd
  ENDIF
  CALL buildjacob_nor(icnt1,jacob)
!*** RES=ARC-LENGTH EQUATION TERMS
  icnt1=(m2-1)*noukYL+(2*nnodesYL+2)
  jacob=+wgpAYL*ph(m)*(rd**2)*ds1                                                                           &
  +wgpAYL*s1*ph(m)*2.d0*rd*drd                                                                              &
  +wgpAYL*ph(m)*(r**2)*(thd**2)*ds1                                                                         &
  +wgpAYL*s1*ph(m)*(r**2)*2.d0*thd*dthd                                                                     &
  -wgpAYL*ph(m)*ds1
  CALL buildjacob_nor(icnt1,jacob)
!*** RES=VOLUME CONSERVATION TERM
  icnt1=((2*nnodesYL+1)-1)*noukYL+(2*nnodesYL+2)
  jacob=+wgpAYL*(r**2)*thd*ds1                                                                              &
  +wgpAYL*(r**2)*s1*dthd
  CALL buildjacob_nor(icnt1,jacob)

!*** D(RES)/D(VltEL) TERMS
  icnt1=(m1-1)*noukYL+(2*nnodesYL+3)
!*** RES=ELECTRIC TERMS
  jacob=-SQRT(r**2*thd**2+rd**2)*VltEL*efieldYLGP(nell)*eoEL*ph(m)*s1*wgpAYL/(roYL*stYL)
  CALL buildjacob_nor(icnt1,jacob)

  DO n=1,2
    n1=ngl(n)
    n2=ngl(n+2)

!*** D(RES)/D(R) TERMS

!*** RES=CURVATURE TERMS
    icnt1=(m1-1)*noukYL+n1
    jacob=+(-phs(m)*rd*s1*thd*wgpAYL/(r**2*thd**2+rd**2))*ph(n)                                             &
    +(phs(m)*r*s1*thd*wgpAYL/(r**2*thd**2+rd**2))*phs(n)
!*** RES=PRESSURE TERMS
    jacob=jacob+(-uYL(2*nnodesYL+1)*ph(m)*r*s1*thd**2*wgpAYL/SQRT(r**2*thd**2+rd**2))*ph(n)                 &
    +(-uYL(2*nnodesYL+1)*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=GRAVITY TERMS
!    jacob=jacob+(GcYL*ph(m)*r**2*s1*thd**2*wgpAYL*COS(th)/SQRT(r**2*thd**2+rd**2)                           &
!    +SQRT(r**2*thd**2+rd**2)*GcYL*ph(m)*s1*wgpAYL*COS(th))*ph(n)+                                           &
!    (GcYL*ph(m)*r**3*s1*thd*wgpAYL*COS(th)/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=ELECTRIC TERMS
     jacob=jacob+(-0.5d0*bne*efieldYLGP(nell)*ph(m)*r*s1*thd**2*wgpAYL/SQRT(r**2*thd**2+rd**2))*ph(n)       &
     +(-0.5d0*bne*efieldYLGP(nell)*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=SOLID-LIQUID INTERACTION TERMS
    IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
      SELECT CASE (flatSolidUp)
        CASE(0)
!*** C1=8, C2=6
!          jacob=jacob+(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*         &
!           LSact*WParamYL*ph(m)*r*s1*thd**2*wgpAYL/SQRT(r**2*thd**2+rd**2))*ph(n)+                          &
!          (-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*                     &
!          LSact*WParamYL*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)+                                &
!          (8.d0*(3.d0*sigmaYL**6/(SmallNumbYL+ueigp)**7-4.d0*sigmaYL**8/(SmallNumbYL+ueigp)**9)*            &
!          SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*dudrei(n)*ph(n)
!*** C1=12, C2=10
          jacob=jacob+(-4.d0*LSact*WParamYL*ph(m)*r*s1*thd**2*wgpAYL*(sigmaYL**10/(SmallNumbYL+ueigp)**10   &
          -sigmaYL**12/(SmallNumbYL+ueigp)**12)/SQRT(r**2*thd**2+rd**2))*ph(n)+(-4.d0*LSact*WParamYL*       &
          ph(m)*rd*s1*wgpAYL*(sigmaYL**10/(SmallNumbYL+ueigp)**10-sigmaYL**12/(SmallNumbYL+                 &
          ueigp)**12)/SQRT(r**2*thd**2+rd**2))*phs(n)+(8.d0*SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*         &
          ph(m)*s1*wgpAYL*(5*sigmaYL**10/(SmallNumbYL+ueigp)**11-6*sigmaYL**12/(SmallNumbYL+                &
          ueigp)**13))*dudrei(n)*ph(n)
        CASE(1)
!*** FLAT SOLID SURFACE	WITH C1=8, C2=6
!          jacob=jacob+(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*         &
!          LSact*WParamYL*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)+                                &
!          (-4.d0*(sigmaYL**6/(r*COS(th)+SmallNumbYL)**6-                                                    &
!          sigmaYL**8/(r*COS(th)+SmallNumbYL)**8)*LSact*WParamYL*ph(m)*r*s1*thd**2*                          &
!          wgpAYL/SQRT(r**2*thd**2+rd**2)+8.d0*(3.d0*sigmaYL**6*COS(th)/(r*COS(th)+                          &
!          SmallNumbYL)**7-4.d0*sigmaYL**8*COS(th)/(r*COS(th)+SmallNumbYL)**9)*SQRT(r**2*thd**2              &
!          +rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*ph(n)
!*** FLAT SOLID SURFACE	WITH C1=12, C2=10
          jacob=jacob+(-4.d0*LSact*WParamYL*ph(m)*r*s1*thd**2*wgpAYL*(sigmaYL**10/(r*COS(th)                &
          +SmallNumbYL)**10-sigmaYL**12/(r*COS(th)+SmallNumbYL)**12)/SQRT(r**2*thd**2+rd**2)                &
          +8.d0*SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL*(5.d0*sigmaYL**10*COS(th)/(r*COS(th) &
          +SmallNumbYL)**11-6.d0*sigmaYL**12*COS(th)/(r*COS(th)+SmallNumbYL)**13))*ph(n)+                   &
          (-4.d0*LSact*WParamYL*ph(m)*rd*s1*wgpAYL*(sigmaYL**10/(r*COS(th)+SmallNumbYL)**10                 &
          -sigmaYL**12/(r*COS(th)+SmallNumbYL)**12)/SQRT(r**2*thd**2+rd**2))*phs(n)
      END SELECT
    ENDIF
    CALL buildjacob_nor(icnt1,jacob)
!*** RES=ARC-LENGTH EQUATION TERMS
    icnt1=(m2-1)*noukYL+n1
    jacob=2.d0*ph(m)*r*s1*thd**2*wgpAYL*ph(n)                                                               &
    +2.d0*ph(m)*rd*s1*wgpAYL*phs(n)
    CALL buildjacob_nor(icnt1,jacob)

!*** D(RES)/D(TH) TERMS

!*** RES=CURVATURE TERMS
    icnt1=(m1-1)*noukYL+n2
    jacob=+(-(phs(m)*r*rd/(r**2*thd**2+rd**2)-ph(m))*s1*wgpAYL)*phs(n)
!*** RES=PRESSURE TERMS
    jacob=jacob+(-uYL(2*nnodesYL+1)*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=GRAVITY TERMS
!    jacob=jacob+(-SQRT(r**2*thd**2+rd**2)*GcYL*ph(m)*r*s1*wgpAYL*SIN(th))*ph(n)                             &
!    +(GcYL*ph(m)*r*rd*s1*wgpAYL*COS(th)/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=ELECTRIC TERMS
     jacob=jacob+(-0.5d0*bne*efieldYLGP(nell)*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
!*** RES=SOLID-LIQUID INTERACTION TERMS
    IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
      SELECT CASE (flatSolidUp)
        CASE(0)
!*** C1=8, C2=6
!          jacob=jacob+(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*         &
!          LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)                           &
!          +(8.d0*(3.d0*sigmaYL**6/(SmallNumbYL+ueigp)**7-4.d0*sigmaYL**8/(SmallNumbYL+ueigp)**9)*           &
!          SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*dudthei(n)*ph(n)
!*** C1=12, C2=10
          jacob=jacob+(-4.d0*LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL*(sigmaYL**10/(SmallNumbYL+ueigp)**10   &
          -sigmaYL**12/(SmallNumbYL+ueigp)**12)/SQRT(r**2*thd**2+rd**2))*phs(n)+                            &
          (8.d0*SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL*(5.d0*sigmaYL**10/(SmallNumbYL+      &
          ueigp)**11-6.d0*sigmaYL**12/(SmallNumbYL+ueigp)**13))*dudthei(n)*ph(n)
        CASE(1)
!*** FLAT SOLID SURFACE	WITH C1=8, C2=6
!          jacob=jacob+(-4.d0*(sigmaYL**6/(SmallNumbYL+ueigp)**6-sigmaYL**8/(SmallNumbYL+ueigp)**8)*         &
!          LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)+                          &
!          (-8.d0*(3.d0*r*sigmaYL**6*SIN(th)/(r*COS(th)+SmallNumbYL)**7-                                     &
!          4.d0*r*sigmaYL**8*SIN(th)/(r*COS(th)+SmallNumbYL)**9)*SQRT(r**2*thd**2+rd**2)*LSact*              &
!          WParamYL*ph(m)*s1*wgpAYL)*ph(n)
!*** FLAT SOLID SURFACE	WITH C1=12, C2=10
          jacob=jacob+(-8.d0*SQRT(r**2*thd**2+rd**2)*LSact*WParamYL*ph(m)*s1*wgpAYL*                        &
          (5.d0*r*sigmaYL**10*SIN(th)/(r*COS(th)+SmallNumbYL)**11-6.d0*r*sigmaYL**12*SIN(th)/(r*COS(th)     &
          +SmallNumbYL)**13))*ph(n)+(-4.d0*LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL*(sigmaYL**10/(r*COS(th)+ &
          SmallNumbYL)**10-sigmaYL**12/(r*COS(th)+SmallNumbYL)**12)/SQRT(r**2*thd**2+rd**2))*phs(n)
      END SELECT
    ENDIF
    CALL buildjacob_nor(icnt1,jacob)
!*** RES=ARC-LENGTH EQUATION TERMS
    icnt1=(m2-1)*noukYL+n2
    jacob=2.d0*ph(m)*r**2*s1*thd*wgpAYL*phs(n)
    CALL buildjacob_nor(icnt1,jacob)

  ENDDO

ENDDO

END SUBROUTINE

!************************************************************************************************************
