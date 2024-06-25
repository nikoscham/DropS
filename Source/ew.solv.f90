!************************************************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com) @ NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!************************************************************************************************************

SUBROUTINE EWSolver

!************************************************************************************************************
!*** DESCRIPTION: SOLVE FORCE BALANCE (AUGMENTED YOUNG-LAPLACE) AND ELECTRIC FIELD DISTRIBUTION EQUATIONS
!************************************************************************************************************

USE CommonVars
USE initsol
USE nodnum
USE defpar
IMPLICIT NONE
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::dcouplflg
INTEGER(KIND=8)::GENcnt
INTEGER(KIND=8)::IntCallCnt
INTEGER(KIND=8)::i
INTEGER(KIND=8)::i1
INTEGER(KIND=8)::i2
INTEGER(KIND=8)::AMcntyl
INTEGER(KIND=8)::inct1
INTEGER(KIND=8)::alloc
INTEGER(KIND=8)::solnum
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::am2yl
DOUBLE PRECISION::normsum
DOUBLE PRECISION::dcp_normsum
DOUBLE PRECISION::lsdistance
DOUBLE PRECISION::VltELinit
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::urYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uthYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uorYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uothYL
CHARACTER(LEN=1)::ElemOrder

!*** SOLVER INITIALIZATION
IF (SolverID.NE.6) THEN
  OPEN(22,FILE='ResCon.txt') !*** OPEN CONTINUATION RESULTS FILE
  WRITE(22,'(a)')'DropS RESULTS'
  WRITE(22,'(a)')'VOLTAGE, MAXIMUM DROPLET HEIGHT'
ENDIF
CALL CHDIR(ADJUSTL(Dir))
CALL defparamYL !*** READ AYL PARAMETERS
CALL defparamEI !*** READ EI PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))

!*** ALLOCATE MEMORY
ALLOCATE (uYL(noukYL),uoYL(noukYL),sptYL(nnodesYL),dsptYL(nnodesYL),xvec(noukYL),nopYL(nellYL,7),           &
res(noukYL),am1(50*noukYL),am2(50*noukYL),xptYLGP(nellYL),yptYLGP(nellYL),efieldYLGP(nellYL),               &
sptuoYL(nnodesYL),urYL(nnodesYL),uthYL(nnodesYL),uorYL(nnodesYL),uothYL(nnodesYL),sptYLGP(nellYL),          &
YLnodes(nnodesYL))

!*** READ AUGMENTED YOUNG-LAPLACE INITIAL SOLUTION
solnum=10
OPEN(solnum,FILE='data_in_uo.txt',ACTION='READ')
CALL Sol1DRead(solnum,noukuoYL,nnodesuoYL,uoYL,sptuoYL)
uYL(2*nnodesYL+1)=uoYL(2*nnodesuoYL+1)
uYL(2*nnodesYL+2)=uoYL(2*nnodesuoYL+2)

!*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE
CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** EVENLY DISTRIBUTED MESH

!*** NODAL NUMBERING
CALL nodnumbAYL

!*** GENERATE AUGMENTED YOUNG-LAPLACE INITIAL SOLUTION FOR CURRENT MESH
IF (nellYL.NE.nelluoYL) THEN
  DO i=1,nnodesuoYL
    uorYL(i)=uoYL(i)
  ENDDO
  inct1=1
  DO i=nnodesuoYL+1,2*nnodesuoYL
    uothYL(inct1)=uoYL(i)
    inct1=inct1+1
  ENDDO
  CALL interpsol1D(uorYL,sptuoYL,nnodesuoYL,urYL,sptYL,nnodesYL)
  CALL interpsol1D(uothYL,sptuoYL,nnodesuoYL,uthYL,sptYL,nnodesYL)
  DO i=1,nnodesYL
    uYL(i)=urYL(i)
  ENDDO
  inct1=1
  DO i=nnodesYL+1,2*nnodesYL
    uYL(i)=uthYL(inct1)
    inct1=inct1+1
  ENDDO
  uYL(1)=uoYL(1)
  uYL(nnodesYL+1)=uoYL(nnodesuoYL+1)
  uYL(nnodesYL)=uoYL(nnodesuoYL)
  uYL(2*nnodesYL)=uoYL(2*nnodesuoYL)
ELSE
  DO i=1,noukYL
    uYL(i)=uoYL(i)
  ENDDO
ENDIF

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
CALL neigElem(nnodesEI,nellEI,nelnodesEI,nopEI,kindEI,neighEI,alloc) !*** FIND NEIGHBORING ELEMENTS

!*** ELECTRIC FIELD DISTRIBUTION INITIALIZATION
CALL CHDIR(ADJUSTL(Dir))
CALL defparamEL !*** READ ELECTRIC FIELD DISTRIBUTION PROBLEM PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))
VltELinit=VltEL
DO i=1,nellYL
  efieldYLGP(i)=0.d0
ENDDO

!*** COUNTERS INITIALIZATION
DCPcnt=0
IntCallCnt=0
GENcnt=0
dcouplflg=0

!*** START ZERO ORDER CONTINUATION
DO WHILE (VltEL<=VltMaxEL)

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
  DO WHILE (inr==0) !*** START NEWTON-RAPHSON ITERATIONS
    DO i=1,noukYL
      uYL(i)=uYL(i)+xvec(i)
    ENDDO
    NRcnt=NRcnt+1
    CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE - EVENLY DISTRIBUTED MESH
    CALL axbYL !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING
    CALL ErrorNormCalc(noukYL,xvec,normsum) !*** CALCULATE ERROR NORM
    IF (normsum<=PPres.OR.GENcnt==1.OR.SolverID==6) THEN
      inr=1
    ENDIF
    WRITE(*,'(a,E10.3)')'AYL ERROR NORM:',normsum
    IF (normsum>1.d3) THEN
      WRITE(*,'(a)')'AYL ERROR NORM > 1.e3 - KILLING PROCESS'
      EXIT
    ENDIF
  ENDDO

!*** END OF TIME STEP
  IF (SolverID==6.AND.DCPcnt==1) THEN
!*** PRINT OUTPUT FILE
    OPEN(200,FILE='data_out.txt')
    DO i=1,noukYL
      WRITE(200,*)res(i)
    ENDDO
    DO i=1,nnodesYL
      WRITE(200,*)sptYL(i)
    ENDDO
    CLOSE(200)
!*** PRINT JACOBIAN MATRIX
    IF (printPrecond==1) THEN
      OPEN(200,FILE='Precond.txt')
      DO i=1,AMcnt
        i1=INT(am2(i)/noukYL+1.d0,8)
        i2=INT(MOD(am2(i),noukYL),8)
        IF (i2==0) THEN
          i1=i1-1
          i2=noukYL
        ENDIF
        WRITE(200,*)i1,i2,am1(i)
      ENDDO
      CLOSE(200)
    ENDIF
    EXIT !*** TERMINATE SOLVER
  ENDIF

!*** ERROR NORM COMPUTATION FOR DECOUPLED PROBLEM
  DO i=1,noukYL
    xvec(i)=uoYL(i)-uYL(i)
  ENDDO
  CALL ErrorNormCalc(noukYL,xvec,dcp_normsum)
  IF (dcp_normsum<=PPres*10.d0) THEN
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
  DO WHILE (inr==0) !*** START NEWTON-RAPHSON ITERATIONS
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
  IF (dcouplflg==1.AND.SolverID.NE.6) THEN
    DCPcnt=DCPcnt+1
    CALL elfieldElms !*** ELECTRIC FIELD CALCULATION
    CALL printres(DCPcnt)
    WRITE(22,*)VltEL,uYL(1)
    DEALLOCATE(efieldELElm)
    WRITE(*,'(a,i4,/)')'DECOUPLED PROBLEM ITERATIONS:',IntCallCnt
    WRITE(*,'(a,1x,f6.2,1x,a,/)')'PARAMETRIC SOLVER PROGRESS:',(1.d0-(VltMaxEL-VltEL)/(VltMaxEL-VltELinit)) &
    *100.d0,'%'
    VltEL=VltEL+VltStepEL !*** UPDATE THE CONTINUATION PARAMETER
    dcouplflg=0
    IntCallCnt=0
  ELSE IF (dcouplflg==1.AND.SolverID==6) THEN
    CALL elfieldElms !*** ELECTRIC FIELD CALCULATION
    CALL printres(DCPcnt)
    DCPcnt=DCPcnt+1
  ENDIF

!*** CLEAR MEMORY
  DEALLOCATE(am1,am2,xvec,res,xptEL,yptEL,nopEL,ncodEL,bcEL,nellidEL,uEL,kindEL,neighEL)

ENDDO

!*** CLOSE FILES
IF (SolverID.NE.6) THEN
  CLOSE(22)
ENDIF

!*** CLEAR MEMORY
DEALLOCATE(uYL,uoYL,sptYL,dsptYL,nopYL,xptYLGP,yptYLGP,efieldYLGP,am2yl,xptEI,yptEI,uEI,nopEI,kindEI,       &
neighEI,sptYLGP)

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE axbEL

!************************************************************************************************************
!*** DESCRIPTION: BOUNDARY CONDIATION IMPLEMENTATION AND SYSTEM SOLVING FOR ELECTRIC FIELD DISTIBUTION EQUATIONS
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::k
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::i1
INTEGER(KIND=8)::icnt1
DOUBLE PRECISION::jacob

!*** MATRICES INITIALIZATION
DO i=1,noukEL
  res(i)=0.d0
ENDDO
DO i=1,SIZE(am1)
  am1(i)=0.d0
ENDDO

!*** ASSEMBLY RESIDUALS (RES) AND JACOBIAN (AM) MATRICES
DO nell=1,nellEL
  CALL abfindEL
ENDDO

!*** IMPOSE DIRICHLET BOUNDARY CONDITIONS
DO i=1,noukEL
  IF (ncodEL(i)==1) THEN
    res(i)=bcEL(i)-uEL(i)
    jacob=0.d0
    DO j=1,kindEL(i)
      DO k=1,nelnodesEL
        i1=nopEL(neighEL(i,j),k)
        icnt1=(i-1)*noukEL+i1
        CALL buildjacob_dir(icnt1,jacob)
      ENDDO
    ENDDO
    jacob=1.d0
    icnt1=(i-1)*noukEL+i
    CALL buildjacob_dir(icnt1,jacob)
  ENDIF
ENDDO

!*** SOLVE THE LINEAR SYSTEM USING MUMPS DIRECT SOLVER
noukLS=noukEL
CALL linsolver

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE abfindEL

!************************************************************************************************************
!*** DESCRIPTION: RESIDUALS & JACOBIAN MATRICES ASSEMBLY FOR ELECTRIC FIELD DISTRIBUTION EQUATIONS
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::k
INTEGER(KIND=8)::l
INTEGER(KIND=8)::n
INTEGER(KIND=8)::m
INTEGER(KIND=8)::m1
INTEGER(KIND=8)::n1
INTEGER(KIND=8)::icnt1
INTEGER(KIND=8)::ngl(6)
DOUBLE PRECISION::x
DOUBLE PRECISION::y
DOUBLE PRECISION::x1
DOUBLE PRECISION::y1
DOUBLE PRECISION::x2
DOUBLE PRECISION::y2
DOUBLE PRECISION::dudx
DOUBLE PRECISION::dudy
DOUBLE PRECISION::dett
DOUBLE PRECISION::jacob
DOUBLE PRECISION::er
DOUBLE PRECISION::phi(6)
DOUBLE PRECISION::phic(6)
DOUBLE PRECISION::phie(6)
DOUBLE PRECISION::tphx(6)
DOUBLE PRECISION::tphy(6)

DO i=1,nelnodesEL
  ngl(i)=nopEL(nell,i)
ENDDO

IF (nellidEL(nell)==DielSurffl) THEN !*** DIELECTRIC DOMAIN
  er=edEL
ELSE IF (nellidEL(nell)==AirSurffl) THEN !*** AIR DOMAIN
  er=esEL
ENDIF

!*** LOOP OVER GAUSS POINTS
DO k=1,NumgpsEL
  DO l=1,NumgpsEL

!*** CALCULATE BASIS FUNCTIONS
    SELECT CASE (eltypeEL)
    CASE(1)
      CALL tsfun2DLIN(gpEL(l),gpEL(l),phi,phic,phie)
    CASE(2)
      CALL tsfun2DQUAD(gpEL(k),gpEL(l),phi,phic,phie)
    END SELECT

!*** INITIALIZE VARIABLES
    x=0.d0; y=0.d0; x1=0.d0; x2=0.d0; y1=0.d0; y2=0.d0; dudx=0.d0; dudy=0.d0

!***ISOPARAMETRIC MAPPING
    DO n=1,nelnodesEL
      x=x+xptEL(ngl(n))*phi(n)
      y=y+yptEL(ngl(n))*phi(n)
      x1=x1+xptEL(ngl(n))*phic(n)
      x2=x2+xptEL(ngl(n))*phie(n)
      y1=y1+yptEL(ngl(n))*phic(n)
      y2=y2+yptEL(ngl(n))*phie(n)
    ENDDO
    dett=x1*y2-x2*y1
    DO n=1,nelnodesEL
      tphx(n)=(y2*phic(n)-y1*phie(n))/dett
      tphy(n)=(x1*phie(n)-x2*phic(n))/dett
      dudx=dudx+uEL(ngl(n))*tphx(n)
      dudy=dudy+uEL(ngl(n))*tphy(n)
    ENDDO

!*** ASSEMBLY RESIDUALS
    DO n=1,nelnodesEL
      n1=ngl(n)
      res(n1)=res(n1)+er*wgpEL(k)*wgpEL(l)*dett*tphx(n)*dudx
      res(n1)=res(n1)+er*wgpEL(k)*wgpEL(l)*dett*tphy(n)*dudy
!*** ASSEMBLY JACOBIAN
      DO m=1,nelnodesEL
        m1=ngl(m)
        icnt1=(n1-1)*noukEL+m1
        jacob=-er*wgpEL(k)*wgpEL(l)*dett*(tphx(n)*tphx(m)+tphy(n)*tphy(m))
        CALL buildjacob_nor(icnt1,jacob)
      ENDDO
    ENDDO

  ENDDO
ENDDO

END SUBROUTINE

!************************************************************************************************************
