!************************************************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com) @ NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!************************************************************************************************************

SUBROUTINE YLSolver

!************************************************************************************************************
!*** DESCRIPTION: AUGMENTED YOUNG-LAPLACE EQUATION SOLVER
!************************************************************************************************************

USE CommonVars
USE initsol
USE nodnum
USE defpar
IMPLICIT NONE
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::GENcnt
INTEGER(KIND=8)::i
INTEGER(KIND=8)::AMcntyl
INTEGER(KIND=8)::inct1
INTEGER(KIND=8)::alloc
INTEGER(KIND=8)::solnum
INTEGER(KIND=8)::TimeCnt
INTEGER(KIND=8)::its
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::am2yl
DOUBLE PRECISION::normsum
DOUBLE PRECISION::normsumInit
DOUBLE PRECISION::lsdistance
DOUBLE PRECISION::WParamYLinit
DOUBLE PRECISION::TStepInit
DOUBLE PRECISION::SChange
DOUBLE PRECISION::NumJacPert
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::urYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uthYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uorYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::uothYL
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::sptYLp
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Tth
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Tr
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Tcurv

!*** SOLVER INITIALIZATION
IF (SolverID.NE.7) THEN
  OPEN(22,FILE='ResCon.txt') !*** OPEN CONTINUATION RESULTS FILE
  WRITE(22,'(a)')'DropS RESULTS'
  WRITE(22,'(a)')'SOLUTION NUMBER, MATERIAL WETTABILITY (YOUNG CONTACT ANGLE), MAXIMUM DROPLET HEIGHT'
ENDIF
CALL CHDIR(ADJUSTL(Dir))
CALL defparamYL !*** READ AYL PARAMETERS
CALL defparamEI !*** READ EI PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))
WParamYLinit=WParamYL
TStepInit=TStep
cVol=1.d0 !*** DROPLET VOLUME MULTIPLIER

!*** ALLOCATE MEMORY
ALLOCATE (uYL(noukYL),uoYL(noukYL),sptYL(nnodesYL),dsptYL(nnodesYL),xvec(noukYL),nopYL(nellYL,7),           &
res(noukYL),am1(50*noukYL),am2(50*noukYL),xptYLGP(nellYL),yptYLGP(nellYL),efieldYLGP(nellYL),               &
sptuoYL(nnodesYL),urYL(nnodesYL),uthYL(nnodesYL),uorYL(nnodesYL),uothYL(nnodesYL),sptYLGP(nellYL),          &
YLnodes(nnodesYL),sptYLp(nnodesYL))

!*** READ INITIAL SOLUTION
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
CALL neigElem(nnodesEI,nellEI,nelnodesEI,nopEI,kindEI,neighEI,alloc) !*** FIND NEIGHBORING ELEMENTS AT EIKONAL'S DOMAIN

!*** ELECTRIC FIELD DISTRIBUTION INITIALIZATION
CALL CHDIR(ADJUSTL(Dir))
CALL defparamEL !*** READ ELECTRIC FIELD DISTIBUTION PROBLEM PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))
DO i=1,nellYL
  efieldYLGP(i)=0.d0
ENDDO

!*** COUNTERS INITIALIZATION
GENcnt=0
DEALLOCATE(uoYL)
ALLOCATE(uoYL(noukYL))
DO i=1,noukYL
  uoYL(i)=uYL(i)
ENDDO

!*** START ZERO ORDER CONTINUATION
DO WHILE (WParamYL<WParamMaxYL) !*** WETTING PARAMETER
!DO WHILE (cVol<10.d0) !*** DROPLET VOLUME

!*** FLAGS & COUNTERS INITIALIZATION
  GENcnt=GENcnt+1
  inr=0
  its=0
  NRcnt=1
  lastelEI=1
  TimeCnt=0
  TStep=TStepInit

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
  DO WHILE (its==0) !*** LOOP OVER PSEUDO-TIME
    DO WHILE (inr==0)!*** START NEWTON-RAPHSON ITERATIONS
      DO i=1,noukYL
        uYL(i)=uYL(i)+xvec(i)
      ENDDO
      NRcnt=NRcnt+1
      SELECT CASE (AdaMesh) !*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE
        CASE(0)
          CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2),nellYL,sptYL,dsptYL) !*** GENERATE MESH FOR THE AUGMENTED YOUNG-LAPLACE - EVENLY DISTRIBUTED MESH
        CASE(1)
          CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2),nellYL,sptYL,dsptYL)
          ALLOCATE(Tth(nnodesYL),Tr(nnodesYL),Tcurv(nnodesYL))
          DO i=1,nnodesYL
            Tr(i)=uYL(i)
            Tth(i)=uYL(nnodesYL+i)
          ENDDO
          CALL NumCurvCalc(nnodesYL,sptYL,Tth,Tr,Tcurv)
          Tcurv=REAL(DCPcnt)*5.e-3*Tcurv
          CALL mesh1DADA(nellYL,sptYL,Tth,Tr,Tcurv) !*** ADAPTIVE MESH
          DO i=1,nnodesYL
            SChange=NumJacPert(sptYL(i))
            CALL mesh1DLIN(0.d0,uYL(2*nnodesYL+2)+SChange,nellYL,sptYLp,dsptYL)
            CALL mesh1DADA(nellYL,sptYLp,Tth,Tr,Tcurv) !*** sptYLp = PERTUBATED ARC-LENGTH COORDINATES
            dsptYL(i)=(sptYLp(i)-sptYL(i))/(SChange) !*** FORWARD DIFFERENCE
            sptYLp(i)=0.d0
          ENDDO
          DEALLOCATE(Tth,Tr,Tcurv)
      END SELECT
      CALL axbYL !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING
      CALL ErrorNormCalc(noukYL,xvec,normsum) !*** CALCULATE ERROR NORM
      IF (normsum<=PPres) THEN
        inr=1
      ENDIF
      WRITE(*,'(a,E10.3)')'AYL ERROR NORM:',normsum
      IF (normsum>1.d3) THEN !*** CHECK INITIAL ERROR NORM
        WRITE(*,'(a)')'AYL ERROR NORM>1.e3 - KILLING PROCESS...'
        STOP
      ENDIF
    ENDDO
    WRITE(*,'(a,i4,/)')'NEWTON ITERATIONS:',NRcnt-1
    SELECT CASE (TraCon)
      CASE(1)
        TimeCnt=TimeCnt+1 !*** UPDATE TIME FOR PSEUDO-TRANSIENT CONTINUATION
        CALL ErrorNormCalc(noukYL,uYL-uoYL,normsum)
        IF (TimeCnt==1) normsumInit=normsum
        WRITE(*,'(a,E10.3)')'PSEUDO-TRANSIENT CONTINUATION ERROR NORM:',normsum
        IF (TStep<=TStepMax) TStep=TStepInit*(normsumInit/normsum) !*** UPDATE TIME STEP FOR PSEUDO-TRANSIENT CONTINUATION
        WRITE(*,'(a,i4,/)')'PSEUDO-TIME STEP:',TimeCnt
    END SELECT
    IF (normsum<=PPres*10.OR.TraCon==0.OR.SolverID==7) THEN
      its=1
      TimeCnt=0
    ELSE
      DO i=1,noukYL
        uoYL(i)=uYL(i)
      ENDDO
      inr=0
      NRcnt=1
    ENDIF
  ENDDO

!*** PRINT RESULTS
  CALL printres(GENcnt)
  CALL thYoungCalc(thYoungYL,WParamYL)
  IF (SolverID==7) THEN !*** END OF TIME STEP IN CASE OF TIME STEP SOLVER
    EXIT !*** TERMINATE SOLVER
  ENDIF
  WRITE(22,*)GENcnt,thYoungYL,uYL(1) !*** WETTING PARAMETER
!  WRITE(22,*)GENcnt,cVol,uYL(1) !*** DROPLET VOLUME
  WRITE(*,'(a,1x,f6.2,1x,a,/)')'PARAMETRIC SOLVER PROGRESS:',(1.d0-(WParamMaxYL-WParamYL)                   &
  /(WParamMaxYL-WParamYLinit))*100.d0,'%'

!*** UPDATE THE CONTINUATION PARAMETER
  WParamYL=WParamYL+WParamStepYL !*** WETTING PARAMETER
!  cVol=cVol+0.1d0 !*** DROPLET VOLUME
!*** UPDATE THE INITIAL GUESS FOR AUGMENTED YOUNG-LAPLACE EQUATION
  DO i=1,noukYL
    uoYL(i)=uYL(i)
  ENDDO

!*** CLEAR MEMORY
  DEALLOCATE(am1,am2,xvec,res)

ENDDO

!*** CLOSE FILES
CLOSE(22)

!*** CLEAR MEMORY
DEALLOCATE(uYL,uoYL,sptYL,dsptYL,nopYL,xptYLGP,yptYLGP,efieldYLGP,am2yl,xptEI,yptEI,uEI,nopEI,kindEI,       &
neighEI,sptYLGP)

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE axbYL

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
  CALL abfindYL
ENDDO

!*** VOLUME CONSERVATION
res(2*nnodesYL+1)=res(2*nnodesYL+1)+cVol*Pi

!*** EXTRA EQUATION FOR THE SMAX UNKNOWN
res(2*nnodesYL+2)=-uYL(2*nnodesYL)
icnt1=((2*nnodesYL+2)-1)*noukYL+2*nnodesYL
jacob=1.d0
CALL buildjacob_nor(icnt1,jacob)

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

SUBROUTINE abfindYL

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
DOUBLE PRECISION::ro
DOUBLE PRECISION::rd
DOUBLE PRECISION::drd
DOUBLE PRECISION::th
DOUBLE PRECISION::tho
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
ro=0.d0; tho=0.d0

!*** ISOPARAMETRIC MAPPING
DO n=1,2
  s=s+sptYL(ngl(n))*ph(n) !*** ARC-LENGTH OF THE EFFECTIVELY ONE-DIMENSIONAL DROPLET SURFACE
  s1=s1+sptYL(ngl(n))*phd(n) !*** S1=THE C DERIVATIVE OF S
  r=r+uYL(ngl(n))*ph(n) !*** R=RADIAL COORDINATE AT THE GAUSS POINT
  ro=ro+uoYL(ngl(n))*ph(n)
  th=th+uYL(ngl(n+2))*ph(n) !*** TH=ANGULAR COORDINATE AT THE GAUSS POINT
  tho=tho+uoYL(ngl(n+2))*ph(n)
  ds=ds+dsptYL(ngl(n))*ph(n) !*** DS=SMAX DERIVATIVE OF S
  ds1=ds1+dsptYL(ngl(n))*phd(n) !*** DS1=SMAX DERIVATIVE OF S1
ENDDO

DO n=1,2
  phs(n)=phd(n)/s1 !*** PHS=S DERIVATIVE OF THE BASIS FUNCTION
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
  CALL eikonval(uYL(ngl(1)),uYL(ngl(3)),dudrei(1),dudthei(1),null1)
  CALL eikonval(uYL(ngl(2)),uYL(ngl(4)),dudrei(2),dudthei(2),null1)
!*** STRUCTURED EIKONAL MESH
!    CALL eikonvalstr(uYL(ngl(1)),uYL(ngl(3)),dudrei(1),dudthei(1),null1)
!    CALL eikonvalstr(uYL(ngl(2)),uYL(ngl(4)),dudrei(2),dudthei(2),null1)
  ENDIF
ENDIF

!*** ELECTRIC BOND NUMBER
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
        CALL eikonval(r,th,null1,null2,ueigp) !*** UNSTRUCTURED EIKONAL MESH
!        CALL eikonvalstr(r,th,null1,null2,ueigp) !*** STRUCTURED EIKONAL MESH
      CASE(1)
        ueigp=r*COS(th) !*** ONLY FOR FLAT SOLID SURFACE
    END SELECT
    fregion=4.d0*WParamYL*((sigmaYL/(ueigp+SmallNumbYL))**C1YL-(sigmaYL/((ueigp+SmallNumbYL)))**C2YL+       &
    flmYL*(sigmaYL/((ueigp+SmallNumbYL)))**C3YL)
    res(m1)=res(m1)-(+LSact*wgpAYL*s1*ph(m)*(fregion*SQRT(r**2*thd**2+rd**2)))
  ENDIF

!*** CURVATURE TERMS
  IF (ATAN2(r*thd,rd).LE.0.d0) THEN
    at=ATAN2(r*thd,rd)+2.d0*Pi !*** CHANGE ATAN2 SIGN
  ELSE
    at=ATAN2(r*thd,rd)
  ENDIF
  res(m1)=res(m1)-(wgpAYL*s1*(ph(m)*thd-phs(m)*at))

!*** PRESSURE TERMS
  res(m1)=res(m1)-(-wgpAYL*s1*ph(m)*uYL(2*nnodesYL+1)*SQRT(r**2*thd**2+rd**2))

!*** ELECTRIC TERMS	
  res(m1)=res(m1)-(-wgpAYL*s1*ph(m)*bne*efieldYLGP(nell)*SQRT(r**2*thd**2+rd**2)/2.)

!*** PSEUDO-TRANSIENT CONTINUATION TERM
  SELECT CASE (TraCon)
    CASE(1)
      res(m1)=res(m1)-(+wgpAYL*s1*ph(m)*(r-ro)/TStep)
  END SELECT

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
  IF (uYL(ngl(1))*COS(uYL(ngl(3)))<HeighLimitYL) THEN
    jacob=jacob+(4.d0*(flmYL*(sigmaYL/(SmallNumbYL+ueigp))**C3YL+                                           &
    (sigmaYL/(SmallNumbYL+ueigp))**C1YL-(sigmaYL/(SmallNumbYL+ueigp))**C2YL)*                               &
    LSact*WParamYL*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*drd+                                         &
    (4.d0*(flmYL*(sigmaYL/(SmallNumbYL+ueigp))**C3YL+(sigmaYL/(SmallNumbYL+ueigp))**C1YL                    &
    -(sigmaYL/(SmallNumbYL+ueigp))**C2YL)*LSact*WParamYL*ph(m)*r**2*s1*thd*                                 &
    wgpAYL/SQRT(r**2*thd**2+rd**2))*dthd
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
          jacob=jacob+(4.d0*(flmYL*(sigmaYL/(SmallNumbYL+ueigp))**C3YL+                                     &
          (sigmaYL/(SmallNumbYL+ueigp))**C1YL-(sigmaYL/(SmallNumbYL+ueigp))**C2YL)*                         &
          LSact*WParamYL*ph(m)*r*s1*thd**2*wgpAYL/SQRT(r**2*thd**2+rd**2))*ph(n)+                           &
          (4.d0*(flmYL*(sigmaYL/(SmallNumbYL+ueigp))**C3YL+(sigmaYL/(SmallNumbYL+ueigp))**C1YL              &
          -(sigmaYL/(SmallNumbYL+ueigp))**C2YL)*LSact*WParamYL*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+         &
          rd**2))*phs(n)+(-4.d0*SQRT(r**2*thd**2+rd**2)*(C3YL*flmYL*sigmaYL*(sigmaYL/(SmallNumbYL+          &
          ueigp))**(C3YL-1)/(SmallNumbYL+ueigp)**2+C1YL*sigmaYL*(sigmaYL/(SmallNumbYL+                      &
          ueigp))**(C1YL-1)/(SmallNumbYL+ueigp)**2-C2YL*sigmaYL*(sigmaYL/(SmallNumbYL+                      &
          ueigp))**(C2YL-1)/(SmallNumbYL+ueigp)**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*dudrei(n)*ph(n)
        CASE(1)
!*** FLAT SOLID SURFACE
          jacob=jacob+(4.d0*(flmYL*(sigmaYL/(r*COS(th)+SmallNumbYL))**C3YL+                                 &  
          (sigmaYL/(r*COS(th)+SmallNumbYL))**C1YL-(sigmaYL/(r*COS(th)+SmallNumbYL))**C2YL)*                 &
          LSact*WParamYL*ph(m)*r*s1*thd**2*wgpAYL/SQRT(r**2*thd**2+rd**2)-                                  &
          4.d0*SQRT(r**2*thd**2+rd**2)*(C3YL*flmYL*sigmaYL*(sigmaYL/(r*COS(th)+                             &
          SmallNumbYL))**(C3YL-1)*COS(th)/(r*COS(th)+SmallNumbYL)**2+                                       &
          C1YL*sigmaYL*(sigmaYL/(r*COS(th)+SmallNumbYL))**(C1YL-1)*COS(th)/(r*COS(th)+                      &
          SmallNumbYL)**2-C2YL*sigmaYL*(sigmaYL/(r*COS(th)+SmallNumbYL))**(C2YL-1)*COS(th)/(r*COS(th)+      &
          SmallNumbYL)**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*ph(n)+(4.d0*(flmYL*(sigmaYL/(r*COS(th)+          &
          SmallNumbYL))**C3YL+(sigmaYL/(r*COS(th)+SmallNumbYL))**C1YL-(sigmaYL/(r*COS(th)+                  &
          SmallNumbYL))**C2YL)*LSact*WParamYL*ph(m)*rd*s1*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
      END SELECT
    ENDIF
!*** RES=PSEUDO-TRANSIENT CONTINUATION TERM
    SELECT CASE (TraCon)
      CASE(1)
        jacob=jacob+(+wgpAYL*s1*ph(m)*ph(n)/TStep)
    END SELECT
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
          jacob=jacob+(4.d0*(flmYL*(sigmaYL/(SmallNumbYL+ueigp))**C3YL+                                     &
          (sigmaYL/(SmallNumbYL+ueigp))**C1YL-(sigmaYL/(SmallNumbYL+                                        &
          ueigp))**C2YL)*LSact*WParamYL*ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)            &
          +(-4.d0*SQRT(r**2*thd**2+rd**2)*(C3YL*flmYL*sigmaYL*(sigmaYL/(SmallNumbYL+                        &
          ueigp))**(C3YL-1)/(SmallNumbYL+ueigp)**2+C1YL*sigmaYL*(sigmaYL/(SmallNumbYL+                      &
          ueigp))**(C1YL-1)/(SmallNumbYL+ueigp)**2-C2YL*sigmaYL*(sigmaYL/(SmallNumbYL+                      &
          ueigp))**(C2YL-1)/(SmallNumbYL+ueigp)**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*dudthei(n)*ph(n)
        CASE(1)
!*** FLAT SOLID SURFACE
          jacob=jacob+(4.d0*SQRT(r**2*thd**2+rd**2)*(C3YL*flmYL*r*sigmaYL*(sigmaYL/(r*COS(th)+              &
          SmallNumbYL))**(C3YL-1)*SIN(th)/(r*COS(th)+SmallNumbYL)**2+C1YL*r*sigmaYL*                        &
          (sigmaYL/(r*COS(th)+SmallNumbYL))**(C1YL-1)*SIN(th)/(r*COS(th)+SmallNumbYL)**2-                   &
          C2YL*r*sigmaYL*(sigmaYL/(r*COS(th)+SmallNumbYL))**(C2YL-1)*SIN(th)/(r*COS(th)+                    &
          SmallNumbYL)**2)*LSact*WParamYL*ph(m)*s1*wgpAYL)*ph(n)+(4.d0*(flmYL*(sigmaYL/(r*COS(th)+          &
          SmallNumbYL))**C3YL+(sigmaYL/(r*COS(th)+SmallNumbYL))**C1YL-                                      &
          (sigmaYL/(r*COS(th)+SmallNumbYL))**C2YL)*LSact*WParamYL*                                          &
          ph(m)*r**2*s1*thd*wgpAYL/SQRT(r**2*thd**2+rd**2))*phs(n)
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
