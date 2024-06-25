!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

SUBROUTINE mesh1DLIN(xfirst,xlast,nell1d,xpt,dxpt)

!*******************************************************************************
!*** DESCRIPTION: 1D EQUIDISTANT MESH GENERATOR (LINEAR ELEMENTS)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nell1d
DOUBLE PRECISION,INTENT(IN)::xfirst
DOUBLE PRECISION,INTENT(IN)::xlast
DOUBLE PRECISION,INTENT(OUT)::xpt(nell1d+1)
DOUBLE PRECISION,INTENT(OUT)::dxpt(nell1d+1)
INTEGER(KIND=8)::i
DOUBLE PRECISION::deltaX

deltaX=(xlast-xfirst)/REAL(nell1d,8)
xpt(1)=xfirst
DO i=2,nell1d+1
  xpt(i)=xpt(i-1)+deltaX
ENDDO
DO i=1,nell1d+1
  dxpt(i)=(REAL(i,8)-1.d0)/REAL(nell1d,8) !*** dxpt = d(x)/d(xlast)
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE mesh1DADA(nell1d,xpt,Tr,Tth,Tw)

!*******************************************************************************
!*** DESCRIPTION: 1D ADAPTIVE MESH GENERATOR
!*******************************************************************************

USE CommonVars
USE nodnum
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nell1d
DOUBLE PRECISION,INTENT(IN)::Tr(nell1d+1)
DOUBLE PRECISION,INTENT(IN)::Tth(nell1d+1)
DOUBLE PRECISION,INTENT(IN)::Tw(nell1d+1)
DOUBLE PRECISION,INTENT(INOUT)::xpt(nell1d+1)
INTEGER(KIND=8)::i
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::alloc
INTEGER(KIND=8)::Sam1
INTEGER(KIND=8)::Sam2
INTEGER(KIND=8)::Sxvec
INTEGER(KIND=8)::Sres
INTEGER(KIND=8)::TAMcnt
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::Tam2
DOUBLE PRECISION::normsum
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Tam1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Txvec
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::Tres

!*** STORE CURRENT MATRICES
Sam1=SIZE(am1); Sam2=SIZE(am2); Sxvec=SIZE(xvec); Sres=SIZE(res)
ALLOCATE(Tam2(Sam2),Tam1(Sam1),Txvec(Sxvec),Tres(Sres))
DO i=1,Sam2
  Tam2(i)=am2(i)
ENDDO
DO i=1,Sam1
  Tam1(i)=am1(i)
ENDDO
DO i=1,Sxvec
  Txvec(i)=xvec(i)
ENDDO
DO i=1,Sres
  Tres(i)=res(i)
ENDDO
TAMcnt=AMcnt

!*** DEALLOCATE MEMORY
DEALLOCATE(am1,am2,xvec,res)

!*** ALLOCATE MEMORY
noukADA=nell1d+1
nellADA=nell1d
nnodesADA=noukADA
nelnodesADA=2
ALLOCATE(am1(50*noukADA),am2(50*noukADA),xvec(noukADA),uADA(noukADA),          &
bcADA(noukADA),res(noukADA),ncodADA(noukADA),nopADA(nellADA,2),                &
cptADA(noukADA),dcptADA(noukADA))

!*** NODAL NUMBERING
CALL nodnumbADA

!*** COUNTERS & MATRICES INITIALIZATION
inr=0
AMcnt=0
DO i=1,SIZE(am2)
  am2(i)=0.d0 !*** INITIALIZE JACOBIAN MATRIX PATTERN
ENDDO
DO i=1,noukADA
  bcADA(i)=xpt(i)
  uADA(i)=bcADA(i) !*** INITIALIZE SOLUTION
  xvec(i)=0.d0 !*** INITIALIZE X VECTOR
  ncodADA(i)=0 !*** INITIALIZE BOUNDARY CONDITIONS FLAG
ENDDO
ncodADA(1)=1 !*** DIRICHLET BC - 1ST NODE
ncodADA(noukADA)=1 !*** DIRICHLET BC - LAST NODE

!*** EVENLY DISTRIBUTED MESH
CALL mesh1DLIN(0.d0,1.d0,nellADA,cptADA,dcptADA)

!*** FIND SPARSE JACOBIAN MATRIX PATTERN FOR ADAPTIVE MESH EQUATION
alloc=100
CALL neigElem(nnodesADA,nellADA,nelnodesADA,nopADA,kindADA,neighADA,alloc)
CALL JacobPtrn(noukADA,nelnodesADA,nopADA,kindADA,neighADA)

!*** SOLVE ADAPTIVE MESH EQUATION
DO WHILE (inr==0)!*** START NEWTON-RAPHSON ITERATIONS
  DO i=1,noukADA
    uADA(i)=uADA(i)+xvec(i)
  ENDDO
  CALL axbADA(nell1d,Tr,Tth,Tw) !*** MATRICES ASSEMBLY AND SYSTEM SOLVING
  CALL ErrorNormCalc(noukADA,xvec,normsum) !*** CALCULATE ERROR NORM
  IF (normsum<=PPres) THEN
    inr=1
  ENDIF
!  WRITE(*,'(a,E10.3)')'ADA ERROR NORM:',normsum
  IF (normsum>1.d3) THEN
    WRITE(*,'(a)')'ADA ERROR NORM>1.e3 - KILLING PROCESS...'
    STOP
  ENDIF
ENDDO
DO i=1,noukADA
  xpt(i)=uADA(i)
ENDDO

!*** CLEAR MEMORY
DEALLOCATE(am1,am2,xvec,res,nopADA,uADA,kindADA,neighADA,ncodADA,bcADA,        &
cptADA,dcptADA)

!*** RESTORE MATRICES
ALLOCATE(am1(Sam1),am2(Sam2),xvec(Sxvec),res(Sres))
DO i=1,Sam2
  am2(i)=Tam2(i)
ENDDO
DO i=1,Sam1
  am1(i)=Tam1(i)
ENDDO
DO i=1,Sxvec
  xvec(i)=Txvec(i)
ENDDO
DO i=1,Sres
  res(i)=Tres(i)
ENDDO
AMcnt=TAMcnt
DEALLOCATE(Tam1,Tam2,Txvec,Tres)

!*******************************************************************************

CONTAINS

  SUBROUTINE axbADA(nell1d,Tr,Tth,Tw)

!*******************************************************************************
!*** DESCRIPTION: BOUNDARY CONDITION IMPLEMENTATION AND 
!*** SYSTEM SOLVING FOR ADAPTIVE MESH EQUATION
!*******************************************************************************

  USE CommonVars
  IMPLICIT NONE
  INTEGER(KIND=8),INTENT(IN)::nell1d
  DOUBLE PRECISION,INTENT(IN)::Tr(nell1d+1)
  DOUBLE PRECISION,INTENT(IN)::Tth(nell1d+1)
  DOUBLE PRECISION,INTENT(IN)::Tw(nell1d+1)
  INTEGER(KIND=8)::k
  INTEGER(KIND=8)::i
  INTEGER(KIND=8)::j
  INTEGER(KIND=8)::i1
  INTEGER(KIND=8)::icnt1
  DOUBLE PRECISION::jacob

!*** MATRICES INITIALIZATION
  DO i=1,noukADA
    res(i)=0.d0
  ENDDO
  DO i=1,SIZE(am1)
    am1(i)=0.d0
  ENDDO

!*** ASSEMBLY RESIDUALS (RES) AND JACOBIAN (AM) MATRICES
  DO nell=1,nellADA
    CALL abfindADA(nell1d,Tr,Tth,Tw)
  ENDDO

!*** IMPOSE DIRICHLET BOUNDARY CONDITIONS
  DO i=1,noukADA
    IF (ncodADA(i)==1) THEN
      res(i)=bcADA(i)-uADA(i)
      jacob=0.d0
      DO j=1,kindADA(i)
        DO k=1,nelnodesADA
          i1=nopADA(neighADA(i,j),k)
          icnt1=(i-1)*noukADA+i1
          CALL buildjacob_dir(icnt1,jacob)
        ENDDO
      ENDDO
      jacob=1.d0
      icnt1=(i-1)*noukADA+i
      CALL buildjacob_dir(icnt1,jacob)
    ENDIF
  ENDDO

!*** SOLVE THE LINEAR SYSTEM USING MUMPS DIRECT SOLVER
  noukLS=noukADA
  CALL linsolver

  END SUBROUTINE

  SUBROUTINE abfindADA(nell1d,Tr,Tth,Tw)

!*******************************************************************************
!*** DESCRIPTION: RESIDUALS & JACOBIAN MATRICES ASSEMBLY FOR ADA EQUATION
!*******************************************************************************

  USE CommonVars
  IMPLICIT NONE
  INTEGER(KIND=8),INTENT(IN)::nell1d
  DOUBLE PRECISION,INTENT(IN)::Tr(nell1d+1)
  DOUBLE PRECISION,INTENT(IN)::Tth(nell1d+1)
  DOUBLE PRECISION,INTENT(IN)::Tw(nell1d+1)
  INTEGER(KIND=8)::i
  INTEGER(KIND=8)::m
  INTEGER(KIND=8)::n
  INTEGER(KIND=8)::m1
  INTEGER(KIND=8)::n1
  INTEGER(KIND=8)::icnt1
  INTEGER(KIND=8)::ngl(2)
  DOUBLE PRECISION::gpADA
  DOUBLE PRECISION::wgpADA
  DOUBLE PRECISION::u
  DOUBLE PRECISION::ud
  DOUBLE PRECISION::c
  DOUBLE PRECISION::c1
  DOUBLE PRECISION::jacob
  DOUBLE PRECISION::ADAweight
  DOUBLE PRECISION::dADAweight
  DOUBLE PRECISION::ph(2)
  DOUBLE PRECISION::phd(2)
  DOUBLE PRECISION::phs(2)

  DO i=1,nelnodesADA
    ngl(i)=nopADA(nell,i)
  ENDDO

!*** GAUSS POINTS AND WEIGHTS
  gpADA=0.5d0
  wgpADA=1.d0

!*** CALCULATE BASIS FUNCTIONS
  CALL tsfun1DLIN(gpADA,ph,phd)

!*** INITIALIZE VARIABLES
  u=0.d0; ud=0.d0; c=0.d0; c1=0.d0; ADAweight=0.d0; dADAweight=0.d0

!***ISOPARAMETRIC MAPPING
  DO n=1,nelnodesADA
    c=c+cptADA(ngl(n))*ph(n)
    c1=c1+cptADA(ngl(n))*phd(n)
    u=u+uADA(ngl(n))*ph(n)
    ADAweight=ADAweight+ABS(1.d0+Tw(ngl(n)))*ph(n)
  ENDDO
  DO n=1,nelnodesADA
    phs(n)=phd(n)/c1
    ud=ud+uADA(ngl(n))*phs(n)
    dADAweight=dADAweight+ABS(1.d0+Tw(ngl(n)))*phs(n)
  ENDDO

!*** ASSEMBLY RESIDUALS
  DO m=1,nelnodesADA
    m1=ngl(m)
    res(m1)=res(m1)-(-wgpADA*c1*ud*phs(m))
    res(m1)=res(m1)-(-wgpADA*c1*(1.d0/ADAweight*dADAweight)*ud*ph(m))
!*** ASSEMBLY JACOBIAN
      DO n=1,nelnodesADA
        n1=ngl(n)
        icnt1=(m1-1)*noukADA+n1
        jacob=+(-wgpADA*c1*phs(m)*phs(n))
        jacob=jacob+(-wgpADA*c1*(1.d0/ADAweight*dADAweight)*ph(m)*phs(n))
        CALL buildjacob_nor(icnt1,jacob)
      ENDDO
  ENDDO

  END SUBROUTINE

!*******************************************************************************

END SUBROUTINE

!*******************************************************************************

SUBROUTINE geogenEL

!*******************************************************************************
!*** DESCRIPTION: GEOMETRY GENERATOR (PRODUCE INPUT FILE FOR THE MESH GENERATOR)
!*** (EL EQUATIONS)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::nnodessolid
INTEGER(KIND=8)::i
DOUBLE PRECISION::PillarFun
DOUBLE PRECISION::pointx(nnodesYL)
DOUBLE PRECISION::pointy(nnodesYL)
DOUBLE PRECISION::deltaxsolid
DOUBLE PRECISION,ALLOCATABLE::xsolid(:)
DOUBLE PRECISION,ALLOCATABLE::ysolid(:)
DOUBLE PRECISION,ALLOCATABLE::xsoliddown(:)
DOUBLE PRECISION,ALLOCATABLE::ysoliddown(:)

!*** LOAD DROPLET SURFACE COORDINATES
DO i=1,nnodesYL
  pointx(i)=uYL(i)*SIN(uYL(nnodesYL+i))
  pointy(i)=uYL(i)*COS(uYL(nnodesYL+i))
ENDDO

!*** CHECK DOMAIN BOUNDARIES
IF (ylimitEL<=pointy(1).OR.xlimitEL<=MAXVAL(pointx)) THEN
  WRITE(*,'(a)')"ERROR AT coord/geogenEL: EXPAND 2D DOMAIN"
  STOP
!*** DYNAMIC DOMAIN MODIFICATION
!  ylimitEL=pointy(1)+0.5d0
!  xlimitEL=MAXVAL(pointx)+0.5d0
ENDIF

!*** EVALUATE/LOAD DIELECTRIC GEOMETRY
nnodessolid=500 !*** NUMBER OF NODES FOR SOLID DIELECTRIC DISCRETISATION
deltaxsolid=(2.d0/3.d0)*xlimitEL/REAL(nnodessolid,8)
ALLOCATE(xsolid(nnodessolid),ysolid(nnodessolid))
DO i=1,nnodessolid
  xsolid(i)=REAL(i-1,8)*deltaxsolid
  ysolid(i)=PillarFun(xsolid(i)) !*** SOLID SURFACE STRUCTURE EQUATION
ENDDO
ALLOCATE(xsoliddown(nnodessolid),ysoliddown(nnodessolid))
SELECT CASE (flatSolidDown)
  CASE(1)
    DO i=1,nnodessolid
      xsoliddown(i)=REAL(i-1,8)*deltaxsolid
      ysoliddown(i)=MINVAL(ysolid)-DielThickEL
    ENDDO
  CASE(0)
    OPEN(101,FILE='diel_geom.txt')
    DO i=1,nnodessolid
      READ(101,*)xsoliddown(i),ysoliddown(i)
    ENDDO
    CLOSE(101)
END SELECT

!*** MESH GENERATOR SOFTWARE
SELECT CASE (meshGen)

  CASE(1) !*** TRIANGLE MESH GENERATOR - (.POLY FILE)
    OPEN(10,FILE='dropletgeom.poly')
    WRITE(10,611)nnodesYL+5+2*nnodessolid !*** NUMBER OF NODES
!*** DROPLET GEOMETRY
    WRITE(10,601)1,pointx(1),ylimitEL
    DO i=2,nnodesYL+1
      WRITE(10,601)i,pointx(i-1),pointy(i-1)
    ENDDO
!*** DIELECTRIC GEOMETRY
    WRITE(10,601)nnodesYL+2,xsolid(1),ysolid(1)
    DO i=2,nnodessolid
      WRITE(10,601)nnodesYL+2+i-1,xsolid(i),ysolid(i)
    ENDDO
    WRITE(10,601)nnodesYL+2+nnodessolid,xlimitEL,ysolid(nnodessolid)
    WRITE(10,601)nnodesYL+3+nnodessolid,xlimitEL,ylimitEL
    WRITE(10,601)nnodesYL+4+nnodessolid,xsoliddown(1),ysoliddown(1)
    DO i=2,nnodessolid
!*** DIELECTRIC BASE
      WRITE(10,601)nnodesYL+4+nnodessolid+i-1,xsoliddown(i),ysoliddown(i)
    ENDDO
    WRITE(10,601)nnodesYL+4+2*nnodessolid,xlimitEL,ysolid(nnodessolid)-        &
    DielThickEL
    WRITE(10,601)nnodesYL+5+2*nnodessolid,xsolid(nnodessolid),ylimitEL
    WRITE(10,612)nnodesYL+6+2*nnodessolid !*** NUMBER OF LINES
    DropLinefl=nnodesYL+15
    WRITE(10,613)1,1,2,1
    DO i=2,nnodesYL
      WRITE(10,613)i,i,i+1,DropLinefl
    ENDDO
    DO i=nnodesYL+1,nnodesYL+2+nnodessolid
      WRITE(10,613)i,i,i+1,1
    ENDDO
    WRITE(10,613)nnodesYL+3+nnodessolid,nnodesYL+3+nnodessolid,1,1
    WRITE(10,613)nnodesYL+4+nnodessolid,nnodesYL+2,nnodesYL+4+nnodessolid,1
    DielLinefl=nnodesYL+17
    WRITE(10,613)nnodesYL+5+nnodessolid,nnodesYL+4+nnodessolid,nnodesYL+4+     &
    nnodessolid+1,DielLinefl
    DO i=1,nnodessolid-1
      WRITE(10,613)nnodesYL+5+nnodessolid+i,nnodesYL+4+nnodessolid+i,nnodesYL+ &
      4+nnodessolid+i+1,DielLinefl
    ENDDO
    WRITE(10,613)nnodesYL+5+2*nnodessolid,nnodesYL+4+2*nnodessolid,nnodesYL+   &
    2+nnodessolid,DielLinefl
	WRITE(10,613)nnodesYL+6+2*nnodessolid,nnodesYL+5+2*nnodessolid,        &
	nnodesYL+1+nnodessolid
	WRITE(10,'(a)')'0' !*** NUMBER OF HOLES
	WRITE(10,'(a)')'3' !*** NUMBER OF DOMAINS
	AirSurffl=nnodesYL+13
    DielSurffl=nnodesYL+14
!*** AIR DOMAIN (LOW MESH DENSITY)
    WRITE(10,614)1,xlimitEL-2.d-3,2.d-3,AirSurffl,ldensEL
!*** DIELECTRIC DOMAIN
    WRITE(10,614)2,2.d-3,-2.d-3,DielSurffl,hdensEL
!*** AIR DOMAIN (HIGH MESH DENSITY)
    WRITE(10,614)3,2.d-3,ylimitEL-2.d-3,AirSurffl,hdensEL
    CLOSE(10)

  CASE(2) !*** GMSH MESH GENERATOR (.GEO FILE)
    OPEN(10,FILE='dropletgeom.geo')
!*** GEOMETRY ENTITIES
!*** DROPLET GEOMETRY
    WRITE(10,602)1,pointx(1),ylimitEL,ldensEL
    DO i=2,nnodesYL+1
      WRITE(10,602)i,pointx(i-1),pointy(i-1),hdensEL
    ENDDO
!*** DIELECTRIC GEOMETRY
    WRITE(10,602)nnodesYL+2,xsolid(1),ysolid(1),hdensEL
    DO i=2,nnodessolid
      WRITE(10,602)nnodesYL+2+i-1,xsolid(i),ysolid(i),hdensEL
    ENDDO
    WRITE(10,602)nnodesYL+2+nnodessolid,xlimitEL,ysolid(nnodessolid),ldensEL
    WRITE(10,602)nnodesYL+3+nnodessolid,xlimitEL,ylimitEL,ldensEL
    DO i=1,nnodesYL+2+nnodessolid
      WRITE(10,603)i,i,i+1
    ENDDO
    WRITE(10,603)nnodesYL+3+nnodessolid,nnodesYL+3+nnodessolid,1
    WRITE(10,604)nnodesYL+4+nnodessolid
    DO i=1,nnodesYL+2+nnodessolid
      WRITE(10,605)i
    ENDDO
    WRITE(10,606)nnodesYL+3+nnodessolid
    WRITE(10,607)nnodesYL+5+nnodessolid,nnodesYL+4+nnodessolid
    WRITE(10,602)nnodesYL+4+nnodessolid,xsoliddown(1),ysoliddown(1),hdensEL*   &
    GrateEL
    DO i=2,nnodessolid
      WRITE(10,602)nnodesYL+4+nnodessolid+i-1,xsoliddown(i),ysoliddown(i),     &
      hdensEL*GrateEL !*** DIELECTRIC BASE
    ENDDO
    WRITE(10,602)nnodesYL+4+2*nnodessolid,xlimitEL,ysolid(nnodessolid)-        &
    DielThickEL,ldensEL
    WRITE(10,603)nnodesYL+4+nnodessolid,nnodesYL+2,nnodesYL+4+nnodessolid
    WRITE(10,603)nnodesYL+5+nnodessolid,nnodesYL+4+nnodessolid,nnodesYL+4+     &
    nnodessolid+1
    DO i=1,nnodessolid-1
      WRITE(10,603)nnodesYL+5+nnodessolid+i,nnodesYL+4+nnodessolid+i,nnodesYL+ &
      4+nnodessolid+i+1
    ENDDO
    WRITE(10,603)nnodesYL+5+2*nnodessolid,nnodesYL+4+2*nnodessolid,nnodesYL+   &
    2+nnodessolid
    WRITE(10,604)nnodesYL+5+nnodessolid
    WRITE(10,605)nnodesYL+5+nnodessolid
    DO i=1,nnodessolid
      WRITE(10,605)nnodesYL+5+nnodessolid+i
    ENDDO
    DO i=nnodessolid+nnodesYL+1,nnodesYL+2,-1
      WRITE(10,605)-i
    ENDDO
    WRITE(10,606)nnodesYL+4+nnodessolid
    WRITE(10,607)nnodesYL+6+nnodessolid,nnodesYL+5+nnodessolid
    DEALLOCATE (xsolid,ysolid,xsoliddown,ysoliddown) !*** CLEAR MEMORY
!*** PHYSICAL ENTITIES
!*** PSYSICAL POINTS
    WRITE(10,608)nnodesYL+12
    dropPointsfl=nnodesYL+12
    DO i=2,nnodesYL
      WRITE(10,605)i
    ENDDO
    WRITE(10,606)nnodesYL+1
    WRITE(10,608)nnodesYL+16
    WRITE(10,605)nnodesYL+4+nnodessolid
    DO i=1,nnodessolid-1
      WRITE(10,605)nnodesYL+4+nnodessolid+i
    ENDDO
    WRITE(10,606)nnodesYL+4+2*nnodessolid
!*** PSYSICAL SURFACES
!*** SURROUNDING MEDIUM DOMAIN
    WRITE(10,609)nnodesYL+13,nnodesYL+5+nnodessolid
    AirSurffl=nnodesYL+13
!*** DIELECTRIC DOMAIN
    WRITE(10,609)nnodesYL+14,nnodesYL+6+nnodessolid
    DielSurffl=nnodesYL+14
!*** PSYSICAL LINES (FOR BOUNDARY CONDITIONS)
    WRITE(10,610)nnodesYL+15 !*** DROPLET SURFACE
    DropLinefl=nnodesYL+15
    DO i=2,nnodesYL-1
      WRITE(10,605)i
    ENDDO
    WRITE(10,606)nnodesYL
    WRITE(10,610)nnodesYL+17 !*** ELECTRODE SURFACE
    DielLinefl=nnodesYL+17
    WRITE(10,605)nnodesYL+5+nnodessolid
    DO i=1,nnodessolid-1
      WRITE(10,605)nnodesYL+5+nnodessolid+i
    ENDDO
    WRITE(10,606)nnodesYL+5+2*nnodessolid
    CLOSE(10)

END SELECT

601 FORMAT(i6,1x,f14.11,1x,f14.11)
602 FORMAT('Point(',i6,')={',f14.11,',',f14.11,',0,',f14.11'};')
603 FORMAT('Line(',i6,')={',i6,',',i6,'}',';')
604 FORMAT('Line Loop(',i6,')={')
605 FORMAT(i6,',')
606 FORMAT(i6,'};')
607 FORMAT('Plane Surface(',i6,') = {',i6,'};')
608 FORMAT('Physical Point(',i6,')={')
609 FORMAT('Physical Surface(',i6,')={',i6,'};')
610 FORMAT('Physical Line(',i6,')={')
611 FORMAT(i6,1x,'2',1x,'0',1x,'0')
612 FORMAT(i6,1x,'1')
613 FORMAT(i6,1x,i6,1x,i6,1x,i6)
614 FORMAT(i1,1x,f14.11,1x,f14.11,1x,i6,1x,f14.11)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE geogenEI

!*******************************************************************************
!*** DESCRIPTION: GEOMETRY GENERATOR (PRODUCE INPUT FILE FOR THE MESH GENERATOR)
!*** (EI EQUATION)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::nnodespillar
INTEGER(KIND=8)::nnodessolid
INTEGER(KIND=8)::i
DOUBLE PRECISION::PillarFun
DOUBLE PRECISION::deltaxsolid
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::xsolid
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::ysolid

!*** EVALUATE SOLID SURFACE GEOMETRY
nnodespillar=500 !*** NUMBER OF NODES FOR SOLID DIELECTRIC DISCRETISATION
IF (geommirEI==1) THEN
  nnodessolid=nnodespillar
ELSE IF (geommirEI==0) THEN
  nnodessolid=nnodespillar*FLOOR(xlimitEI/(wpillar+dpillar))
ENDIF
deltaxsolid=xlimitEI/REAL(nnodessolid,8)
ALLOCATE(xsolid(nnodessolid),ysolid(nnodessolid))
DO i=1,nnodessolid
  xsolid(i)=REAL(i,8)*deltaxsolid
  ysolid(i)=PillarFun(xsolid(i))
ENDDO

!*** MESH GENERATOR SOFTWARE
SELECT CASE (meshGen)

  CASE(1) !*** TRIANGLE MESH GENERATOR (.POLY FILE)
    OPEN(10,FILE='eikongeom.poly')
    WRITE(10,701)nnodessolid+3 !*** NUMBER OF NODES
    WRITE(10,702)1,xfirstEI,0.d0
    DO i=1,nnodessolid
      WRITE(10,702)i+1,xsolid(i),ysolid(i)
    ENDDO
    WRITE(10,702)nnodessolid+2,xlimitEI,ylimitEI
    WRITE(10,702)nnodessolid+3,xfirstEI,ylimitEI
    WRITE(10,703)nnodessolid+3 !*** NUMBER OF LINES
    SolidLinefl=nnodessolid+3
    UpperLinefl=nnodessolid+4
    DO i=1,nnodessolid
      WRITE(10,704)i,i,i+1,SolidLinefl
    ENDDO
    WRITE(10,704)nnodessolid+1,nnodessolid+1,nnodessolid+2,1
    WRITE(10,704)nnodessolid+2,nnodessolid+2,nnodessolid+3,UpperLinefl
    WRITE(10,704)nnodessolid+3,nnodessolid+3,1,1
    WRITE(10,'(a)')'0' !*** NUMBER OF HOLES
    WRITE(10,'(a)')'1' !*** NUMBER OF DOMAINS
    WRITE(10,705)1,xlimitEI/2.d0,ylimitEI/2.d0,1,hdensEI
    CLOSE(10)

  CASE(2) !*** GMSH MESH GENERATOR (.GEO FILE)
    OPEN(10,FILE='eikongeom.geo')
    WRITE(10,706)1,xfirstEI,0.d0,hdensEI
    DO i=1,nnodessolid
      WRITE(10,706)i+1,xsolid(i),ysolid(i),hdensEI
    ENDDO
    WRITE(10,706)nnodessolid+2,xlimitEI,ylimitEI,ldensEI
    WRITE(10,706)nnodessolid+3,xfirstEI,ylimitEI,ldensEI
    DO i=1,nnodessolid
      WRITE(10,707)i,i,i+1
    ENDDO
    WRITE(10,707)nnodessolid+1,nnodessolid+1,nnodessolid+2
    WRITE(10,707)nnodessolid+2,nnodessolid+2,nnodessolid+3
    WRITE(10,707)nnodessolid+3,nnodessolid+3,1
    WRITE(10,708)nnodessolid+4
    DO i=1,nnodessolid+2
      WRITE(10,709)i
    ENDDO
    WRITE(10,710)nnodessolid+3
    WRITE(10,711)nnodessolid+5,nnodessolid+4
    WRITE(10,712)nnodessolid+6,nnodessolid+5
    AirSurffl=nnodessolid+6
    WRITE(10,713)nnodessolid+7
    DO i=1,nnodessolid-1
      WRITE(10,709)i
    ENDDO
    WRITE(10,710)nnodessolid
    SolidLinefl=nnodessolid+7
    WRITE(10,713)nnodessolid+8
    WRITE(10,710)nnodessolid+2
    UpperLinefl=nnodessolid+8
    CLOSE(10)

END SELECT

701 FORMAT(i6,1x,'2',1x,'0',1x,'0')
702 FORMAT(i6,1x,f14.11,1x,f14.11)
703 FORMAT(i6,1x,'1')
704 FORMAT(i6,1x,i6,1x,i6,1x,i6)
705 FORMAT(i1,1x,f14.11,1x,f14.11,1x,i1,1x,f14.11)
706 FORMAT('Point(',i6,')={',f14.11,',',f14.11,',0,',f14.11'};')
707 FORMAT('Line(',i6,')={',i6,',',i6,'}',';')
708 FORMAT('Line Loop(',i6,')={')
709 FORMAT(i6,',')
710 FORMAT(i6,'};')
711 FORMAT('Plane Surface(',i6,') = {',i6,'};')
712 FORMAT('Physical Surface(',i6,')={',i6,'};')
713 FORMAT('Physical Line(',i6,')={')

END SUBROUTINE

!*******************************************************************************

SUBROUTINE readmeshEL

!*******************************************************************************
!*** DESCRIPTION: IMPORT 2D MESH (EL EQUATIONS)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::null
INTEGER(KIND=8)::i
INTEGER(KIND=8)::ndnum
INTEGER(KIND=8)::nellcnt
INTEGER(KIND=8)::YLncnt
DOUBLE PRECISION::zcordi
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::elmtype
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::physent
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::geoment
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:)::totnopEL

SELECT CASE (meshGen) !*** MESH GENERATOR

  CASE(1) !***  TRIANGLE MESH GENERATOR (.NODE & .ELE FILES)
    OPEN(10,FILE='dropletgeom.1.node')
	READ(10,*)nnodesEL,null,null,null
	ALLOCATE(ncodEL(nnodesEL),bcEL(nnodesEL),xptEL(nnodesEL),              &
	yptEL(nnodesEL),physent(nnodesEL))
    DO i=1,nnodesEL
      READ(10,*)null,xptEL(i),yptEL(i),physent(i)
    ENDDO
    CLOSE(10)
    OPEN(10,FILE='dropletgeom.1.ele')
    READ(10,*)nellEL,null,null
!*** READ NOP ARRAY
	ALLOCATE(nopEL(nellEL,6),nellidEL(nellEL))
    SELECT CASE (eltypeEL)
      CASE(1)
        DO i=1,nellEL
          READ(10,*)null,nopEL(i,1),nopEL(i,2),nopEL(i,3),nellidEL(i)
        ENDDO
      CASE(2)
        DO i=1,nellEL
          READ(10,*)null,nopEL(i,1),nopEL(i,2),nopEL(i,3),nopEL(i,5),          &
          nopEL(i,6),nopEL(i,4),nellidEL(i)
        ENDDO
    END SELECT
    CLOSE(10)
!*** IMPOSE BOUNDARY CONDITIONS
    DO i=1,nnodesEL
      ncodEL(i)=0
      bcEL(i)=0.d0
    ENDDO
    YLncnt=0
    DO i=2,nnodesYL+1
      YLncnt=YLncnt+1
      YLnodes(YLncnt)=i
    ENDDO
    DO i=1,nnodesEL
      IF (physent(i)==DropLinefl) THEN
        ncodEL(i)=1
        bcEL(i)=1.d0
      ELSE IF (physent(i)==DielLinefl) THEN
        ncodEL(i)=1
        bcEL(i)=0.d0
      ENDIF
    ENDDO
    noukEL=nnodesEL
    DEALLOCATE(physent)

  CASE(2) !*** GMSH MESH GENERATOR (.MSH FILE)
    OPEN(10,FILE='dropletgeom.msh')
    DO i=1,4
      READ(10,*)
    ENDDO
    READ(10,*)nnodesEL !*** READ NUMBER OF NODES
    ALLOCATE(xptEL(nnodesEL),yptEL(nnodesEL))
    DO i=1,nnodesEL
      READ(10,*)ndnum,xptEL(i),yptEL(i),zcordi !*** READ NODE COORDINATES
    ENDDO
    DO i=1,2
      READ(10,*)
    ENDDO
    READ(10,*)nellEL !*** READ NUMBER OF ELEMENTES
    ALLOCATE(totnopEL(nellEL,6),ncodEL(nnodesEL),bcEL(nnodesEL),               &
    elmtype(nellEL),physent(nellEL),geoment(nellEL),nellidEL(nellEL))
    DO i=1,nellEL
!*** READ PHYSICAL ENTITIES
      READ(10,*)null,elmtype(i),null,physent(i),geoment(i)
    ENDDO
    REWIND(10)
    DO i=1,5
      READ(10,*)
    ENDDO
    DO i=1,nnodesEL
      READ(10,*)
    ENDDO
    DO i=1,3
      READ(10,*)
    ENDDO
    nellcnt=0
!*** READ NOP ARRAY
    SELECT CASE (eltypeEL)
      CASE(1)
        DO i=1,nellEL
          SELECT CASE(elmtype(i))
            CASE(15)
              READ(10,*)null,null,null,null,null,totnopEL(i,1)
            CASE(1)
              READ(10,*)null,null,null,null,null,totnopEL(i,1),totnopEL(i,2)
            CASE(2)
              READ(10,*)null,null,null,null,null,totnopEL(i,1),totnopEL(i,2),  &
              totnopEL(i,3)
              nellcnt=nellcnt+1
          END SELECT
        ENDDO
      CASE(2)
        DO i=1,nellEL
          SELECT CASE(elmtype(i))
            CASE(15)
              READ(10,*)null,null,null,null,null,totnopEL(i,1)
            CASE(8)
              READ(10,*)null,null,null,null,null,totnopEL(i,1),totnopEL(i,2),  &
              totnopEL(i,3)
            CASE(9)
              READ(10,*)null,null,null,null,null,totnopEL(i,1),totnopEL(i,2),  &
              totnopEL(i,3),totnopEL(i,4),totnopEL(i,5),totnopEL(i,6)
              nellcnt=nellcnt+1
          END SELECT
        ENDDO
      END SELECT
    CLOSE(10)
    ALLOCATE(nopEL(nellcnt,6))
    nellcnt=0
    SELECT CASE (eltypeEL)
      CASE(1)
        DO i=1,nellEL
          IF (elmtype(i)==2) THEN
            nellcnt=nellcnt+1
            nopEL(nellcnt,1)=totnopEL(i,1)
            nopEL(nellcnt,2)=totnopEL(i,2)
            nopEL(nellcnt,3)=totnopEL(i,3)
            nellidEL(nellcnt)=physent(i)
          ENDIF
        ENDDO
      CASE(2)
        DO i=1,nellEL
          IF (elmtype(i)==9) THEN
            nellcnt=nellcnt+1
            nopEL(nellcnt,1)=totnopEL(i,1)
            nopEL(nellcnt,2)=totnopEL(i,2)
            nopEL(nellcnt,3)=totnopEL(i,3)
            nopEL(nellcnt,4)=totnopEL(i,4)
            nopEL(nellcnt,5)=totnopEL(i,5)
            nopEL(nellcnt,6)=totnopEL(i,6)
            nellidEL(nellcnt)=physent(i)
          ENDIF
        ENDDO
    END SELECT
!*** IMPOSE BOUNDARY CONDITIONS
    DO i=1,nnodesEL
      ncodEL(i)=0
      bcEL(i)=0.d0
    ENDDO
    YLncnt=0
    SELECT CASE (eltypeEL)
      CASE(1)
        DO i=1,nellEL
          IF (elmtype(i)==1.AND.physent(i)==DropLinefl) THEN
            ncodEL(totnopEL(i,1))=1; ncodEL(totnopEL(i,2))=1
            bcEL(totnopEL(i,1))=1.d0; bcEL(totnopEL(i,2))=1.d0
          ELSE IF (elmtype(i)==1.AND.physent(i)==DielLinefl) THEN
            ncodEL(totnopEL(i,1))=1; ncodEL(totnopEL(i,2))=1
            bcEL(totnopEL(i,1))=0.d0; bcEL(totnopEL(i,2))=0.d0
          ENDIF
          IF (elmtype(i)==15.AND.physent(i)==dropPointsfl) THEN
            YLncnt=YLncnt+1
            YLnodes(YLncnt)=totnopEL(i,1)
          ENDIF
        ENDDO
      CASE(2)
        DO i=1,nellEL
          IF (elmtype(i)==8.AND.physent(i)==DropLinefl) THEN
            ncodEL(totnopEL(i,1))=1
            ncodEL(totnopEL(i,2))=1
            ncodEL(totnopEL(i,3))=1
            bcEL(totnopEL(i,1))=1.d0
            bcEL(totnopEL(i,2))=1.d0
            bcEL(totnopEL(i,3))=1.d0
          ELSE IF (elmtype(i)==8.AND.physent(i)==DielLinefl) THEN
            ncodEL(totnopEL(i,1))=1
            ncodEL(totnopEL(i,2))=1
            ncodEL(totnopEL(i,3))=1
            bcEL(totnopEL(i,1))=0.d0
            bcEL(totnopEL(i,2))=0.d0
            bcEL(totnopEL(i,3))=0.d0
          ENDIF
          IF (elmtype(i)==15.AND.physent(i)==dropPointsfl) THEN
            YLncnt=YLncnt+1
            YLnodes(YLncnt)=totnopEL(i,1)
          ENDIF
        ENDDO
    END SELECT
    nellEL=nellcnt
    noukEL=nnodesEL
    DEALLOCATE(totnopEL,elmtype,physent,geoment)

END SELECT

END SUBROUTINE

!*******************************************************************************

SUBROUTINE readmeshEI

!*******************************************************************************
!*** DESCRIPTION: IMPORT 2D MESH (EI EQUATION)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::null
INTEGER(KIND=8)::i
INTEGER(KIND=8)::ndnum
INTEGER(KIND=8)::nellcnt
INTEGER(KIND=8)::YLncnt
DOUBLE PRECISION::zcordi
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::elmtype
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::physent
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::geoment
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:)::totnopEI

SELECT CASE (meshGen) !*** MESH GENERATOR

  CASE(1) !***  TRIANGLE MESH GENERATOR (.NODE & .ELE FILES)
    OPEN(10,FILE='eikongeom.1.node')
    READ(10,*)nnodesEI,null,null,null
    ALLOCATE(ncodEI(nnodesEI),bcEI(nnodesEI),xptEI(nnodesEI),yptEI(nnodesEI),  &
    physent(nnodesEI))
    DO i=1,nnodesEI
      READ(10,*)null,xptEI(i),yptEI(i),physent(i)
    ENDDO
    CLOSE(10)
    OPEN(10,FILE='eikongeom.1.ele')
    READ(10,*)nellEI,null,null
!*** READ NOP ARRAY
    ALLOCATE(nopEI(nellEI,6))
    SELECT CASE (eltypeEI)
      CASE(1)
        DO i=1,nellEI
          READ(10,*)null,nopEI(i,1),nopEI(i,2),nopEI(i,3),null
        ENDDO
      CASE(2)
        DO i=1,nellEI
          READ(10,*)null,nopEI(i,1),nopEI(i,2),nopEI(i,3),nopEI(i,5),          &
          nopEI(i,6),nopEI(i,4),null
        ENDDO
    END SELECT
    CLOSE(10)
!*** IMPOSE BOUNDARY CONDITIONS
    DO i=1,nnodesEI
      ncodEI(i)=0
      bcEI(i)=1.d0
    ENDDO
    DO i=1,nnodesEI
      IF (physent(i)==SolidLinefl) THEN
        ncodEI(i)=1
        bcEI(i)=0.d0
      ELSE IF (physent(i)==UpperLinefl) THEN
        ncodEI(i)=1
        bcEI(i)=1.d0
      ENDIF
    ENDDO
    noukEI=nnodesEI
    DEALLOCATE(physent)

  CASE(2) !*** GMSH MESH GENERATOR (.MSH FILE)
    OPEN(10,FILE='eikongeom.msh')
    DO i=1,4
      READ(10,*)
    ENDDO
    READ(10,*)nnodesEI !*** READ NUMBER OF NODES
    ALLOCATE(xptEI(nnodesEI),yptEI(nnodesEI))
    DO i=1,nnodesEI
      READ(10,*)ndnum,xptEI(i),yptEI(i),zcordi !*** READ NODE COORDINATES
    ENDDO
    DO i=1,2
      READ(10,*)
    ENDDO
    READ(10,*)nellEI !*** READ NUMBER OF ELEMENTES
    ALLOCATE(totnopEI(nellEI,6),ncodEI(nnodesEI),bcEI(nnodesEI),               &
    elmtype(nellEI),physent(nellEI),geoment(nellEI))
    DO i=1,nellEI
!*** READ PHYSICAL ENTITIES
      READ(10,*)null,elmtype(i),null,physent(i),geoment(i)
    ENDDO
    REWIND(10)
    DO i=1,5
      READ(10,*)
    ENDDO
    DO i=1,nnodesEI
      READ(10,*)
    ENDDO
    DO i=1,3
      READ(10,*)
    ENDDO
    nellcnt=0
!*** READ NOP ARRAY
    SELECT CASE (eltypeEI)
      CASE(1)
        DO i=1,nellEI
          SELECT CASE(elmtype(i))
            CASE(15)
              READ(10,*)null,null,null,null,null,totnopEI(i,1)
            CASE(1)
              READ(10,*)null,null,null,null,null,totnopEI(i,1),totnopEI(i,2)
            CASE(2)
              READ(10,*)null,null,null,null,null,totnopEI(i,1),totnopEI(i,2),  &
              totnopEI(i,3)
              nellcnt=nellcnt+1
          END SELECT
        ENDDO
      CASE(2)
        DO i=1,nellEI
          SELECT CASE(elmtype(i))
            CASE(15)
              READ(10,*)null,null,null,null,null,totnopEI(i,1)
            CASE(8)
              READ(10,*)null,null,null,null,null,totnopEI(i,1),totnopEI(i,2),  &
              totnopEI(i,3)
            CASE(9)
              READ(10,*)null,null,null,null,null,totnopEI(i,1),totnopEI(i,2),  &
              totnopEI(i,3),totnopEI(i,4),totnopEI(i,5),totnopEI(i,6)
              nellcnt=nellcnt+1
          END SELECT
        ENDDO
      END SELECT
    CLOSE(10)
    ALLOCATE(nopEI(nellcnt,6))
    nellcnt=0
    SELECT CASE (eltypeEI)
      CASE(1)
        DO i=1,nellEI
          IF (elmtype(i)==2) THEN
            nellcnt=nellcnt+1
            nopEI(nellcnt,1)=totnopEI(i,1)
            nopEI(nellcnt,2)=totnopEI(i,2)
            nopEI(nellcnt,3)=totnopEI(i,3)
          ENDIF
        ENDDO
      CASE(2)
        DO i=1,nellEI
          IF (elmtype(i)==9) THEN
            nellcnt=nellcnt+1
            nopEI(nellcnt,1)=totnopEI(i,1)
            nopEI(nellcnt,2)=totnopEI(i,2)
            nopEI(nellcnt,3)=totnopEI(i,3)
            nopEI(nellcnt,4)=totnopEI(i,4)
            nopEI(nellcnt,5)=totnopEI(i,5)
            nopEI(nellcnt,6)=totnopEI(i,6)
          ENDIF
        ENDDO
    END SELECT
!*** IMPOSE BOUNDARY CONDITIONS
    DO i=1,nnodesEI
      ncodEI(i)=0
      bcEI(i)=0.d0
    ENDDO
    YLncnt=0
    SELECT CASE (eltypeEI)
      CASE(1)
        DO i=1,nellEI
          IF (elmtype(i)==1.AND.physent(i)==SolidLinefl) THEN
            ncodEI(totnopEI(i,1))=1
            ncodEI(totnopEI(i,2))=1
            bcEI(totnopEI(i,1))=0.d0
            bcEI(totnopEI(i,2))=0.d0
          ELSE IF (elmtype(i)==1.AND.physent(i)==UpperLinefl) THEN
            ncodEI(totnopEI(i,1))=1
            ncodEI(totnopEI(i,2))=1
            bcEI(totnopEI(i,1))=1.d0
            bcEI(totnopEI(i,2))=1.d0
          ENDIF
        ENDDO
      CASE(2)
        DO i=1,nellEI
          IF (elmtype(i)==8.AND.physent(i)==SolidLinefl) THEN
            ncodEI(totnopEI(i,1))=1
            ncodEI(totnopEI(i,2))=1
            ncodEI(totnopEI(i,3))=1
            bcEI(totnopEI(i,1))=0.d0
            bcEI(totnopEI(i,2))=0.d0
            bcEI(totnopEI(i,3))=0.d0
          ELSE IF (elmtype(i)==8.AND.physent(i)==UpperLinefl) THEN
            ncodEI(totnopEI(i,1))=1
            ncodEI(totnopEI(i,2))=1
            ncodEI(totnopEI(i,3))=1
            bcEI(totnopEI(i,1))=1.d0
            bcEI(totnopEI(i,2))=1.d0
            bcEI(totnopEI(i,3))=1.d0
          ENDIF
        ENDDO
    END SELECT
    nellEI=nellcnt
    noukEI=nnodesEI
    DEALLOCATE(totnopEI,elmtype,physent,geoment)

END SELECT

END SUBROUTINE

!*******************************************************************************

DOUBLE PRECISION FUNCTION PillarFun(x)

!*******************************************************************************
!*** DESCRIPTION: SOLID SURFACE GEOMETRY FUNCTION
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
INTEGER(KIND=8)::ip
DOUBLE PRECISION::rend
DOUBLE PRECISION::rstart

!*** NUMBER OF PILLARS
SELECT CASE (SolverID)
  CASE (1)
    IF (geommirEI==1) THEN
      nopillars=1
    ELSE IF (geommirEI==0) THEN
      nopillars=FLOOR(xlimitEI/(wpillar+dpillar))
    ENDIF
  CASE DEFAULT
    nopillars=FLOOR((2.d0/3.d0)*xlimitEL/(wpillar+dpillar))
END SELECT
!*** MULTIPILLAR FUNCTION
SELECT CASE (flatSolidUp)
  CASE(1)
    PillarFun=0.d0
  CASE(0)
     PillarFun=0.d0
     DO ip=1,nopillars
       rstart=(REAL(ip)-1.d0)*(wpillar+dpillar)+wpillar*0.5d0
       rend=rstart+dpillar
       PillarFun=PillarFun-0.5d0*hpillar*(TANH((x-rstart)/smpillar)-           &
       TANH((x-rend)/smpillar))
     ENDDO
     IF (x>((wpillar+dpillar)*nopillars-dpillar)) THEN
       PillarFun=-hpillar
     ENDIF
!     PillarFun=(1.d0-erf((x-1.3d0)*10.d0))/5.d0-0.4d0 !*** SINGLE PILLAR
END SELECT

END FUNCTION

!*******************************************************************************
