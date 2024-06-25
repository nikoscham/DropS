!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

SUBROUTINE EISolver

!*******************************************************************************
!*** DESCRIPTION: EIKONAL EQUATION SOLVER
!*******************************************************************************

USE CommonVars
USE nodnum
USE defpar
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::inr
INTEGER(KIND=8)::alloc
INTEGER(KIND=8)::eikcnt
DOUBLE PRECISION::normsum
CHARACTER(LEN=1)::ElemOrder

!*** SOLVER INITIALIZATION
CALL CHDIR(ADJUSTL(Dir))
CALL defparamYL !*** READ AYL PARAMETERS
CALL defparamEI !*** READ EI PARAMETERS
CALL CHDIR(ADJUSTL(ResDir))

!*** GENERATE SOLID DIELECTRIC GEOMETRY
CALL geogenEI

!*** 2D MESHING
SELECT CASE (meshGen)
  CASE(1) !*** GENERATE UNSTRUCTURED MESH USING TRIANGLE
    SELECT CASE (eltypeEI)
      CASE(1)
        CALL SYSTEM("triangle -A -a -q25 -Q eikongeom.poly")
      CASE(2)
        CALL SYSTEM("triangle -A -a -q25 -Q -o2 eikongeom.poly")
    END SELECT
  CASE(2) !*** GENERATE UNSTRUCTURED MESH USING GMSH
    WRITE(ElemOrder,'(i1)')eltypeEI
    CALL SYSTEM("gmsh eikongeom.geo -2 -v 2 -algo del2d -rand 1.e-9 -clmax 1.e-&
    &1 -order " //ElemOrder)
  CASE(3) !*** STRUCTURED MESH
END SELECT

!*** READ 2D MESH FILE
CALL readmeshEI

!*** ALLOCATE MEMORY
ALLOCATE(am1(50*noukEI),am2(50*noukEI),xvec(noukEI),uEI(noukEI),res(noukEI))

!*** COUNTERS & MATRICES INITIALIZATION
inr=0
NRcnt=1
AMcnt=0
DO i=1,SIZE(am2)
  am2(i)=0.d0 !*** INITIALIZE JACOBIAN MATRIX PATTERN
ENDDO
DO i=1,noukEI
  uEI(i)=bcEI(i) !*** INITIALIZE SOLUTION
  xvec(i)=0.d0 !*** INITIALIZE X VECTOR
ENDDO

!*** FIND SPARSE JACOBIAN MATRIX PATTERN FOR EIKONAL EQUATION
alloc=100
CALL neigElem(nnodesEI,nellEI,nelnodesEI,nopEI,kindEI,neighEI,alloc)
CALL JacobPtrn(noukEI,nelnodesEI,nopEI,kindEI,neighEI)

!*** SOLVE EIKONAL EQUATION
actFlagEI=0.d0
DO WHILE (actFlagEI<=1.d0)
  effViscTermEI=(1.d0-actFlagEI)+ViscTermEI
  DO WHILE (inr==0) !*** START NEWTON-RAPHSON ITERATIONS
    DO i=1,noukEI
      uEI(i)=uEI(i)+xvec(i)
    ENDDO
    NRcnt=NRcnt+1
    CALL axbEI !*** JACOBIAN & RESIDUAL MATRICES ASSEMBLY AND SYSTEM SOLVING 
    CALL ErrorNormCalc(noukEI,xvec,normsum) !*** CALCULATE ERROR NORM
    IF (normsum<=PPres) THEN
      inr=1
    ENDIF
    WRITE(*,'(a,E10.3)')'EI ERROR NORM:',normsum
    IF (normsum>1.d5) THEN
      WRITE(*,'(a)')'EI ERROR NORM > 1.e5 - KILLING PROCESS...'
      STOP
    ENDIF
  ENDDO
  WRITE (*,'(a,i4,/)')'NEWTON ITERATIONS:',NRcnt-1
  WRITE(*,'(a,1x,f6.2,1x,a,/)')'PARAMETRIC SOLVER PROGRESS:',actFlagEI*100.d0,'&
  &%'
!*** RESET COUNTERS AND MATRICES
  inr=0
  NRcnt=1
  DO i=1,noukEI
    xvec(i)=0.d0 !*** RESET X VECTOR
  ENDDO
!*** REMOVE UPPER BOUNDARY CONDITION
  IF (actFlagEI==0.d0) THEN
    DO i=1,nnodesEI
      IF (bcEI(i)==1.d0) THEN
        ncodEI(i)=0
      ENDIF
    ENDDO
  ENDIF
  actFlagEI=actFlagEI+1.d0/REAL(ExtItEI,8)
ENDDO

!*** PRINT RESULTS
eikcnt=0
CALL printres(eikcnt)

!*** CLEAR MEMORY
DEALLOCATE(am1,am2,xvec,res,xptEI,yptEI,nopEI,ncodEI,bcEI,uEI,kindEI,neighEI)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE axbEI

!*******************************************************************************
!*** DESCRIPTION: BOUNDARY CONDITION IMPLEMENTATION AND
!*** SYSTEM SOLVING FOR EIKONAL EQUATION
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::k
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::i1
INTEGER(KIND=8)::icnt1
DOUBLE PRECISION::jacob

!*** MATRICES INITIALIZATION
DO i=1,noukEI
  res(i)=0.d0
ENDDO
DO i=1,SIZE(am1)
  am1(i)=0.d0
ENDDO

!*** ASSEMBLY RESIDUALS (RES) AND JACOBIAN (AM) MATRICES
DO nell=1,nellEI
  CALL abfindEI
ENDDO

!*** IMPOSE DIRICHLET BOUNDARY CONDITIONS
DO i=1,noukEI
  IF (ncodEI(i)==1) THEN
    res(i)=bcEI(i)-uEI(i)
    jacob=0.d0
    DO j=1,kindEI(i)
      DO k=1,nelnodesEI
        i1=nopEI(neighEI(i,j),k)
        icnt1=(i-1)*noukEI+i1
        CALL buildjacob_dir(icnt1,jacob)
      ENDDO
    ENDDO
    jacob=1.d0
    icnt1=(i-1)*noukEI+i
    CALL buildjacob_dir(icnt1,jacob)
  ENDIF
ENDDO

!*** SOLVE THE LINEAR SYSTEM USING MUMPS DIRECT SOLVER
noukLS=noukEI
CALL linsolver

END SUBROUTINE

!*******************************************************************************

SUBROUTINE abfindEI

!*******************************************************************************
!*** DESCRIPTION: RESIDUALS & JACOBIAN MATRICES ASSEMBLY FOR EIKONAL EQUATION
!*******************************************************************************

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
INTEGER(KIND=8)::ngl(9)
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
DOUBLE PRECISION::phi(9)
DOUBLE PRECISION::phic(9)
DOUBLE PRECISION::phie(9)
DOUBLE PRECISION::tphx(9)
DOUBLE PRECISION::tphy(9)

DO i=1,nelnodesEI
  ngl(i)=nopEI(nell,i)
ENDDO

!*** LOOP OVER GAUSS POINTS
DO k=1,NumgpsEI
  DO l=1,NumgpsEI

!*** CALCULATE BASIS FUNCTIONS
    SELECT CASE (eltypeEI)
    CASE(1)
      CALL tsfun2DLIN(gpEI(l),gpEI(l),phi,phic,phie)
    CASE(2)
      CALL tsfun2DQUAD(gpEI(k),gpEI(l),phi,phic,phie)
    END SELECT

!*** INITIALIZE VARIABLES
    x=0.d0; y=0.d0; x1=0.d0; x2=0.d0; y1=0.d0; y2=0.d0; dudx=0.d0; dudy=0.d0

!***ISOPARAMETRIC MAPPING
    DO n=1,nelnodesEI
      x=x+xptEI(ngl(n))*phi(n)
      y=y+yptEI(ngl(n))*phi(n)
      x1=x1+xptEI(ngl(n))*phic(n)
      x2=x2+xptEI(ngl(n))*phie(n)
      y1=y1+yptEI(ngl(n))*phic(n)
      y2=y2+yptEI(ngl(n))*phie(n)
    ENDDO
    dett=x1*y2-x2*y1
    DO n=1,nelnodesEI
      tphx(n)=(y2*phic(n)-y1*phie(n))/dett
      tphy(n)=(x1*phie(n)-x2*phic(n))/dett
      dudx=dudx+uEI(ngl(n))*tphx(n)
      dudy=dudy+uEI(ngl(n))*tphy(n)
    ENDDO

!*** ASSEMBLY RESIDUALS
    DO n=1,nelnodesEI
      n1=ngl(n)
!*** VISCOUS TERM
      res(n1)=res(n1)+effViscTermEI*wgpEI(k)*wgpEI(l)*dett*tphx(n)*dudx+       &
      effViscTermEI*wgpEI(k)*wgpEI(l)*dett*tphy(n)*dudy
!*** EIKONAL EQUATION TERM
      IF (actFlagEI.NE.0.d0) THEN
        res(n1)=res(n1)+actFlagEI*(wgpEI(k)*wgpEI(l)*dett*phi(n)*              &
        SQRT(dudx**2+dudy**2)-wgpEI(k)*wgpEI(l)*dett*phi(n))
      ENDIF
!*** ASSEMBLY JACOBIAN
      DO m=1,nelnodesEI
        m1=ngl(m)
        icnt1=(n1-1)*noukEI+m1
!*** VISCOUS TERM
        jacob=-effViscTermEI*wgpEI(k)*wgpEI(l)*dett*(tphx(n)*tphx(m)+          &
        tphy(n)*tphy(m))
!*** EIKONAL EQUATION TERM
        IF (actFlagEI.NE.0.d0) THEN
          jacob=jacob+actFlagEI*(-((dett*dudx*phi(n)*wgpEI(k)*wgpEI(l))/       &
          (SQRT(dudx**2+dudy**2)))*tphx(m)-((dett*dudy*phi(n)*wgpEI(k)*        &
          wgpEI(l))/(SQRT(dudx**2+dudy**2)))*tphy(m))
        ENDIF
        CALL buildjacob_nor(icnt1,jacob)
      ENDDO
    ENDDO

  ENDDO
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE eikonval(r,th,dudr,dudth,ueii)

!*******************************************************************************
!*** DESCRIPTION: CALCULATE EIKONAL'S VALUE AND DERIVATIVE AT
!*** A GIVEN POINT FOR UNSTRUCTURED MESH
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::k
INTEGER(KIND=8)::n
INTEGER(KIND=8)::inflag
DOUBLE PRECISION::r
DOUBLE PRECISION::th
DOUBLE PRECISION::vec0(2)
DOUBLE PRECISION::vec1(2)
DOUBLE PRECISION::vec2(2)
DOUBLE PRECISION::dot00
DOUBLE PRECISION::dot01
DOUBLE PRECISION::dot02
DOUBLE PRECISION::dot11
DOUBLE PRECISION::dot12
DOUBLE PRECISION::vDenom
DOUBLE PRECISION::uch
DOUBLE PRECISION::vch
DOUBLE PRECISION::gp1
DOUBLE PRECISION::gp2
DOUBLE PRECISION::x12d
DOUBLE PRECISION::x22d
DOUBLE PRECISION::y12d
DOUBLE PRECISION::y22d
DOUBLE PRECISION::dudx
DOUBLE PRECISION::dudy
DOUBLE PRECISION::dett
DOUBLE PRECISION::tphx(6)
DOUBLE PRECISION::tphy(6)
DOUBLE PRECISION::dudr
DOUBLE PRECISION::dudth
DOUBLE PRECISION::phi(6)
DOUBLE PRECISION::phic(6)
DOUBLE PRECISION::phie(6)
DOUBLE PRECISION::xstar
DOUBLE PRECISION::ystar
DOUBLE PRECISION::dsign
DOUBLE PRECISION::ueii

SELECT CASE (geommirEI)
  CASE(1)
    IF (MOD(r*SIN(th),scycleEI)<=scycleEI/2.d0) THEN
      xstar=MOD(r*SIN(th),scycleEI/2.d0)
      dsign=1.d0
    ELSE IF ((MOD(r*SIN(th),scycleEI))>scycleEI/2.d0) THEN
      xstar=scycleEI/2.-MOD(r*SIN(th),scycleEI/2.d0)
      dsign=-1.d0
    ENDIF
  CASE(0)
    xstar=r*SIN(th)
    dsign=1.d0
END SELECT
ystar=r*COS(th)

!*** QUICK SEARCH
inflag=0
DO k=1,nelnodesEI
  DO j=1,kindEI(nopEI(lastelEI,k))
    i=neighEI(nopEI(lastelEI,k),j)
!*** POINT IN A TRIAGLE
!*** COMPUTE VECTORS
    vec0=(/xptEI(nopEI(i,2))-xptEI(nopEI(i,1)),yptEI(nopEI(i,2))-              &
    yptEI(nopEI(i,1))/)
    vec1=(/xptEI(nopEI(i,3))-xptEI(nopEI(i,1)),yptEI(nopEI(i,3))-              &
    yptEI(nopEI(i,1))/)
    vec2=(/xstar-xptEI(nopEI(i,1)),ystar-yptEI(nopEI(i,1))/)
!*** COMPUTE DOT PRODUCTS
    dot00=DOT_PRODUCT(vec0,vec0)
    dot01=DOT_PRODUCT(vec0,vec1)
    dot02=DOT_PRODUCT(vec0,vec2)
    dot11=DOT_PRODUCT(vec1,vec1)
    dot12=DOT_PRODUCT(vec1,vec2)
!*** COMPUTE BARYCENTRIC COORDINATES
    vDenom=1.d0/(dot00*dot11-dot01*dot01)
    uch=(dot11*dot02-dot01*dot12)*vDenom
    vch=(dot00*dot12-dot01*dot02)*vDenom
!*** CHECK IF POINT IS IN TRIANGLE
    IF (uch>=0.AND.vch>=0.AND.uch+vch<1) THEN
      gp1=-((yptEI(nopEI(i,1))-yptEI(nopEI(i,3)))*xstar-(ystar-                &
      yptEI(nopEI(i,3)))*xptEI(nopEI(i,1))+(ystar-yptEI(nopEI(i,1)))*          &
      xptEI(nopEI(i,3)))/((yptEI(nopEI(i,2))-yptEI(nopEI(i,3)))*               &
      xptEI(nopEI(i,1))-(yptEI(nopEI(i,1))-yptEI(nopEI(i,3)))*                 &
      xptEI(nopEI(i,2))+(yptEI(nopEI(i,1))-yptEI(nopEI(i,2)))*                 &
      xptEI(nopEI(i,3)))    
      gp2=((yptEI(nopEI(i,1))-yptEI(nopEI(i,2)))*xstar-(ystar-                 &
      yptEI(nopEI(i,2)))*xptEI(nopEI(i,1))+(ystar-yptEI(nopEI(i,1)))*          &
      xptEI(nopEI(i,2)))/((yptEI(nopEI(i,2))-yptEI(nopEI(i,3)))*               &
      xptEI(nopEI(i,1))-(yptEI(nopEI(i,1))-yptEI(nopEI(i,3)))*                 &
      xptEI(nopEI(i,2))+(yptEI(nopEI(i,1))-yptEI(nopEI(i,2)))*                 &
      xptEI(nopEI(i,3)))
      SELECT CASE (eltypeEI)
        CASE(1)
          CALL tsfun2DLIN(gp1,gp2,phi,phic,phie)
        CASE(2)
          CALL tsfun2DQUAD(gp1,gp2,phi,phic,phie)
      END SELECT
      x12d=0.d0; x22d=0.d0; y12d=0.d0; y22d=0.d0
      dudx=0.d0; dudy=0.d0; ueii=0.d0
      DO n=1,nelnodesEI
        x12d=x12d+xptEI(nopEI(i,n))*phic(n)
        x22d=x22d+xptEI(nopEI(i,n))*phie(n)
        y12d=y12d+yptEI(nopEI(i,n))*phic(n)
        y22d=y22d+yptEI(nopEI(i,n))*phie(n)
        ueii=ueii+uEI(nopEI(i,n))*phi(n)
      ENDDO
      dett=x12d*y22d-x22d*y12d
      DO n=1,nelnodesEI
        tphx(n)=(y22d*phic(n)-y12d*phie(n))/dett
        tphy(n)=(x12d*phie(n)-x22d*phic(n))/dett
      ENDDO
      DO n=1,nelnodesEI
        dudx=dudx+uEI(nopEI(i,n))*tphx(n)
        dudy=dudy+uEI(nopEI(i,n))*tphy(n)
      ENDDO
      dudx=dsign*dudx
      dudr=SIN(th)*dudx+COS(th)*dudy
      dudth=r*COS(th)*dudx-r*SIN(th)*dudy
      lastelEI=i
      inflag=1 !*** SEARCH COMPLETE
      EXIT
    ENDIF
  ENDDO
ENDDO

!*** SLOW SEARCH
IF (inflag==0) THEN
  DO i=1,nellEI
!*** POINT IN A TRIAGLE
!*** COMPUTE VECTORS
    vec0=(/xptEI(nopEI(i,2))-xptEI(nopEI(i,1)),yptEI(nopEI(i,2))-              &
    yptEI(nopEI(i,1))/)
    vec1=(/xptEI(nopEI(i,3))-xptEI(nopEI(i,1)), yptEI(nopEI(i,3))-             &
    yptEI(nopEI(i,1))/)
    vec2=(/xstar-xptEI(nopEI(i,1)),ystar-yptEI(nopEI(i,1))/)
!*** COMPUTE DOT PRODUCTS
    dot00=DOT_PRODUCT(vec0,vec0)
    dot01=DOT_PRODUCT(vec0,vec1)
    dot02=DOT_PRODUCT(vec0,vec2)
    dot11=DOT_PRODUCT(vec1,vec1)
    dot12=DOT_PRODUCT(vec1,vec2)
!*** COMPUTE BARYCENTRIC COORDINATES
    vDenom=1.d0/(dot00*dot11-dot01*dot01)
    uch=(dot11*dot02-dot01*dot12)*vDenom
    vch=(dot00*dot12-dot01*dot02)*vDenom
!*** CHECK IF POINT IS IN TRIANGLE
    IF (uch>=0.AND.vch>=0.AND.uch+vch<1) THEN
      gp1=-((yptEI(nopEI(i,1))-yptEI(nopEI(i,3)))*xstar-                       &
      (ystar-yptEI(nopEI(i,3)))*xptEI(nopEI(i,1))+(ystar-                      &
      yptEI(nopEI(i,1)))*xptEI(nopEI(i,3)))/((yptEI(nopEI(i,2))-               &
      yptEI(nopEI(i,3)))*xptEI(nopEI(i,1))-(yptEI(nopEI(i,1))-                 &
      yptEI(nopEI(i,3)))*xptEI(nopEI(i,2))+(yptEI(nopEI(i,1))-                 &
      yptEI(nopEI(i,2)))*xptEI(nopEI(i,3)))                                                                 
      gp2=((yptEI(nopEI(i,1))-yptEI(nopEI(i,2)))*xstar-                        &
      (ystar-yptEI(nopEI(i,2)))*xptEI(nopEI(i,1))+(ystar-                      &
      yptEI(nopEI(i,1)))*xptEI(nopEI(i,2)))/((yptEI(nopEI(i,2))-               &
      yptEI(nopEI(i,3)))*xptEI(nopEI(i,1))-(yptEI(nopEI(i,1))-                 &
      yptEI(nopEI(i,3)))*xptEI(nopEI(i,2))+(yptEI(nopEI(i,1))-                 &
      yptEI(nopEI(i,2)))*xptEI(nopEI(i,3)))
      SELECT CASE (eltypeEI)
        CASE(1)
          CALL tsfun2DLIN(gp1,gp2,phi,phic,phie)
        CASE(2)
          CALL tsfun2DQUAD(gp1,gp2,phi,phic,phie)
      END SELECT
      x12d=0.d0; x22d=0.d0; y12d=0.d0; y22d=0.d0
      dudx=0.d0; dudy=0.d0; ueii=0.d0
      DO n=1,nelnodesEI
        x12d=x12d+xptEI(nopEI(i,n))*phic(n)
        x22d=x22d+xptEI(nopEI(i,n))*phie(n)
        y12d=y12d+yptEI(nopEI(i,n))*phic(n)
        y22d=y22d+yptEI(nopEI(i,n))*phie(n)
        ueii=ueii+uEI(nopEI(i,n))*phi(n)
      ENDDO
      dett=x12d*y22d-x22d*y12d
      DO n=1,nelnodesEI
        tphx(n)=(y22d*phic(n)-y12d*phie(n))/dett
        tphy(n)=(x12d*phie(n)-x22d*phic(n))/dett
      ENDDO
      DO n=1,nelnodesEI
        dudx=dudx+uEI(nopEI(i,n))*tphx(n)
        dudy=dudy+uEI(nopEI(i,n))*tphy(n)
      ENDDO
      dudx=dsign*dudx
      dudr=SIN(th)*dudx+COS(th)*dudy
      dudth=r*COS(th)*dudx-r*SIN(th)*dudy
      lastelEI=i
      inflag=1 !*** SEARCH COMPLETE
      EXIT
    ENDIF
  ENDDO
ENDIF

END SUBROUTINE

!*******************************************************************************

SUBROUTINE eikonvalstr(r,th,dudr,dudth,ueii)

!*******************************************************************************
!*** DESCRIPTION: CALCULATE EIKONAL'S VALUE AND DERIVATIVE AT 
!*** A GIVEN POINT FOR STRUCTURED MESH
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::k
INTEGER(KIND=8)::n
INTEGER(KIND=8)::inflag
DOUBLE PRECISION::r
DOUBLE PRECISION::th
DOUBLE PRECISION::gp1
DOUBLE PRECISION::gp2
DOUBLE PRECISION::x12d
DOUBLE PRECISION::x22d
DOUBLE PRECISION::y12d
DOUBLE PRECISION::y22d
DOUBLE PRECISION::dudx
DOUBLE PRECISION::dudy
DOUBLE PRECISION::dett
DOUBLE PRECISION::tphx(9)
DOUBLE PRECISION::tphy(9)
DOUBLE PRECISION::dudr
DOUBLE PRECISION::dudth
DOUBLE PRECISION::phi(9)
DOUBLE PRECISION::phic(9)
DOUBLE PRECISION::phie(9)
DOUBLE PRECISION::xstar
DOUBLE PRECISION::ystar
DOUBLE PRECISION::dsign
DOUBLE PRECISION::ueii
DOUBLE PRECISION::xptr
DOUBLE PRECISION::xptl
DOUBLE PRECISION::yptu
DOUBLE PRECISION::yptd
DOUBLE PRECISION::du
DOUBLE PRECISION::dd

SELECT CASE (geommirEI)
  CASE(1)
    IF (MOD(r*SIN(th),scycleEI)<=scycleEI/2.d0) THEN
      xstar=MOD(r*SIN(th),scycleEI/2.d0)
      dsign=1.d0
    ELSE IF ((MOD(r*SIN(th),scycleEI))>scycleEI/2.d0) THEN
      xstar=scycleEI/2.d0-MOD(r*SIN(th),scycleEI/2.d0)
      dsign=-1.d0
    ENDIF
  CASE(0)
    xstar=r*SIN(th)
    dsign=1.d0
END SELECT
ystar=r*COS(th)

!*** QUICK SEARCH
inflag=0
DO k=1,nelnodesEI
  DO j=1,kindEI(nopEI(lastelEI,k))
    i=neighEI(nopEI(lastelEI,k),j)
    xptr=xptEI(nopEI(i,8))
    xptl=xptEI(nopEI(i,2))
    yptu=yptEI(nopEI(i,3))+((yptEI(nopEI(i,9))-yptEI(nopEI(i,3)))/             &
    (xptEI(nopEI(i,9))-xptEI(nopEI(i,3))))*(r*SIN(th)-xptEI(nopEI(i,3)))
    yptd=yptEI(nopEI(i,1))+((yptEI(nopEI(i,7))-yptEI(nopEI(i,1)))/             &
    (xptEI(nopEI(i,7))-xptEI(nopEI(i,1))))*(r*SIN(th)-xptEI(nopEI(i,1)))
!*** CHECK IF POINT IS IN QUADRAGLE
    IF ((xstar<xptr.AND.xstar>xptl).AND.(ystar<yptu.AND.ystar>yptd)) THEN
      du=abs((xptEI(nopEI(i,9))-(xptEI(nopEI(i,3))))*((yptEI(nopEI(i,3))-      &
      r*COS(th)))-(xptEI(nopEI(i,3))-r*SIN(th))*((yptEI(nopEI(i,9))-           &
      yptEI(nopEI(i,3)))))/SQRT((xptEI(nopEI(i,9))-(xptEI(nopEI(i,3))))**2+    &
      (yptEI(nopEI(i,9))-(yptEI(nopEI(i,3))))**2)
      dd=abs((xptEI(nopEI(i,7))-(xptEI(nopEI(i,1))))*((yptEI(nopEI(i,1))-      &
      r*COS(th)))-(xptEI(nopEI(i,1))-r*SIN(th))*((yptEI(nopEI(i,7))-           &
      yptEI(nopEI(i,1)))))/SQRT((xptEI(nopEI(i,7))-(xptEI(nopEI(i,1))))**2+    &
      (yptEI(nopEI(i,7))-(yptEI(nopEI(i,1))))**2)
      gp1=(xstar-xptl)/(xptr-xptl)
      gp2=dd/(dd+du)
      CALL tsfun2DQUADSTR(gp1,gp2,phi,phic,phie)
      x12d=0.d0; x22d=0.d0; y12d=0.d0; y22d=0.d0
      dudx=0.d0; dudy=0.d0; ueii=0.d0
      DO n=1,nelnodesEI
        x12d=x12d+xptEI(nopEI(i,n))*phic(n)
        x22d=x22d+xptEI(nopEI(i,n))*phie(n)
        y12d=y12d+yptEI(nopEI(i,n))*phic(n)
        y22d=y22d+yptEI(nopEI(i,n))*phie(n)
        ueii=ueii+uEI(nopEI(i,n))*phi(n)
      ENDDO
      dett=x12d*y22d-x22d*y12d
      DO n=1,nelnodesEI
        tphx(n)=(y22d*phic(n)-y12d*phie(n))/dett
        tphy(n)=(x12d*phie(n)-x22d*phic(n))/dett
      ENDDO
      DO n=1,nelnodesEI
        dudx=dudx+uEI(nopEI(i,n))*tphx(n)
        dudy=dudy+uEI(nopEI(i,n))*tphy(n)
      ENDDO
      dudx=dsign*dudx
      dudr=SIN(th)*dudx+COS(th)*dudy
      dudth=r*COS(th)*dudx-r*SIN(th)*dudy
      lastelEI=i
      inflag=1 !*** SEARCH COMPLETE
      EXIT
    ENDIF
  ENDDO
ENDDO

!*** SLOW SEARCH
IF (inflag==0) THEN
  DO i=1,nellEI
    xptr=xptEI(nopEI(i,8))
    xptl=xptEI(nopEI(i,2))
    yptu=yptEI(nopEI(i,3))+((yptEI(nopEI(i,9))-yptEI(nopEI(i,3)))/             &
    (xptEI(nopEI(i,9))-xptEI(nopEI(i,3))))*(r*SIN(th)-xptEI(nopEI(i,3)))
    yptd=yptEI(nopEI(i,1))+((yptEI(nopEI(i,7))-yptEI(nopEI(i,1)))/             &
    (xptEI(nopEI(i,7))-xptEI(nopEI(i,1))))*(r*SIN(th)-xptEI(nopEI(i,1)))
!*** CHECK IF POINT IS IN QUADRAGLE
    IF ((xstar<xptr.AND.xstar>xptl).AND.(ystar<yptu.AND.ystar>yptd)) THEN
      du=abs((xptEI(nopEI(i,9))-(xptEI(nopEI(i,3))))*((yptEI(nopEI(i,3))-      &
      r*COS(th)))-(xptEI(nopEI(i,3))-r*SIN(th))*((yptEI(nopEI(i,9))-           &
      yptEI(nopEI(i,3)))))/SQRT((xptEI(nopEI(i,9))-(xptEI(nopEI(i,3))))**2+    &
      (yptEI(nopEI(i,9))-(yptEI(nopEI(i,3))))**2)
      dd=abs((xptEI(nopEI(i,7))-(xptEI(nopEI(i,1))))*((yptEI(nopEI(i,1))-      &
      r*COS(th)))-(xptEI(nopEI(i,1))-r*SIN(th))*((yptEI(nopEI(i,7))-           &
      yptEI(nopEI(i,1)))))/SQRT((xptEI(nopEI(i,7))-(xptEI(nopEI(i,1))))**2+    &
      (yptEI(nopEI(i,7))-(yptEI(nopEI(i,1))))**2)
      gp1=(xstar-xptl)/(xptr-xptl)
      gp2=dd/(dd+du)
      CALL tsfun2DQUADSTR(gp1,gp2,phi,phic,phie)
      x12d=0.d0; x22d=0.d0; y12d=0.d0; y22d=0.d0
      dudx=0.d0; dudy=0.d0; ueii=0.d0
      DO n=1,nelnodesEI
        x12d=x12d+xptEI(nopEI(i,n))*phic(n)
        x22d=x22d+xptEI(nopEI(i,n))*phie(n)
        y12d=y12d+yptEI(nopEI(i,n))*phic(n)
        y22d=y22d+yptEI(nopEI(i,n))*phie(n)
        ueii=ueii+uEI(nopEI(i,n))*phi(n)
      ENDDO
      dett=x12d*y22d-x22d*y12d
      DO n=1,nelnodesEI
        tphx(n)=(y22d*phic(n)-y12d*phie(n))/dett
        tphy(n)=(x12d*phie(n)-x22d*phic(n))/dett
      ENDDO
      DO n=1,nelnodesEI
        dudx=dudx+uEI(nopEI(i,n))*tphx(n)
        dudy=dudy+uEI(nopEI(i,n))*tphy(n)
      ENDDO
      dudx=dsign*dudx
      dudr=SIN(th)*dudx+COS(th)*dudy
      dudth=r*COS(th)*dudx-r*SIN(th)*dudy
      lastelEI=i
      inflag=1 !*** SEARCH COMPLETE
      EXIT
    ENDIF
  ENDDO
ENDIF

END SUBROUTINE

!*******************************************************************************
