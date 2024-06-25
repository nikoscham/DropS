!************************************************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com) @ NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!************************************************************************************************************

MODULE initsol

CONTAINS

!************************************************************************************************************

SUBROUTINE Sol1DRead(FileID,nouk,nnodes,u,spt)

!************************************************************************************************************
!*** DESCRIPTION: READ 1D SOLUTION FOR THE AUGMENTED YL EQUATION
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::FileID
INTEGER(KIND=8),INTENT(IN)::nouk
INTEGER(KIND=8),INTENT(IN)::nnodes
DOUBLE PRECISION,INTENT(OUT)::u(nouk)
DOUBLE PRECISION,INTENT(OUT)::spt(nnodes)
DOUBLE PRECISION::null1
INTEGER(KIND=8)::i

DO i=1,nouk
  READ(FileID,*)u(i)
ENDDO
IF (SolverID==2) THEN !*** SKIP PSEUDO ARC-LENGTH VALUE
  READ(FileID,*)null1
ENDIF
DO i=1,nnodes
  READ(FileID,*)spt(i)
ENDDO
IF (SolverID==6) THEN
  READ(FileID,*)VltEL
ENDIF
CLOSE(FileID)

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE Sol2DRead(FileID,xpt,ypt,u,nop,nnodes,nelli,eltype)

!************************************************************************************************************
!*** DESCRIPTION: READ 2D SOLUTION
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::FileID
INTEGER(KIND=8),INTENT(IN)::eltype
INTEGER(KIND=8),INTENT(OUT)::nnodes
INTEGER(KIND=8),INTENT(OUT)::nelli
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)::nop
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(OUT)::xpt
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(OUT)::ypt
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(OUT)::u
INTEGER(KIND=8)::i

READ(FileID,*)nnodes !*** READ NUMBER OF NODES
ALLOCATE(xpt(nnodes),ypt(nnodes),u(nnodes))
DO i=1,nnodes !*** READ NODES COORDINATES
  READ(FileID,*)xpt(i),ypt(i),u(i)
ENDDO
READ(FileID,*)nelli !*** READ NUMBER OF ELEMENTS
SELECT CASE (eltype)
  CASE(1)
    ALLOCATE(nop(nelli,3))
    DO i=1,nelli !*** READ NOP ARRAY FOR LINEAR TRIANGULAR ELEMENTS
      READ(FileID,*)nop(i,1),nop(i,2),nop(i,3)
    ENDDO
  CASE(2)
    ALLOCATE(nop(nelli,6))
    DO i=1,nelli !*** READ NOP ARRAY FOR QUADRATIC TRIANGULAR ELEMENTS
      READ(FileID,*)nop(i,1),nop(i,2),nop(i,3),nop(i,4),nop(i,5),nop(i,6)
    ENDDO
  CASE(3)
    ALLOCATE(nop(nelli,9))
    DO i=1,nelli !*** READ NOP ARRAY FOR QUADRATIC RECTANGULAR ELEMENTS
      READ(FileID,*)nop(i,1),nop(i,2),nop(i,3),nop(i,4),nop(i,5),nop(i,6),nop(i,7),nop(i,8),nop(i,9)
    ENDDO
END SELECT
CLOSE(FileID)

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE interpsol1D(uinit,xcoordinit,nnodesinit,u,xcoord,nnodes)

!************************************************************************************************************
!*** DESCRIPTION: PERFORM INTERPOLATION FOR THE 1D AUGMENTED YL SOLUTION
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nnodes
INTEGER(KIND=8),INTENT(IN)::nnodesinit
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(IN)::xcoord
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(IN)::xcoordinit
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(IN)::uinit
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:),INTENT(OUT)::u
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
DOUBLE PRECISION::gpx
DOUBLE PRECISION::ph(2)
DOUBLE PRECISION::phd(2)

ALLOCATE(u(nnodes))
DO i=1,nnodes
  DO j=2,nnodesinit
    IF (xcoord(i)<xcoordinit(j).AND.xcoord(i)>xcoordinit(j-1)) THEN
      gpx=(xcoord(i)-xcoordinit(j-1))/(xcoordinit(j)-xcoordinit(j-1))
      CALL tsfun1DLIN(gpx,ph,phd)
      u(i)=uinit(j-1)*ph(1)
      u(i)=u(i)+uinit(j)*ph(2)
    ELSE IF (xcoord(i)==xcoordinit(j)) THEN
      u(i)=uinit(j)
    ELSE IF (xcoord(i)==xcoordinit(j-1)) THEN
      u(i)=uinit(j-1)
    ENDIF
  ENDDO
ENDDO

END SUBROUTINE

!************************************************************************************************************

END MODULE

!************************************************************************************************************

SUBROUTINE efextr

!************************************************************************************************************
!*** DESCRIPTION: CALCULATE ELECTRID FIELD STRENGTH AT THE AUGMENTED YL GAUSS POINTS
!************************************************************************************************************

USE CommonVars
USE qsort_moduleDP
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::k
INTEGER(KIND=8)::n
INTEGER(KIND=8)::elem
INTEGER(KIND=8)::inflag
INTEGER(KIND=8)::lastel
INTEGER(KIND=8)::eleband
INTEGER,ALLOCATABLE,DIMENSION(:)::boundele
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::ellnum
DOUBLE PRECISION::xstar
DOUBLE PRECISION::ystar
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
DOUBLE PRECISION::x
DOUBLE PRECISION::y
DOUBLE PRECISION::x12
DOUBLE PRECISION::x22
DOUBLE PRECISION::y12
DOUBLE PRECISION::y22
DOUBLE PRECISION::dudx
DOUBLE PRECISION::dudy
DOUBLE PRECISION::dett
DOUBLE PRECISION::tphx(6)
DOUBLE PRECISION::tphy(6)
DOUBLE PRECISION::phi(6)
DOUBLE PRECISION::phic(6)
DOUBLE PRECISION::phie(6)
DOUBLE PRECISION::xpos(nellYL)
DOUBLE PRECISION::ypos(nellYL)
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::sptelle
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::barydist

!*** FIND THE GAUSS POINTS NEIGHBORING ELEMENTS
ALLOCATE(boundele(10*nnodesYL))
k=0
DO i=1,nnodesYL
  DO j=1,kindEL(YLnodes(i))
    k=k+1
    boundele(k)=neighEL(YLnodes(i),j)
  ENDDO
ENDDO
!*** DELETE COMMON ELEMENTS
DO i=1,k
  DO j=i,k
    IF (boundele(i)==boundele(j).AND.i.NE.j) THEN
      boundele(j)=0
    ENDIF
  ENDDO
ENDDO
j=0
!*** REORDER THE NEIGHBORING ELEMENTS MATRIX
DO i=1,k
  IF (boundele(i).NE.0) THEN
    j=j+1
    boundele(j)=boundele(i)
  ENDIF
ENDDO
ALLOCATE(barydist(j),ellnum(j))

eleband=10
lastel=eleband+1
DO i=1,nellYL
  inflag=0
!*** LOAD GAUSS POINT COORDINATES
  xstar=xptYLGP(i)
  ystar=yptYLGP(i)

!*** QUICK SEARCH
  DO k=lastel-eleband,lastel+eleband
!*** POINT IN A TRIAGLE
!*** COMPUTE VECTORS
    elem=boundele(k)
    vec0=(/xptEL(nopEL(elem,2))-xptEL(nopEL(elem,1)),yptEL(nopEL(elem,2))-yptEL(nopEL(elem,1))/)
    vec1=(/xptEL(nopEL(elem,3))-xptEL(nopEL(elem,1)),yptEL(nopEL(elem,3))-yptEL(nopEL(elem,1))/)
    vec2=(/xstar-xptEL(nopEL(elem,1)),ystar-yptEL(nopEL(elem,1))/)
!*** COMPUTE DOT PRODUCTS
    dot00=DOT_PRODUCT(vec0,vec0)
    dot01=DOT_PRODUCT(vec0,vec1)
    dot02=DOT_PRODUCT(vec0,vec2)
    dot11=DOT_PRODUCT(vec1,vec1)
    dot12=DOT_PRODUCT(vec1,vec2)
!*** COMPUTE BARYCENTRIC COORDINATES
    vDenom=1./(dot00*dot11-dot01*dot01)
    uch=(dot11*dot02-dot01*dot12)*vDenom
    vch=(dot00*dot12-dot01*dot02)*vDenom
!*** CHECK IF POINT IS IN TRIANGLE
    IF (uch>=0.AND.vch>=0.AND.uch+vch<1) THEN
      gp1=0.5d0
      gp2=0.d0
      SELECT CASE (eltypeEL)
        CASE(1)
          CALL tsfun2DLIN(gp1,gp2,phi,phic,phie)
        CASE(2)
          CALL tsfun2DQUAD(gp1,gp2,phi,phic,phie)
      END SELECT
      x=0.d0; y=0.d0; x12=0.d0; x22=0.d0; y12=0.d0; y22=0.d0; dudx=0.d0; dudy=0.d0
      DO n=1,nelnodesEL
        x=x+xptEL(nopEL(elem,n))*phi(n)
        y=y+yptEL(nopEL(elem,n))*phi(n)
        x12=x12+xptEL(nopEL(elem,n))*phic(n)
        x22=x22+xptEL(nopEL(elem,n))*phie(n)
        y12=y12+yptEL(nopEL(elem,n))*phic(n)
        y22=y22+yptEL(nopEL(elem,n))*phie(n)
      ENDDO
      dett=x12*y22-x22*y12
      DO n=1,nelnodesEL
        tphx(n)=(y22*phic(n)-y12*phie(n))/dett
        tphy(n)=(x12*phie(n)-x22*phic(n))/dett      
      ENDDO
      DO n=1,nelnodesEL
        dudx=dudx+uEL(nopEL(elem,n))*tphx(n)
        dudy=dudy+uEL(nopEL(elem,n))*tphy(n)
      ENDDO
      efieldYLGP(i)=SQRT(dudx**2+dudy**2)**2
      xpos(i)=x; ypos(i)=y
      inflag=1
      IF (k<j-(eleband+1).AND.k>(eleband+1)) THEN
        lastel=k
      ENDIF
      EXIT
    ENDIF
  ENDDO

!*** SLOW SEARCH
  IF (inflag==0) THEN
    gp1=1.d0/3.d0
    gp2=1.d0/3.d0
    SELECT CASE (eltypeEL)
      CASE(1)
        CALL tsfun2DLIN(gp1,gp2,phi,phic,phie)
      CASE(2)
        CALL tsfun2DQUAD(gp1,gp2,phi,phic,phie)
    END SELECT
    DO k=1,j
      barydist(k)=10.d0
      ellnum(k)=k
    ENDDO
    DO k=1,j
      elem=boundele(k)
      x=0.d0
      y=0.d0
      DO n=1,nelnodesEL
        x=x+xptEL(nopEL(elem,n))*phi(n)
        y=y+yptEL(nopEL(elem,n))*phi(n)
      ENDDO
      barydist(k)=SQRT((x-xstar)**2+(y-ystar)**2)
    ENDDO
    CALL qsortDP(barydist,ellnum)
    elem=boundele(ellnum(1))
    gp1=0.5d0
    gp2=0.d0
    SELECT CASE (eltypeEL)
      CASE(1)
        CALL tsfun2DLIN(gp1,gp2,phi,phic,phie)
      CASE(2)
        CALL tsfun2DQUAD(gp1,gp2,phi,phic,phie)
    END SELECT
    x=0.d0; y=0.d0; x12=0.d0; x22=0.d0; y12=0.d0; y22=0.d0; dudx=0.d0; dudy=0.d0
    DO n=1,nelnodesEL
      x=x+xptEL(nopEL(elem,n))*phi(n)
      y=y+yptEL(nopEL(elem,n))*phi(n)
      x12=x12+xptEL(nopEL(elem,n))*phic(n)
      x22=x22+xptEL(nopEL(elem,n))*phie(n)
      y12=y12+yptEL(nopEL(elem,n))*phic(n)
      y22=y22+yptEL(nopEL(elem,n))*phie(n)
    ENDDO
    dett=x12*y22-x22*y12
    DO n=1,nelnodesEL
      tphx(n)=(y22*phic(n)-y12*phie(n))/dett
      tphy(n)=(x12*phie(n)-x22*phic(n))/dett
    ENDDO
    DO n=1,nelnodesEL
      dudx=dudx+uEL(nopEL(elem,n))*tphx(n)
      dudy=dudy+uEL(nopEL(elem,n))*tphy(n)
    ENDDO
    efieldYLGP(i)=SQRT(dudx**2+dudy**2)**2
    xpos(i)=x; ypos(i)=y
    inflag=1
    IF (ellnum(1)<j-(eleband+1).AND.ellnum(1)>(eleband+1)) THEN
      lastel=ellnum(1)
    ENDIF
  ENDIF

ENDDO

DEALLOCATE(boundele,barydist,ellnum) !*** CLEAR MEMORY

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE wlsCalc(thYoungExtYL,WParamExtYL)

!************************************************************************************************************
!*** DESCRIPTION: CALCULATE WETTING PARAMETER VALUE (FRUMKIN-DERJAGUIN APPROACH)
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::thYoungExtYL
DOUBLE PRECISION,INTENT(OUT)::WParamExtYL

WParamExtYL=(COS(thYoungExtYL*Pi/180)+1.d0)/((sigmaYL/(C2YL-1)-sigmaYL/(C1YL-1)-(flmYL*sigmaYL)/(C3YL-1)))
WParamExtYL=WParamExtYL/4.d0

END SUBROUTINE

!************************************************************************************************************

SUBROUTINE thYoungCalc(thYoungExtYL,WParamExtYL)

!************************************************************************************************************
!*** DESCRIPTION: CALCULATE YOUNG CONTACT ANGLE (FRUMKIN-DERJAGUIN APPROACH)
!************************************************************************************************************

USE CommonVars
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::WParamExtYL
DOUBLE PRECISION,INTENT(OUT)::thYoungExtYL
DOUBLE PRECISION::WParamIntYL

WParamIntYL=4.d0*WParamExtYL
thYoungExtYL=180.d0*ACOS((sigmaYL/(C2YL-1)-sigmaYL/(C1YL-1)-(flmYL*sigmaYL)/(C3YL-1))*WParamIntYL-1.d0)/Pi

END SUBROUTINE

!************************************************************************************************************
