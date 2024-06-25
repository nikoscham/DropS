!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

SUBROUTINE printres(fcnt)

!*******************************************************************************
!*** DESCRIPTION: PRINT RESULTS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::fcnt
DOUBLE PRECISION::zcordi
CHARACTER(LEN=5)::FileNumber
CHARACTER(LEN=18)::FileName
CHARACTER(LEN=9)::dataout
CHARACTER(LEN=9)::parav

!**** INITIALIZATION
zcordi=0.d0
WRITE(FileNumber,'(i5.5)')fcnt

!*** PRINT .VTK FILE
IF (SolverID==4.OR.SolverID==5.OR.SolverID==6) THEN
  parav='res_parav'
  FileName=parav//FileNumber//'.vtk'
  OPEN(10,FILE=FileName)
  SELECT CASE (eltypeEL)
    CASE(1) !*** LINEAR TRIAGULAR ELEMENTS
      WRITE(10,'(a)')'# vtk DataFile Version 3.0'
      WRITE(10,'(a)')'DropS RESULTS'
      WRITE(10,'(a,/)')'ASCII'
      WRITE(10,'(a)')'DATASET UNSTRUCTURED_GRID'
      WRITE(10,3010)'POINTS',noukEL
      DO i=1,noukEL
        WRITE(10,3020)xptEL(i),yptEL(i),zcordi
      ENDDO
      WRITE(10,3030)'CELLS',nellEL,nellEL*4
      DO i=1,nellEL
        WRITE(10,3040)nopEL(i,1)-1,nopEL(i,2)-1,nopEL(i,3)-1
      ENDDO
      WRITE(10,3050)'CELL_TYPES',nellEL
      DO i =1,nellEL
        WRITE(10,'(i1)')5
      ENDDO
      WRITE(10,3060)'POINT_DATA',noukEL
      WRITE(10,'(a)')'SCALARS ELECTRIC_POTENTIAL float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,noukEL
        WRITE(10,'(f14.10)')uEL(i)
      ENDDO
      WRITE(10,'(a,i7)')'CELL_DATA ',nellEL
      WRITE(10,'(a)')'SCALARS ELECTRIC_FIELD float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,nellEL
        WRITE(10,'(f14.10)')efieldELElm(i)
      ENDDO
    CASE(2) !*** QUADRATIC TRIAGULAR ELEMENTS
      WRITE(10,'(a)')'# vtk DataFile Version 3.0'
      WRITE(10,'(a)')'DropS RESULTS'
      WRITE(10,'(a,/)')'ASCII'
      WRITE(10,'(a)')'DATASET UNSTRUCTURED_GRID'
      WRITE(10,3010)'POINTS',noukEL
      DO i=1,noukEL
        WRITE(10,3020)xptEL(i),yptEL(i),zcordi
      ENDDO
      WRITE(10,3030)'CELLS',nellEL,nellEL*7
      DO i=1,nellEL
        WRITE(10,3070)nopEL(i,1)-1,nopEL(i,2)-1,nopEL(i,3)-1,nopEL(i,4)-1,     &
        nopEL(i,5)-1,nopEL(i,6)-1
      ENDDO
      WRITE(10,3050)'CELL_TYPES',nellEL
      DO i=1,nellEL
        WRITE(10,'(i2)')22
      ENDDO
      WRITE(10,3060)'POINT_DATA',noukEL
      WRITE(10,'(a)')'SCALARS ELECTRIC_POTENTIAL float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,noukEL
        WRITE(10,'(f14.10)')uEL(i)
      ENDDO
      WRITE(10,'(a,i7)')'CELL_DATA ',nellEL
      WRITE(10,'(a)')'SCALARS ELECTRIC_FIELD float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,nellEL
        WRITE(10,'(f14.10)')efieldELElm(i)
      ENDDO
  END SELECT
  CLOSE(10)
ELSE IF (SolverID==1) THEN
  parav='res_parav'
  FileName=parav//FileNumber//'.vtk'
  OPEN(10,FILE=FileName)
  SELECT CASE (eltypeEI)
    CASE(1)
      WRITE(10,'(a)')'# vtk DataFile Version 3.0'
      WRITE(10,'(a)')'DropS RESULTS'
      WRITE(10,'(a,/)')'ASCII'
      WRITE(10,'(a)')'DATASET UNSTRUCTURED_GRID'
      WRITE(10,3010)'POINTS',noukEI
      DO i=1,nnodesEI
        WRITE(10,3020)xptEI(i),yptEI(i),zcordi
      ENDDO
      WRITE(10,3030)'CELLS',nellEI,nellEI*4
      DO i=1,nellEI
        WRITE(10,3040)nopEI(i,1)-1,nopEI(i,2)-1,nopEI(i,3)-1
      ENDDO
      WRITE(10,3050)'CELL_TYPES',nellEI
      DO i=1,nellEI
        WRITE(10,'(i1)')5
      ENDDO
      WRITE(10,3060)'POINT_DATA',noukEI
      WRITE(10,'(a)')'SCALARS DISTANCE_FROM_SOLID float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,noukEI
        WRITE(10,'(f14.10)')uEI(i)
      ENDDO
    CASE(2)
      WRITE(10,'(a)')'# vtk DataFile Version 3.0'
      WRITE(10,'(a)')'DropS RESULTS'
      WRITE(10,'(a,/)')'ASCII'
      WRITE(10,'(a)')'DATASET UNSTRUCTURED_GRID'
      WRITE(10,3010)'POINTS',noukEI
      DO i=1,nnodesEI
        WRITE(10,3020)xptEI(i),yptEI(i),zcordi
      ENDDO
      WRITE(10,3030)'CELLS',nellEI,nellEI*7
      DO i=1,nellEI
        WRITE(10,3070)nopEI(i,1)-1,nopEI(i,2)-1,nopEI(i,3)-1,nopEI(i,4)-1,     &
        nopEI(i,5)-1,nopEI(i,6)-1
      ENDDO
      WRITE(10,3050)'CELL_TYPES',nellEI
      DO i=1,nellEI
        WRITE(10,'(i2)')22
      ENDDO
      WRITE(10,3060)'POINT_DATA',noukEI
      WRITE(10,'(a)')'SCALARS DISTANCE_FROM_SOLID float 1'
      WRITE(10,'(a)')'LOOKUP_TABLE default'
      DO i=1,noukEI
        WRITE(10,'(f14.10)')uEI(i)
      ENDDO
  END SELECT
  CLOSE(10)
ENDIF

!*** PRINT OUTPUT FILE
IF (SolverID==2.OR.SolverID==3.OR.SolverID==4.OR.SolverID==5.OR.SolverID==6.OR.SolverID==7)   &
  THEN
  dataout='res_dtout'
  FileName=dataout//FileNumber//'.txt'
  OPEN(20,FILE=FileName)
  DO i=1,noukYL
    WRITE(20,'(f14.9)')uYL(i)
  ENDDO
!  WRITE(20,'(f14.9)')WParamYL
  DO i=1,nnodesYL
    WRITE(20,'(f14.11)')sptYL(i)
!    WRITE(20,'(f14.9,f14.9)')uYL(i)*sin(uYL(nnodesYL+i)),uYL(i)*cos(uYL(nnodesYL+i))
  ENDDO
  CLOSE(20)
ELSE IF (SolverID==8) THEN
  dataout='res_dtout'
  FileName=dataout//FileNumber//'.txt'
  OPEN(20,FILE=FileName)
  DO i=1,noukYL
    WRITE(20,'(f14.9)')uYL(i)
  ENDDO
  DO i=1,nnodesYL
    WRITE(20,'(f14.9)')thptYL(i)
!    WRITE(20,'(f14.9,f14.9)')uYL(i)*sin(thptYL(i)),uYL(i)*cos(thptYL(i))
  ENDDO
  CLOSE(20)
ENDIF

!*** PRINT EIKONAL SOLUTION FILE
IF (SolverID==1) THEN
  OPEN (30,FILE='EikonalSolution.txt')
  WRITE(30,'(i7)')noukEI
  DO i=1,noukEI
    WRITE(30,3020)xptEI(i),yptEI(i),uEI(i)
  ENDDO
  WRITE(30,'(i7)')nellEI
  SELECT CASE (eltypeEI)
    CASE(1)
      DO i=1,nellEI
        WRITE(30,3080)nopEI(i,1),nopEI(i,2),nopEI(i,3)
      ENDDO
    CASE(2)
      DO i=1,nellEI
        WRITE(30,3090)nopEI(i,1),nopEI(i,2),nopEI(i,3),nopEI(i,4),nopEI(i,5),  &
        nopEI(i,6)
      ENDDO
  END SELECT
  CLOSE(30)
ENDIF

3010  FORMAT(a6,1x,i8,1x,'float')
3020  FORMAT(f14.11,1x,f14.11,1x,f14.11)
3030  FORMAT(/,a5,1x,i7,1x,i8)
3040  FORMAT('3',1x,i7,1x,i7,1x,i7)
3050  FORMAT(/,a10,1x,i7)
3060  FORMAT(/,a10,1x,i7)
3070  FORMAT('6',1x,i7,1x,i7,1x,i7,1x,i7,1x,i7,1x,i7)
3080  FORMAT(i7,1x,i7,1x,i7)
3090  FORMAT(i7,1x,i7,1x,i7,1x,i7,1x,i7,1x,i7)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE elfieldElms

!*******************************************************************************
!*** DESCRIPTION: ELECTRIC FIELD STRENGTH CALCULATION
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::n
INTEGER(KIND=8)::nelle
INTEGER(KIND=8)::ngl(6)
DOUBLE PRECISION::gpx
DOUBLE PRECISION::dett
DOUBLE PRECISION::x
DOUBLE PRECISION::y
DOUBLE PRECISION::x1
DOUBLE PRECISION::y1
DOUBLE PRECISION::x2
DOUBLE PRECISION::y2
DOUBLE PRECISION::dudx
DOUBLE PRECISION::dudy
DOUBLE PRECISION::tphx(6)
DOUBLE PRECISION::tphy(6)
DOUBLE PRECISION::phi(6)
DOUBLE PRECISION::phic(6)
DOUBLE PRECISION::phie(6)

!*** INITIALIZE
ALLOCATE(efieldELElm(nellEL))

!*** BASIS FUNCTIONS
gpx=1.d0/3.d0
SELECT CASE (eltypeEL)
CASE(1)
  CALL tsfun2DLIN(gpx,gpx,phi,phic,phie)
CASE(2)
  CALL tsfun2DQUAD(gpx,gpx,phi,phic,phie)
END SELECT

DO nelle=1,nellEL
  DO i=1,nelnodesEL
    ngl(i)=nopEL(nelle,i)
  ENDDO
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
  efieldELElm(nelle)=SQRT(dudx**2+dudy**2) !*** ELECTRIC FIELD STRENGTH
ENDDO

END SUBROUTINE

!*******************************************************************************
