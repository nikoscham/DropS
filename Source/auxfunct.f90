!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

FUNCTION binarySearch(value,amflag)

!*******************************************************************************
!*** DESCRIPTION: BINARY SEARCH
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(OUT)::amflag
INTEGER(KIND=8),INTENT(IN)::value
INTEGER(KIND=8)::left
INTEGER(KIND=8)::right
INTEGER(KIND=8)::mid
INTEGER(KIND=8)::binarySearch

amflag=0
left=1
right=AMcnt
DO WHILE (left<=right)
  mid=FLOOR(REAL(right-left,8)/2.d0)+left
  IF (am2(mid)==value) THEN
    binarySearch=mid
    amflag=1
    EXIT
  ELSE IF (value<am2(mid)) THEN
    right=mid-1
  ELSE
    left=mid+1
  ENDIF
ENDDO

END FUNCTION

!*******************************************************************************

MODULE qsort_module

!*******************************************************************************
!*** DESCRIPTION: QUICKSORT ALGORITHM FOR INTEGER MATRICES
!*******************************************************************************

IMPLICIT NONE

CONTAINS

!*******************************************************************************
 
RECURSIVE SUBROUTINE qsort(a,b)
 
INTEGER(KIND=8),INTENT(IN OUT)::a(:)
DOUBLE PRECISION,INTENT(IN OUT)::b(:)
INTEGER(KIND=8)::split

IF (SIZE(a)>1) THEN
  CALL partition_par(a,b,split)
  CALL qsort(a(:split-1),b(:split-1))
  CALL qsort(a(split:),b(split:))
END IF

END SUBROUTINE

!*******************************************************************************

SUBROUTINE partition_par(a,b,marker)

INTEGER(KIND=8),INTENT(IN OUT)::a(:)
DOUBLE PRECISION,INTENT(IN OUT)::b(:)
INTEGER(KIND=8),INTENT(OUT)::marker
INTEGER(KIND=8)::left
INTEGER(KIND=8)::right
INTEGER(KIND=8)::Pivot
INTEGER(KIND=8)::temp
DOUBLE PRECISION::temp2

Pivot=(a(1)+a(SIZE(a)))/2
left=0     
right=SIZE(a)+1

DO WHILE (left<right)
  right=right-1
  DO WHILE (a(right)>Pivot)
    right=right-1
  END DO
  left=left+1
  DO WHILE (a(left)<Pivot)
    left=left+1
  END DO
  IF (left<right) THEN
    temp=a(left)
    a(left)=a(right)
    a(right)=temp
    temp2=b(left)
    b(left)=b(right)
    b(right)=temp2
  END IF
END DO

IF (left==right) THEN
  marker=left+1
ELSE
  marker=left
END IF

END SUBROUTINE

END MODULE

!*******************************************************************************

MODULE qsort_moduleDP

!*******************************************************************************
!*** DESCRIPTION: QUICKSORT ALGORITHM FOR DOUBLE PRECISION MATRICES
!*******************************************************************************

IMPLICIT NONE

CONTAINS

!*******************************************************************************
 
RECURSIVE SUBROUTINE qsortDP(a,b)

DOUBLE PRECISION,INTENT(IN OUT)::a(:)
INTEGER(KIND=8),INTENT(IN OUT)::b(:)
INTEGER(KIND=8)::split

IF(SIZE(a)>1) THEN
  CALL partition_parDP(a,b,split)
  CALL qsortDP(a(:split-1),b(:split-1))
  CALL qsortDP(a(split:),b(split:))
END IF

END SUBROUTINE

!*******************************************************************************
 
SUBROUTINE partition_parDP(a,b,marker)
 
DOUBLE PRECISION,INTENT(IN OUT)::a(:)
INTEGER(KIND=8),INTENT(IN OUT)::b(:)
INTEGER(KIND=8),INTENT(OUT)::marker
INTEGER(KIND=8)::left
INTEGER(KIND=8)::right
INTEGER(KIND=8)::temp2
DOUBLE PRECISION::Pivot
DOUBLE PRECISION::temp

Pivot=(a(1)+a(SIZE(a)))/2
left=0     
right=SIZE(a)+1

DO WHILE (left<right)
  right=right-1
  DO WHILE (a(right)>Pivot)
    right=right-1
  END DO
  left=left+1
  DO WHILE (a(left)<Pivot)
    left=left+1
  END DO
  IF (left<right) THEN
    temp=a(left)
    a(left)=a(right)
    a(right)=temp
    temp2=b(left)
    b(left)=b(right)
    b(right)=temp2
  END IF
END DO
IF (left==right) THEN
  marker=left+1
ELSE
  marker=left
END IF

END SUBROUTINE

END MODULE

!*******************************************************************************

SUBROUTINE ErrorNormCalc(nouk,x,error)

!*******************************************************************************
!*** DESCRIPTION: ERROR NORM CALCULATION
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nouk
DOUBLE PRECISION,INTENT(IN)::x(nouk)
DOUBLE PRECISION,INTENT(OUT)::error
INTEGER(KIND=8)::i

error=0.d0
DO i=1,nouk
  error=error+(x(i))**2
ENDDO
error=SQRT(error)
!error=SQRT(error)/nouk

END SUBROUTINE

!*******************************************************************************

SUBROUTINE buildjacob_nor(icnt1,jacob)

!*******************************************************************************
!*** DESCRIPTION: JACOBIAN MATRIX BUILDER (COORDINATES FORMAT)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::amflag
INTEGER(KIND=8)::binarySearch
INTEGER(KIND=8)::k
INTEGER(KIND=8),INTENT(IN)::icnt1
DOUBLE PRECISION,INTENT(IN)::jacob

SELECT CASE (NRcnt)
  CASE(1)
    amflag=0
    DO k=AMcnt,1,-1
      IF (icnt1==am2(k)) THEN
        am1(k)=am1(k)+jacob
        amflag=1
        EXIT
      ENDIF
    ENDDO
    IF (amflag==0) THEN
      AMcnt=AMcnt+1
      am1(AMcnt)=jacob
      am2(AMcnt)=icnt1
    ENDIF
  CASE DEFAULT
    k=binarySearch(icnt1,amflag)
    am1(k)=am1(k)+jacob
END SELECT

END SUBROUTINE

!*******************************************************************************

SUBROUTINE buildjacob_dir(icnt1,jacob)

!*******************************************************************************
!*** DESCRIPTION: JACOBIAN MATRIX BUILDER (FOR DIRICHLET BOUNDARY CONDITIONS)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::amflag
INTEGER(KIND=8)::binarySearch
INTEGER(KIND=8)::k
INTEGER(KIND=8),INTENT(IN)::icnt1
DOUBLE PRECISION,INTENT(IN)::jacob

k=binarySearch(icnt1,amflag)
IF (amflag==1) THEN
  am1(k)=jacob
ENDIF

END SUBROUTINE

!*******************************************************************************

SUBROUTINE linsolver

!*******************************************************************************
!*** DESCRIPTION: LINEAR SYSTEM SOLVER
!*******************************************************************************

USE CommonVars
INCLUDE 'dmumps_struc.h'
TYPE (dmumps_struc) mumps_par
INTEGER(KIND=8)::i
INTEGER(KIND=8)::i1
INTEGER(KIND=8)::i2

!*** MUMPS SOLVER INITIALIZATION
mumps_par%JOB=-1
mumps_par%SYM=0
mumps_par%PAR=1
CALL dmumps(mumps_par)
mumps_par%N=noukLS
mumps_par%NZ=AMcnt
ALLOCATE(mumps_par%IRN(mumps_par%NZ))
ALLOCATE(mumps_par%JCN(mumps_par%NZ))
ALLOCATE(mumps_par%A(mumps_par%NZ))
ALLOCATE(mumps_par%RHS(mumps_par%N))
DO i=1,mumps_par%NZ
  i1=INT(am2(i)/mumps_par%N+1.d0,8)
  i2=INT(MOD(am2(i),mumps_par%N),8)
  IF (i2==0) THEN
    i1=i1-1
    i2=mumps_par%N
  ENDIF
  mumps_par%IRN(i)=i1
  mumps_par%JCN(i)=i2
  mumps_par%A(i)=am1(i)
ENDDO
DO i=1,mumps_par%N
  mumps_par%RHS(i)=res(i)
ENDDO

!*** SOLVE SYSTEM
mumps_par%ICNTL(3)=0
mumps_par%ICNTL(4)=1
mumps_par%JOB=6
CALL dmumps(mumps_par)
DO i=1,mumps_par%N
  xvec(i)=mumps_par%RHS(i)
ENDDO

!*** CLEAR MEMORY
DEALLOCATE(mumps_par%IRN)
DEALLOCATE(mumps_par%JCN)
DEALLOCATE(mumps_par%A)
DEALLOCATE(mumps_par%RHS)
mumps_par%JOB=-2
CALL dmumps(mumps_par)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE NumCurvCalc(Tnodes,Ts,Tth,Tr,Tcurv)

!*******************************************************************************
!*** DESCRIPTION: MEAN CURVATURE CALCULATION (CENTRAL DIFFERENCES)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE

INTEGER(KIND=8),INTENT(IN)::Tnodes
DOUBLE PRECISION,INTENT(IN)::Ts(Tnodes)
DOUBLE PRECISION,INTENT(IN)::Tth(Tnodes)
DOUBLE PRECISION,INTENT(IN)::Tr(Tnodes)
DOUBLE PRECISION,INTENT(OUT)::Tcurv(Tnodes)
INTEGER(KIND=8)::i
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::dthds
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::drds
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::dthdsds
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::drdsds
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::P1
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::P2
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::P3
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:)::P4

!***MATRICES INITIALIZATION
ALLOCATE(P1(Tnodes),P2(Tnodes),P3(Tnodes),P4(Tnodes),dthds(Tnodes),            &
drds(Tnodes),dthdsds(Tnodes),drdsds(Tnodes))
DO i=1,Tnodes
  dthds(i)=0.d0; drds(i)=0.d0; dthdsds(i)=0.d0; drdsds(i)=0.d0
  P1(i)=0.d0; P2(i)=0.d0; P3(i)=0.d0; P4(i)=0.d0; Tcurv(i)=0.d0
ENDDO

DO i=2,Tnodes-1
  dthds(i)=(Tth(i+1)-Tth(i-1))/((Ts(i+1)-Ts(i-1)))
  drds(i)=(Tr(i+1)-Tr(i-1))/((Ts(i+1)-Ts(i-1)))
ENDDO

DO i=4,Tnodes-3
  dthdsds(i)=(dthds(i+1)-dthds(i-1))/((Ts(i+1)-Ts(i-1)))
  drdsds(i)=(drds(i+1)-drds(i-1))/((Ts(i+1)-Ts(i-1)))
  P1(i)=1.d0/(SQRT(Tr(i)**2*dthds(i)**2+drds(i)**2))
  P2(i)=dthds(i)
  P3(i)=drds(i)**2*dthds(i)/(Tr(i)**2*dthds(i)**2+drds(i)**2)
  P4(i)=(Tr(i)*dthdsds(i)*drds(i)-drdsds(i)*Tr(i)*dthds(i))/(dthds(i)**2*      &
  (Tr(i)**2+drds(i)**2/dthds(i)**2))
  Tcurv(i)=P1(i)*(P2(i)+P3(i)+P4(i)) !*** MEAN CURVATURE
ENDDO
Tcurv(1)=Tcurv(3)
Tcurv(2)=Tcurv(3)
Tcurv(3)=Tcurv(4)
Tcurv(Tnodes-2)=Tcurv(Tnodes-3)
Tcurv(Tnodes-1)=Tcurv(Tnodes-3)
Tcurv(Tnodes)=Tcurv(Tnodes-3)

DEALLOCATE(P1,P2,P3,P4,dthds,drds,dthdsds,drdsds)

END SUBROUTINE

!*******************************************************************************

DOUBLE PRECISION FUNCTION NumJacPert(x)

!*******************************************************************************
!*** DESCRIPTION: COMPUTE PERTUBATION FOR NUMERICAL JACOBIAN
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
DOUBLE PRECISION,INTENT(IN)::x
DOUBLE PRECISION::eps
DOUBLE PRECISION::Temp
DOUBLE PRECISION::Meps

Meps=1.e-16
Temp=X
eps=SQRT(Meps)
NumJacPert=eps*ABS(Temp)
IF (NumJacPert<eps) NumJacPert=eps

END FUNCTION

!*******************************************************************************

