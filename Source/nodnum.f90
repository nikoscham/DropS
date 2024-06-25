!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

MODULE nodnum

CONTAINS

!*******************************************************************************

SUBROUTINE nodnumbAYL

!*******************************************************************************
!*** DESCRIPTION: NODAL NUMBERING (AUGMENTED YOUNG-LAPLACE COMPUTATIONAL DOMAIN)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j

DO i=1,nellYL
  DO j=1,2
    nopYL(i,j)=i+(j-1)
  ENDDO
  DO j=3,4
    nopYL(i,j)=nopYL(i,j-2)+nnodesYL
  ENDDO
ENDDO
DO i=1,nellYL
  nopYL(i,5)=2*nnodesYL+1
ENDDO
DO i=1,nellYL
  nopYL(i,6)=2*nnodesYL+2
ENDDO
DO i=1,nellYL
  nopYL(i,7)=2*nnodesYL+3
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE nodnumbAdaAYL

!*******************************************************************************
!*** DESCRIPTION: NODAL NUMBERING (AUGMENTED YOUNG-LAPLACE WITH ADAPTIVE MESH)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j

DO i=1,nellYL
  DO j=1,2
    nopYL(i,j)=i+(j-1)
  ENDDO
  DO j=3,4
    nopYL(i,j)=nopYL(i,j-2)+nnodesYL
  ENDDO
  DO j=5,6
    nopYL(i,j)=nopYL(i,j-4)+2*nnodesYL
  ENDDO
ENDDO
DO i=1,nellYL
  nopYL(i,7)=3*nnodesYL+1
ENDDO
DO i=1,nellYL
  nopYL(i,8)=3*nnodesYL+2
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE nodnumbADA

!*******************************************************************************
!*** DESCRIPTION: NODAL NUMBERING (ADAPTIVE MESH COMPUTATIONAL DOMAIN)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j

DO i=1,nellADA
  DO j=1,2
    nopADA(i,j)=i+(j-1)
  ENDDO
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE nodnumbCYL

!*******************************************************************************
!*** DESCRIPTION: NODAL NUMBERING (CONVENTIONAL YOUNG-LAPLACE)
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j

DO i=1,nellYL
  DO j=1,2
    nopYL(i,j)=i+(j-1)
  ENDDO
ENDDO
DO i=1,nellYL
  nopYL(i,3)=nnodesYL+1
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE neigElem(nnodes,nellem,nellnodes,nop,kindi,neigh,alloc)

!*******************************************************************************
!*** DESCRIPTION: FIND NEIGHBORING ELEMENTS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nellnodes
INTEGER(KIND=8),INTENT(IN)::alloc
INTEGER(KIND=8),INTENT(IN)::nnodes
INTEGER(KIND=8),INTENT(IN)::nellem
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::nop
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(OUT)::kindi
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:),INTENT(OUT)::neigh
INTEGER(KIND=8)::j
INTEGER(KIND=8)::i

ALLOCATE(kindi(nnodes),neigh(nnodes,alloc))
DO i=1,nellem
  DO j=1,nellnodes
    kindi(nop(i,j))=0
  ENDDO
ENDDO
DO i=1,nellem
  DO j=1,nellnodes
!*** THE TOTAL NUMBER OF NEIGHBORING ELEMENTS OF A NODE
    kindi(nop(i,j))=kindi(nop(i,j))+1
!*** NEIGH(NODE,#NTH)=#NTH NEIGHBORING ELEMENT OF THAT NODE
    neigh(nop(i,j),kindi(nop(i,j)))=i
  ENDDO
ENDDO

END SUBROUTINE

!*******************************************************************************

SUBROUTINE JacobPtrn(nouk,nellnodes,nop,kindi,neigh)

!*******************************************************************************
!*** DESCRIPTION: COMPUTE JACOBIAN MATRIX SPARSITY PATTERN
!*******************************************************************************

USE CommonVars
USE qsort_module
IMPLICIT NONE
INTEGER(KIND=8),INTENT(IN)::nellnodes
INTEGER(KIND=8),INTENT(IN)::nouk
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:),INTENT(IN)::kindi
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::nop
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:,:),INTENT(IN)::neigh
INTEGER(KIND=8)::k
INTEGER(KIND=8)::i
INTEGER(KIND=8)::j
INTEGER(KIND=8)::l
INTEGER(KIND=8)::m
INTEGER(KIND=8)::i1
INTEGER(KIND=8)::i2
INTEGER(KIND=8)::inct2
INTEGER(KIND=8),ALLOCATABLE,DIMENSION(:)::icnt1

ALLOCATE(icnt1(4*noukYL))
DO i=1,nouk
  inct2=0
  DO j=1,kindi(i)
    DO k=1,nellnodes
      inct2=inct2+1
      i1=i
      i2=nop(neigh(i,j),k)
      icnt1(inct2)=(i1-1)*nouk+i2
    ENDDO
  ENDDO
  DO l=1,inct2
    DO m=l,inct2
      IF (icnt1(l)==icnt1(m).AND.l.NE.m) THEN
        icnt1(m)=0
     ENDIF
   ENDDO
  ENDDO
  DO l=1,inct2
    IF (icnt1(l).NE.0) THEN
      AMcnt=AMcnt+1
      am2(AMcnt)=icnt1(l)
    ENDIF
  ENDDO
  DO l=1,SIZE(icnt1)
    icnt1(l)=0
  ENDDO
ENDDO
DO i=1,AMcnt
  am1(i)=am2(i)
ENDDO
DEALLOCATE(am2)
ALLOCATE(am2(AMcnt))
DO i=1,AMcnt
  am2(i)=INT(am1(i),8)
ENDDO
DEALLOCATE(am1,icnt1)
ALLOCATE(am1(AMcnt))
CALL qsort(am2,am1)

END SUBROUTINE

!*******************************************************************************

END MODULE

!*******************************************************************************
