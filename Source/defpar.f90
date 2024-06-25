!*******************************************************************************
!*** DROPlet Simulation (DropS)
!*** BY NIKOLAOS CHAMAKOS (nikoscham@gmail.com)
!*** NATIONAL TECHNICAL UNIVERSITY OF ATHENS, GREECE
!*******************************************************************************

MODULE defpar

CONTAINS

!*******************************************************************************

SUBROUTINE defparam

!*******************************************************************************
!*** DESCRIPTION: GLOBAL PARAMETERS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE

OPEN(10,FILE='GlobalParameters.txt')
READ(10,4010)SolverID !*** SELECT SOLVER
READ(10,4020)JacobMatP !*** JACOBIAN MATRIX PATTERN INPUT
PPres=1.d-6 !*** NEWTON ITERATIONS TOLERANCE
Pi=16.d0*ATAN(1.d0/5.d0)-4.d0*ATAN(1.d0/239.d0)
CLOSE(10)

4010 FORMAT(/,/,i2)
4020 FORMAT(/,/,/,/,/,/,/,/,/,i2)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE defparamYL

!*******************************************************************************
!*** DESCRIPTION: AUGMENTED YOUNG-LAPLACE PARAMETERS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE
DOUBLE PRECISION::thYoungMinYL

OPEN(10,FILE='AYLParameters.txt')
!*** NUMBER OF ELEMENTS
READ(10,5010)nellYL
!*** NUMBER OF ELEMENTS (INITIAL SOLUTION)
READ(10,5020)nelluoYL
!*** VOLUME OF THE DROPLET (m^3)
READ(10,5030)voYL
!*** SOLID-LIQUID INTERACTION CONSTANT
READ(10,5030)sigmaYL
!*** SOLID-LIQUID INTERACTION CONSTANT
READ(10,5030)SmallNumbYL
!*** DISJOINING PRESSURE EXPONENT C1
READ(10,5020)C1YL
!*** DISJOINING PRESSURE EXPONENT C2
READ(10,5020)C2YL
!*** DISJOINING PRESSURE EXPONENT C3
READ(10,5020)C3YL
!*** DISJOINING PRESSURE PRECURSOR FILM PARAMETER
READ(10,5030)flmYL
!*** GRAVITY EFFECT (1=ON, 0=OFF)
GcYL=0.d0
!*** HEIGH LIMIT OF EIKONAL EQUATION
HeighLimitYL=0.14d0
!*** MATERIAL WETTABILITY (YOUNG'S CONTACT ANGLE)
READ(10,5030)thYoungYL
!*** WETTING PARAMETER
CALL wlsCalc (thYoungYL,WParamYL)
!*** WETTING PARAMETER STEP
WParamStepYL=2.d0
!*** MAXIMUM MATERIAL WETTABILITY (MINIMUM YOUNG CONTACT ANGLE)
READ(10,5030)thYoungMinYL
!*** MAXIMUM WETTING PARAMETER
CALL wlsCalc (thYoungMinYL,WParamMaxYL)
IF (WParamYL<0.OR.WParamMaxYL<0) THEN
  WRITE(*,'(a)')"ERROR AT defpar: WLS<0"
  STOP
ENDIF
!*** MINIMUM DROPLET HEIGHT
READ(10,5030)MinDropletHeight
!*** NUMBER OF NODES
nnodesYL=nellYL+1
!*** NUMBER OF NODES (INITIAL SOLUTION)
nnodesuoYL=nelluoYL+1
!*** NUMBER OF UNKNOWNS
noukYL=2*nnodesYL+2
!*** NUMBER OF UNKNOWNS (INITIAL SOLUTION)
noukuoYL=2*nnodesuoYL+2
!*** NUMBER OF UNKNOWNS AT EACH ELEMENT
nelnodesYL=6
!*** INTERFACIAL TENSION (N/m)
READ(10,5030)stYL
!*** CHARACTERISTIC LENGTH
roYL=((3.d0*voYL)/(4.d0*Pi))**(1.d0/3.d0)
!*** FLAT SOLID SURFACE (1=ON, 0=OFF)
READ(10,5020)flatSolidUp
!*** PRINT JACOBIAN MATRIX (1=ON, 0=OFF)
printPrecond=0
!*** ADAPTIVE MESH REFINEMENT (DECOUPLED WITH YL EQUATION) (1=ON, 0=OFF)
READ(10,5020)AdaMesh
!*** PSEUDO-TRANSIENT CONTINUATION
READ(10,5020)TraCon
TStep=5.d-4 !*** PSEUDO-TIME STEP
TStepMax=8.d-4
CLOSE(10)

5010 FORMAT(/,/,i8)
5020 FORMAT(/,i7)
5030 FORMAT(/,f14.11)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE defparamEI

!*******************************************************************************
!*** DESCRIPTION: EIKONAL PARAMETERS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE

OPEN(10,FILE='EIParameters.txt')
!*** GEOMETRY MIRROR FUNCTION (1=ON, 0=OFF)
READ(10,6010)geommirEI
!*** PILLAR WIDTH
READ(10,6030)wpillar
!*** DISTANCE BETWEEN PILLARS
READ(10,6030)dpillar
!*** PILLAR HEIGHT
READ(10,6030)hpillar
!*** PILLAR SMOOTHNESS - SMALLER VALUES CORRESPOND TO SHARPER PILLARS
READ(10,6030)smpillar
!*** REPETITION PERIOD OF THE SOLID SURFACE
scycleEI=wpillar+dpillar
!*** VISCOUS TERM
ViscTermEI=5.d-3
!*** NUMBER OF EXTERNAL ITERATIONS FOR EI SOLUTION
ExtItEI=10
!*** MINIMUM X COORDINATE OF THE 2D DOMAIN
xfirstEI=-1.d-3
!*** MAXIMUM WIDTH OF THE 2D DOMAIN
IF (geommirEI==1) THEN
  xlimitEI=((wpillar+dpillar)/2.d0)
ELSE IF (geommirEI==0) THEN
  xlimitEI=2.d0
ENDIF
!*** MAXIMUM HEIGHT OF THE 2D DOMAIN
ylimitEI=HeighLimitYL+1.d-2
!*** MESH GENERATOR SOFTWARE (1=TRIANGLE, 2=GMSH) 
!*** TRIANGLE INFO: "http://www.cs.cmu.edu/~quake/triangle.html"
!*** GMSH INFO: "http://geuz.org/gmsh/"
READ(10,6020)meshGen
!*** ELEMENT SIZE (DENSE REGION)
READ(10,6030)hdensEI
!*** ELEMENT SIZE (SPARSE REGION)
READ(10,6030)ldensEI
!*** ELEMENT TYPE
!*** (1=LINEAR TRIANGULAR, 2=QUADRATIC TRIANGULAR, 3=QUADRATIC RECTANGULAR)
READ(10,6020)eltypeEI
!*** GAUSS POINTS AND WEIGHTS
CALL gauPointsWeights(eltypeEI,nelnodesEI,gpEI,wgpEI,NumgpsEI)
CLOSE(10)

6010 FORMAT(/,/,i7)
6020 FORMAT(/,i7)
6030 FORMAT(/,f14.11)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE defparamEL

!*******************************************************************************
!*** DESCRIPTION: ELECTRIC FIELD DISTRIBUTION PARAMETERS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE

OPEN(10,FILE='ELParameters.txt')
eoEL=8.854187817d-12
!*** MAXIMUM WIDTH OF THE 2D DOMAIN
READ(10,7010)xlimitEL
!*** MAXIMUM HEIGHT OF THE 2D DOMAIN
READ(10,7030)ylimitEL
!*** MESH GENERATOR SOFTWARE (1=TRIANGLE, 2=GMSH)
READ(10,7020)meshGen
!*** ELEMENT SIZE (DENSE REGION)
READ(10,7030)hdensEL
!*** ELEMENT SIZE (SPARSE REGION)
READ(10,7030)ldensEL
!*** MESH GROWTH RATE
READ(10,7030)GrateEL
!*** SOLID MATERIAL DIELECTRIC CONSTANT
READ(10,7030)edEL
!*** SURROUNDING MEDIUM DIELECTRIC CONSTANT
READ(10,7030)esEL
!*** DIELECTRIC THICKNESS
READ(10,7030)DielThickEL; DielThickEL=(DielThickEL/roYL)
!*** FLAT SOLID SURFACE AT THE BOTTOM OF THE DIELECTRIC (1=ON, 0=OFF)
READ(10,7020)flatSolidDown
!*** APPLIED VOLTAGE
IF (SolverID.NE.2.AND.SolverID.NE.3.AND.SolverID.NE.6) VltEL=0.d0
!*** MAXIMUM APPLIED VOLTAGE
READ(10,7030)VltMaxEL
!*** APPLIED VOLTAGE STEP
READ(10,7030)VltStepEL
!*** ELEMENT TYPE (1=LINEAR TRIANGULAR, 2=QUADRATIC TRIANGULAR)
READ(10,7020)eltypeEL
!*** GAUSS POINTS AND WEIGHTS
CALL gauPointsWeights(eltypeEL,nelnodesEL,gpEL,wgpEL,NumgpsEL)
CLOSE(10)

7010 FORMAT(/,/,f14.11)
7020 FORMAT(/,i7)
7030 FORMAT(/,f14.11)

END SUBROUTINE

!*******************************************************************************

SUBROUTINE defparamCYL

!*******************************************************************************
!*** DESCRIPTION: CONVENTIONAL YOUNG-LAPLACE PARAMETERS
!*******************************************************************************

USE CommonVars
IMPLICIT NONE

OPEN(10,FILE='CYLParameters.txt')
!*** NUMBER OF ELEMENTS
READ(10,8010)nellYL
!*** NUMBER OF ELEMENTS (INITIAL SOLUTION)
READ(10,8020)nelluoYL
!*** VOLUME OF THE DROPLET (m^3)
READ(10,8030)voYL
!*** INTERFACIAL TENSION (N/m)
READ(10,8030)stYL
!*** GRAVITATIONAL ACCELERATION (m/s^2)
READ(10,8030)gYL
!*** LIQUID DENSITY (Kg/m^3)
READ(10,8030)denYL
!*** MATERIAL WETTABILITY (YOUNG CONTACT ANGLE)
READ(10,8030)thYoungYL
thYoungYL=thYoungYL*Pi/180.
!*** NUMBER OF NODES
nnodesYL=nellYL+1
!*** NUMBER OF NODES (INITIAL SOLUTION)
nnodesuoYL=nelluoYL+1
!*** NUMBER OF UNKNOWNS
noukYL=nnodesYL+1
!*** NUMBER OF UNKNOWNS (INITIAL SOLUTION)
noukuoYL=nnodesuoYL+1
!*** NUMBER OF UNKNOWNS AT EACH ELEMENT
nelnodesYL=3
CLOSE(10)
!*** CHARACTERISTIC LENGTH
roYL=((3.d0*voYL)/(4.d0*Pi))**(1.d0/3.d0)

8010 FORMAT(/,/,i7)
8020 FORMAT(/,i7)
8030 FORMAT(/,f14.11)

END SUBROUTINE

!*******************************************************************************

END MODULE

!*******************************************************************************
