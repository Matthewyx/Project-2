!!DIFFUSION PROBLEM WITH LATTICE BOLTZMANN METHOD
!!XU YU (PH.D CANDIDATE), PROF. KLAUS REGENAUER-LIEB, DR. FANGBAO TIAN
!!SCHOOL OF MINERALS AND ENERGY RESOURCES ENGINEERING
!!UNIVERSITY OF NEW SOUTH WALES
!!TIME:2018 AUG.30    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!         MODULES OF PARA DEFIN          !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE LB_mass_transfer
IMPLICIT NONE
!PUBLIC CONSTANTS
INTEGER,PARAMETER::X_MAX    =201            !DOMAIN SIZE IN X AXIS
INTEGER,PARAMETER::X_MIN    =1              !DOMAIN SIZE IN X AXIS
INTEGER,PARAMETER::Y_MAX    =101            !DOMAIN SIZE IN Y AXIS
INTEGER,PARAMETER::Y_MIN    =1              !DOMAIN SIZE IN Y AXIS
INTEGER,PARAMETER::Y_BOUND  =15             !FRACTURE BOTTOM Y POSITION
INTEGER,PARAMETER::DEPTH  	=10             !PORE DEPTH ON THE WALL
INTEGER,PARAMETER::WIDTH  	=10             !PORE WIDTH ON THE WALL
INTEGER,PARAMETER::P_NUM  	=0              !NUMBER OF PORES
INTEGER,PARAMETER::Q_S      =9              !D2Q9 LATTICE MODEL
INTEGER,PARAMETER::T_max    =500000         !MAXIMUM TIME STEP
INTEGER,PARAMETER::T_PLOT   =100          !INTERVAL OF PLOT
REAL(KIND=8),PARAMETER::LX  =2.0            !DOMAIN SIZE IN Y AXIS
REAL(KIND=8),PARAMETER::LY  =1.0            !DOMAIN SIZE IN Y AXIS
REAL(KIND=8),PARAMETER::L_REF   =5.0D-5     !LENGTH REFERENCE
REAL(KIND=8),PARAMETER::U_REF   =0.01       !VELOCITY REFERENCE
REAL(KIND=8),PARAMETER::RHO_REF =0.6        !DENSITY REFERENCE
REAL(KIND=8),PARAMETER::NIU_REF =2.97D-1    !KINETIC VISICOSITY REFERENCE
REAL(KIND=8),PARAMETER::T_REF 	=8.42D-4    !TIME REFERENCE
REAL(KIND=8),PARAMETER::Kr      =1.00D-2    
!!PARAMETER DEFINITION
REAL(KIND=8) G(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S), G_EQ(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S)
REAL(KIND=8) G_POST_COLLID(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S)
REAL(KIND=8) W(Q_S),CX(Q_S),CY(Q_S),OPP(Q_S)
REAL(KIND=8) X(X_MIN:X_MAX,Y_MIN:Y_MAX),Y(X_MIN:X_MAX,Y_MIN:Y_MAX)
REAL(KIND=8) DX,DY,DT,NIU,PI,C0,U_TEM
REAL(KIND=8) C(0:X_MAX,0:Y_MAX)
REAL(KIND=8) TAU1,OMEGA1,D_coef,CK
!
REAL(KIND=8) F(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S),FEQ(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S)
REAL(KIND=8) F_POST(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S)
REAL(KIND=8) RHO(X_MIN:X_MAX,Y_MIN:Y_MAX),U(X_MIN:X_MAX,Y_MIN:Y_MAX),V(X_MIN:X_MAX,Y_MIN:Y_MAX)
!
REAL(KIND=8) TAU,OMEGA,CI,CS,DEN
REAL(KIND=8) U_ANALY(Y_MIN:Y_MAX)
REAL(KIND=8) FIX(X_MIN:X_MAX,Y_MIN:Y_MAX),FIY(X_MIN:X_MAX,Y_MIN:Y_MAX)
REAL(KIND=8) RHO_OUT,RHO_IN,U_MAX,GRADP,PE
REAL(KIND=8) U0_IN,RE
!READ ROUGH SURFACE DATA TO THE MATRIX
REAL(KIND=8) SURF(X_MIN:X_MAX),NUM(X_MIN:X_MAX),Y_B_IN,Y_B_OUT
REAL(KIND=8) N_SORP(X_MIN:X_MAX,Y_MIN:Y_MAX),DS(X_MIN:X_MAX),SS(X_MIN:X_MAX)

REAL(KIND=8) ERROR,TIME
INTEGER(KIND=4) IMAGE(X_MIN:X_MAX,Y_MIN:Y_MAX),IS_WALL(X_MIN:X_MAX,Y_MIN:Y_MAX)

INTEGER BOUND_FLAG
INTEGER(KIND=4) NSTEP,I,J,K
END MODULE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!         MAIN PROGRAM          !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM MAIN
USE LB_mass_transfer
IMPLICIT NONE
REAL(KIND=8) TEM1,TEM2,START,FINISH
CHARACTER*6 TITLE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     IMPORT THE SURFACE DATA        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(12,FILE='sawtooth_h0.1_mesh100.txt')
READ(12,*) TITLE
WRITE(*,*) TITLE
DO I=X_MIN,X_MAX
    READ(12,*) NUM(I),SURF(I)
ENDDO
!!SMOOTH THE SURFACE WITH 10% Neibouring Difference
!DO I=X_MIN,X_MAX-1
!IF(ABS(SURF(I+1)-SURF(I))/SURF(I).GT.0.10)THEN
!    SURF(I+1)=SURF(I)*(1.0+0.18*(SURF(I+1)-SURF(I))/ABS(SURF(I+1)-SURF(I)))
!ENDIF
!ENDDO
!DO I=X_MIN+1,X_MAX-1
!    TEM1=SQRT((SURF(I)-SURF(I-1))**2+DX**2)
!    TEM2=SQRT((SURF(I+1)-SURF(I))**2+DX**2)
!    DS(I)=(TEM1+TEM2)/2.0
!ENDDO
!DS(X_MIN)=SQRT((SURF(X_MIN+1)-SURF(X_MIN))**2+DX**2)
!DS(X_MAX)=SQRT((SURF(X_MAX)-SURF(X_MAX-1))**2+DX**2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!        Definition of params        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!boundary flag: 0-Dirichlet, 1-pressure

BOUND_FLAG  =   1

!!time
CALL CPU_TIME(START)

!!initial
CALL INITIAL

WRITE(*,*) 'THE DOMAIN SIZE:'
WRITE(*,*) 'X_MAX=',X_MAX,";",'Y_MAX=',Y_MAX
WRITE(*,*) 'PARAMETERS OF LATTICE BOLTZMANN METHOD:'
WRITE(*,*) 'Fluid dynamics relaxation time TAU=',TAU
WRITE(*,*) 'Diffusion relaxation time TAU=',TAU1
WRITE(*,*) 'Initial velocity=',U_TEM
WRITE(*,*) 'Peclet number=', PE
WRITE(*,*) 'START_TIME=',START

OPEN(5,FILE='OUTPUT//LOGNOTE.TXT')
WRITE(5,*) 'THE DOMAIN SIZE:'
WRITE(5,*) 'X_MAX=',X_MAX,";",'Y_MAX=',Y_MAX
WRITE(5,*) 'PARAMETERS OF LATTICE BOLTZMANN METHOD:'
WRITE(5,*) 'Fluid dynamics relaxation time TAU=',TAU
WRITE(5,*) 'Diffusion relaxation time TAU=',TAU1
WRITE(5,*) 'Initial velocity=',U_TEM
WRITE(5,*) 'Peclet number=', PE
WRITE(5,*) 'START_TIME=',START
!pause
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!     OUTPUT BREAKTHROUTH CURVE      !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

OPEN(11,FILE='OUTPUT//BREAK_CURVE.PLT',ACCESS='APPEND')
WRITE(11,*) 'VARIABLES=TIME, C/C0'
WRITE(11,*) 'ZONE T='//TITLE
NSTEP=0
CALL OUTPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!            MAIN LOOP START         !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!DO WHILE(NSTEP*DT*D_coef/(Y_MAX*DY)**2.LE.100)
DO WHILE(NSTEP.LE.T_max) 
NSTEP=NSTEP+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    BREAKTHROUGH CUR AND MONITOR    !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    WRITE(11,*) NSTEP, C(X_MAX,Y_MAX)
    WRITE(*, *) NSTEP, C(X_MAX,Y_MAX),U(X_MAX,Y_MAX)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!COMPUTING MACROSCOPIC PARAMETERS !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL MACRO_C
    CALL MACRO_FLOW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! COLLISION_STREAMING_MASS_DIFFUSION  !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL COLLIDE_STREAM_C
    CALL COLLID_STREAM_FLOW
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!     BOUNDARY CONDITIONS         !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CALL BOUND_C
	IF(BOUND_FLAG.EQ.1)THEN
		CALL PRESSURE_BOUND
    ELSEIF(BOUND_FLAG.EQ.0)THEN
		CALL VELOCITY_BOUND
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!    CALCULATION OF  SORPTION AMOUNT     !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO I=X_MIN,X_MAX
DO J=Y_MIN,Y_MAX
    IF(IS_WALL(I,J).EQ.1)THEN
    N_SORP(I,J)	=	N_SORP(I,J)+SS(I)*DT
    ENDIF
ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!       OUTPUT FILES/T_PLOT TIME        !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	IF(MOD(NSTEP,T_PLOT).EQ.0)THEN
	    CALL OUTPUT
	ENDIF
ENDDO
CALL CPU_TIME(FINISH)
WRITE(5,*) 'END_TIME=',FINISH
WRITE(5,*) 'TOTAL_TIME=',FINISH-START
END PROGRAM MAIN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SEGMENTATION OF THE DOMAIN !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
SUBROUTINE SEGMENTATION
USE LB_mass_transfer
!SOLID=0, FLUID=1, INTERFACE
!L_INT	=	(X_MAX-1)/(P_NUM+1)		! PORE INTERVAL
IMAGE(:,:)  =   0
IS_WALL(:,:)=   0
DO I=X_MIN, X_MAX
DO J=Y_MIN, Y_MAX
IF(Y(I,J).GT.SURF(I))THEN
    IMAGE(I,J)=1
ENDIF
ENDDO
ENDDO
!DO I=1,P_NUM
!    T1	=	L_INT*I
!    T2	=	L_INT*I+WIDTH
!    IMAGE(T1:T2,Y_BOUND-DEPTH:Y_BOUND)=   1
!    IS_WALL(T1:T2, Y_BOUND-DEPTH:Y_BOUND)=   1
!ENDDO
! CHECK THE NEIGHBOURS OF SOLID NODES
DO I    =   X_MIN,X_MAX
DO J    =   Y_MIN+1,Y_MAX-1
IF(IMAGE(I,J).EQ.1)THEN
!NORTH(N),SOUTH(S),WEST(W),EAST(E)
IF(IMAGE(I,J+1).EQ.0)THEN 		!N
	IS_WALL(I,J)   =   1
ELSEIF(IMAGE(I,J-1).EQ.0)THEN	!S
	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I-1,J).EQ.1)THEN	!W
!	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I+1,J).EQ.1)THEN	!E
!	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I+1,J+1).EQ.1)THEN	!NE
!	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I-1,J+1).EQ.1)THEN	!NW
!	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I-1,J-1).EQ.1)THEN	!SW
!	IS_WALL(I,J)   =   1
!ELSEIF(IMAGE(I+1,J-1).EQ.1)THEN	!SE
!	IS_WALL(I,J)   =   1
ENDIF
ENDIF
ENDDO
ENDDO
IS_WALL(X_MIN,CEILING(SURF(X_MIN)/DX)) =   1
IS_WALL(X_MAX,CEILING(SURF(X_MAX)/DX)) =   1
!!================================================================!!
!!POSITIONS OF BOUNDARY NODES
OPEN(1,FILE='SURF_DATA.PLT')
DO I=X_MIN,X_MAX
WRITE(1,*) (I-1)*DX, SURF(I)
ENDDO
CLOSE(1)
OPEN(1,FILE='OUTPUT//IMAGE.PLT')
WRITE(1,*)'VARIABLES=X,Y,IMG,IS_WALL'
WRITE(1,*) 'ZONE T="ZONE 1"'
WRITE(1,*) 'I=',X_MAX, ' J=',Y_MAX
WRITE(1,*) 'F=POINT'
DO J=Y_MIN,Y_MAX
DO I=X_MIN,X_MAX
    WRITE(1,*) X(I,J),Y(I,J),IMAGE(I,J),IS_WALL(I,J)
ENDDO
ENDDO
CLOSE(1)
RETURN
    END SUBROUTINE SEGMENTATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!============    INPUT PARAMETERS    ===============!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INPUT
USE LB_mass_transfer
REAL(KIND=8) L_R

    L_R     =   LY-SURF(X_MIN)          !REFERENCE LENGTH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    INPUT THE DIFFUSION PARAMS  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    D_COEF  =   1.0D-4                   !! mass diffusion coeficient
    TAU1    =   D_COEF/CS**2/DT+0.5		!! relaxation time
    OMEGA1  =   1.0/TAU1                !! relaxation parameter
    C0      =   1.0                     !! INITIAL CONCENTRATION
    DEN     =   1.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! THE FOLLOWING PROCESS IS TO INPUT THE VALUES FOR THE FLUID DYNAMICS
!!IF THE FLAG OF BOUNDARY CONDITION IS, O FOR THE VELOCITY BOUND
!! 1 FOR PRESSURE BOUNDARY CONDITION
IF(BOUND_FLAG.EQ.0)THEN     !!VELOCITY BOUND INITIAL     !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    RE  	=   2.0D0                   !Reynolds NUMBER
    U0_IN   =   0.1                    !INLET VELOCITY
    PE      =   U0_IN*L_R/D_COEF        !PECLET NUMBER
    NIU		=	U0_IN*L_R/RE            !KINEMATIC NUMBER
    TAU		=	NIU/DT/CS**2+0.5        !RELAXATION TIME
    OMEGA   =   1.0/TAU 
    U_TEM   =   U0_IN

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!     	INITIAL STATUS      !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    U(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   U_TEM*IMAGE(X_MIN:X_MAX,Y_MIN:Y_MAX)
    V(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   0.0
    RHO(X_MIN:X_MAX,Y_MIN:Y_MAX)=   DEN
    C(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   0.0
    DO J=Y_MIN,Y_MAX
        C(X_MIN,J)  =   C0*IMAGE(X_MIN,J)
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELSEIF(BOUND_FLAG.EQ.1)THEN     !! PRESSURE BOUND INITIAL    !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    TAU     =   0.8
    NIU     =   (TAU-0.5)*DT*CS**2
    OMEGA   =   1.0/TAU
    !！
    PE      =   1000                  !PECLET NUMBER
    U_MAX   =   PE*D_COEF/L_R
    GRADP   =   8.0*DEN*NIU*U_MAX/(4.0*L_R**2)   !1.0D-4, PRESSURE GRADIENT, LB UNIT
    RHO_OUT =   DEN                             !-0.5*GRADP*LX/CS**2
    RHO_IN  =   DEN+GRADP*LX/CS**2
    !U_MAX  =   4.0*LY**2*GRADP/(8.0*NIU*DEN)
    !PE     =   U_MAX*2.0*LY/D_COEF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!     ANALYTICAL SOLUTION     !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO J=Y_MIN,Y_MAX
        U_ANALY(J)  = 4.0*U_MAX*(2.0*L_R*Y(1,J)-Y(1,J)*Y(1,J))/L_R**2
    ENDDO
        U_TEM   =   U_MAX
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!     	INITIAL STATUS      !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    U(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   0.0
    V(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   0.0
    RHO(X_MIN:X_MAX,Y_MIN:Y_MAX)=   DEN
    C(X_MIN:X_MAX,Y_MIN:Y_MAX)  =   0.0
ENDIF

RETURN
END SUBROUTINE INPUT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!        INITIALISATION	PARAMETERS      !!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE INITIAL
USE LB_mass_transfer
IMPLICIT NONE
REAL(KIND=8) TEMP1,CU,U2,T1,T2
INTEGER(KIND=8) L_INT		!INTERVAL BETWEEN TWO PORE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!        PUBLIC PARAMETERS        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DX      =   LX/(X_MAX-1)            !X SPACE STEP OF CARTISIAN COORDINATE
DY      =   DX                      !Y SPACE STEP OF CARTISIAN COORDINATE
DT      =   DX                      !TIME STEP LENGTH
CK      =   DX/DT                   !lattice speed
CS      =   1.0/SQRT(3.0)           !SPEED OF SOUND
PI      =   4.0*ATAN(1.0)           ! pi
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MESH GRID GENERATOR        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO I=X_MIN,X_MAX
DO J=Y_MIN,Y_MAX
    X(I,J)  = (I-1)*DX
    Y(I,J)  = (J-1)*DY
ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  LATTICE SPEED&WEIGHT FACTOR   !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
DO I=1,4
CX(I)   =  COS((I-1)*0.5*PI)*DX/DT
CY(I)   =  SIN((I-1)*0.5*PI)*DY/DT
W(I)    =   1./9.
ENDDO
DO I=5,8
CX(I)   =  SQRT(2.0)*COS(0.25*PI+0.5*PI*(I-5))*DX/DT
CY(I)   =  SQRT(2.0)*SIN(0.25*PI+0.5*PI*(I-5))*DY/DT
W(I)    =   1./36.
ENDDO
CX(9)   =  0.0
CY(9)   =  0.0
W(9)    =  4./9.
OPP(:)  = (/3,4,1,2,7,8,5,6,9/)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  SEGMENTATION OF THE DOMAIN !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL SEGMENTATION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!  INPUT SIMULATION PARAMETERS !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

CALL INPUT

!!======================================================================!!
!!THE INITIAL STATE OF THE DISTRIBUTION FUCNTION F &G= EQUILIBRIUM DISTRIBUTION FUNCTION
DO I=X_MIN,X_MAX
DO J=Y_MIN,Y_MAX
DO K=1,Q_S
    CU  = CX(K)*U(I,J)+CY(K)*V(I,J)
    U2  = U(I,J)**2+V(I,J)**2
    FEQ(I,J,K)  = W(K)*RHO(I,J)*(1.0+3.0*CU+4.5*CU**2-1.5*U2)     
ENDDO
ENDDO
ENDDO
F   =  FEQ
DO I=X_MIN,X_MAX
DO J=Y_MIN,Y_MAX
DO K=1,Q_S
    TEMP1=CX(K)*U(I,J)+CY(K)*V(I,J)
	G_EQ(I,J,K)=W(K)*C(I,J)*(1+3.0*TEMP1)
ENDDO
ENDDO
ENDDO
G   =  G_EQ

RETURN
END SUBROUTINE INITIAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     	    MACROSCOPIC PARAMETERS              !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE MACRO_C
    USE LB_mass_transfer
    IMPLICIT NONE
    REAL(KIND=8) CC
    DO I=X_MIN,X_MAX
    DO J=Y_MIN,Y_MAX
    !IF(IMAGE(I,J).NE.0)THEN
        CC=0.0
        DO K=1,Q_S
            CC=CC+G(I,J,K)
        ENDDO
        C(I,J)=CC
    !ELSE
    IF(image(i,j).eq.0)then
        C(I,J)=0.0
    ENDIF
    ENDDO
    ENDDO
    RETURN
END SUBROUTINE MACRO_C

SUBROUTINE MACRO_FLOW
    USE LB_mass_transfer
    IMPLICIT NONE
    REAL(KIND=8) DD,UU,VV
    DO I=X_MIN,X_MAX
    DO J=Y_MIN,Y_MAX
    IF(IMAGE(I,J).NE.0)THEN
        DD=0.0
        UU=0.0
        VV=0.0
        DO K=1,Q_S
            DD=DD+F(I,J,K)
            UU=UU+F(I,J,K)*CX(K)
            VV=VV+F(I,J,K)*CY(K)
        ENDDO
        RHO(I,J)= DD
        U(I,J)  = UU/DD 
        V(I,J)  = VV/DD 
    ELSE
        U(I,J)  =   0.0
        V(I,J)  =   0.0
        RHO(I,J)=   DEN
    ENDIF
    ENDDO
    ENDDO
    RETURN
END SUBROUTINE MACRO_FLOW    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!	    STREMING & COLLIDING          !!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!For fluid dynamics,   COLLID-STREAM PROCESS
SUBROUTINE COLLID_STREAM_FLOW
USE LB_mass_transfer
IMPLICIT NONE
REAL(KIND=8) CU,U2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      COLLISION PROCESS       !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO I=   X_MIN, X_MAX
    DO J=   Y_MIN, Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    DO K=   1, Q_S
        CU  =   CX(K)*U(I,J)+CY(K)*V(I,J)
        U2  =   U(I,J)**2+V(I,J)**2
        FEQ(I,J,K)  =   W(K)*RHO(I,J)*(1.0+3.0*CU+4.5*CU*CU-1.5*U2)     
        F(I,J,K)    =   (1-OMEGA)*F(I,J,K)+OMEGA*FEQ(I,J,K)
        ENDDO
        ENDIF
    ENDDO
    ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!      STREAMING PROCESS       !!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    F_POST  =   F
    do I=   X_MIN+1,X_MAX
    do J=   Y_MIN,Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,1)=   F_POST(I-1,J,1)
        ENDIF
    enddo
    enddo
    do I=   X_MIN,X_MAX
    do J=   Y_MIN+1,Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,2)=   F_POST(I,J-1,2)
        ENDIF
    enddo
    enddo
    do I=   X_MIN,X_MAX-1
    do J=   Y_MIN,Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,3)=   F_POST(I+1,J,3)
        ENDIF
    enddo
    enddo
    do I=   X_MIN,X_MAX
    do J=   Y_MIN,Y_MAX-1
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,4)=   F_POST(I,J+1,4)
        ENDIF
    enddo
    enddo
    do I=   X_MIN+1,X_MAX
    do J=   Y_MIN+1,Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,5)=   F_POST(I-1,J-1,5)
        ENDIF
    enddo
    enddo
    do I=   X_MIN,X_MAX-1
    do J=   Y_MIN+1,Y_MAX
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,6)=   F_POST(I+1,J-1,6)
        ENDIF
    enddo
    enddo
    do I=   X_MIN,X_MAX-1
    do J=   Y_MIN,Y_MAX-1
        IF(IMAGE(I,J).NE.0)THEN
	    F(I,J,7)=   F_POST(I+1,J+1,7)
        ENDIF
    enddo
    enddo
    do I=   X_MIN+1,X_MAX
    do J=   Y_MIN,Y_MAX-1
        IF(IMAGE(I,J).NE.0)THEN
	    f(I,J,8)=   f_post(I-1,J+1,8)
        ENDIF
    enddo
    enddo
    RETURN
END SUBROUTINE COLLID_STREAM_FLOW
!!COLLISION AND STREAMING PROCESSES ARE FUNCTIONING ABOUT THE Y_BOUND
SUBROUTINE COLLIDE_STREAM_C
USE LB_mass_transfer
IMPLICIT NONE
REAL GG(X_MIN:X_MAX,Y_MIN:Y_MAX,Q_S), TEMP1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!      COLLID PROCESSING       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO I=   X_MIN, X_MAX
    DO J=   Y_MIN, Y_MAX
    !IF(IMAGE(I,J).NE.0)THEN
    DO K =   1, Q_S
        TEMP1=CX(K)*U(I,J)+CY(K)*V(I,J)
        G_EQ(I,J,K)=W(K)*C(I,J)*(1+3.0*TEMP1)
        G(I,J,K)=(1-OMEGA1)*G(I,J,K)+OMEGA1*G_EQ(I,J,K)
    ENDDO
    !ENDIF
    ENDDO
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!      STREAM PROCESSING       !!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    GG  =   G
    DO I = X_MIN+1, X_MAX
    DO J = Y_MIN, Y_MAX
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,1)=GG(I-1,J,1)
        !ENDIF
    ENDDO
    ENDDO
    DO I = X_MIN, X_MAX
    DO J = Y_MIN+1, Y_MAX
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,2)=GG(I,J-1,2)
        !ENDIF
    ENDDO
    ENDDO
    DO I =  X_MIN, X_MAX-1
    DO J =  Y_MIN, Y_MAX
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,3)=GG(I+1,J,3)
        !ENDIF
    ENDDO
    ENDDO
    DO I =   X_MIN, X_MAX
    DO J =  Y_MIN, Y_MAX-1
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,4)=GG(I,J+1,4)
        !ENDIF
    ENDDO
    ENDDO
    DO I =   X_MIN+1, X_MAX
    DO J =   Y_MIN+1, Y_MAX
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,5)=GG(I-1,J-1,5)
        !ENDIF
    ENDDO
    ENDDO
    DO I =   X_MIN, X_MAX-1
    DO J =   Y_MIN+1, Y_MAX
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,6)=GG(I+1,J-1,6)
        !ENDIF
    ENDDO
    ENDDO
    DO I=   X_MIN, X_MAX-1
    DO J=   Y_MIN, Y_MAX-1
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,7)=GG(I+1,J+1,7)
        !ENDIF
    ENDDO
    ENDDO
    DO I =   X_MIN+1,X_MAX
    DO J =   Y_MIN,Y_MAX-1
        !IF(IMAGE(I,J).NE.0)THEN
        G(I,J,8)=GG(I-1,J+1,8)
        !ENDIF
    ENDDO
    ENDDO
    RETURN 
END SUBROUTINE COLLIDE_STREAM_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     PRESSURE BOUNDARY CONDITION     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE PRESSURE_BOUND
USE LB_mass_transfer
IMPLICIT NONE
!FOR IS_WALL=1 on the WALL
!WEST: INLET DENSITY=RHO_IN
DO J = Y_MIN,Y_MAX
    IF(IMAGE(X_MIN,J).NE.0)THEN
    U(X_MIN,J)  =   1.0-(F(X_MIN,J,2)+F(X_MIN,J,4)+F(X_MIN,J,9)+2.0*(F(X_MIN,J,3) &
                    +F(X_MIN,J,6)+F(X_MIN,J,7)))/RHO_IN
    V(X_MIN,J)  =   0.0
    F(X_MIN,J,1)=   F(X_MIN,J,3)+(2.0/3.0)*RHO_IN*U(X_MIN,J)
    F(X_MIN,J,5)=   F(X_MIN,J,7)-0.5*(F(X_MIN,J,2)-F(X_MIN,J,4))+RHO_IN*U(X_MIN,J)/6.0
    F(X_MIN,J,8)=   F(X_MIN,J,6)+0.5*(F(X_MIN,J,2)-F(X_MIN,J,4))+RHO_IN*U(X_MIN,J)/6.0
    ENDIF
!EAST: DENSITY=RHO_OUT, GRANDP=8*NIU*U0_IN/LY^2
    IF(IMAGE(X_MAX,J).NE.0)THEN
    U(X_MAX,J)  =   -1.0+(F(X_MAX,J,2)+F(X_MAX,J,4)+F(X_MAX,J,9)+2.0*(F(X_MAX,J,1) &
                    +F(X_MAX,J,5)+F(X_MAX,J,8)))/RHO_OUT
    V(X_MAX,J)  =   0.0
    F(X_MAX,J,3)=   F(X_MAX,J,1)-(2.0/3.0)*RHO_OUT*U(X_MAX,J)
    F(X_MAX,J,6)=   F(X_MAX,J,8)-0.5*(F(X_MAX,J,2)-F(X_MAX,J,4))-RHO_OUT*U(X_MAX,J)/6.0
    F(X_MAX,J,7)=   F(X_MAX,J,5)+0.5*(F(X_MAX,J,2)-F(X_MAX,J,4))-RHO_OUT*U(X_MAX,J)/6.0
    ENDIF
ENDDO
!NORTH:
DO K = 1, Q_S
DO I = X_MIN, X_MAX
    F(I,Y_MAX,K)    =   F(I,Y_MAX-1,K)
ENDDO
ENDDO
!SOUTH:
DO K = 1, Q_S
DO I =  X_MIN,X_MAX
DO J =  Y_MIN,Y_MAX
    IF(IS_WALL(I,J).EQ.1)THEN
        F(I,J,K)   =   F(I,J,OPP(K))
    ENDIF
!IF(IMAGE(I,J).EQ.0)THEN
!	IF(IMAGE(I,J+1).EQ.1)THEN 		!N
!		G(I,J,K)  	=	G(I,J,OPP(K))    
!    ELSEIF(IMAGE(I,J-1).EQ.1)THEN	!S
!		G(I,J,K)  	=	G(I,J,OPP(K))
!	ELSEIF(IMAGE(I-1,J).EQ.1)THEN	!W
!		G(I,J,K)  	=	G(I,J,OPP(K))
!	ELSEIF(IMAGE(I+1,J).EQ.1)THEN	!E
!		G(I,J,K)  	=	G(I,J,OPP(K))
	!ELSEIF(IMAGE(I+1,J+1).EQ.1)THEN	!NE
	!	G(I,J,K)  	=	G(I+1,J+1,K)
	!ELSEIF(IMAGE(I-1,J+1).EQ.1)THEN	!NW
	!	G(I,J,K)  	=	G(I-1,J+1,K)
	!ELSEIF(IMAGE(I-1,J-1).EQ.1)THEN	!SW
	!	G(I,J,K)  	=	G(I-1,J-1,K)
	!ELSEIF(IMAGE(I+1,J-1).EQ.1)THEN	!SE
	!	G(I,J,K)  	=	G(I+1,J-1,K)
!	ENDIF
!ENDIF
ENDDO
ENDDO
ENDDO
RETURN    
END SUBROUTINE PRESSURE_BOUND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     PERIODIC BOUNDARY CONDITION     !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE VELOCITY_BOUND
USE LB_mass_transfer
IMPLICIT NONE
!WEST&EAST
U(X_MIN,Y_MIN:Y_MAX)  =   U0_IN*IMAGE(X_MIN,Y_MIN:Y_MAX)
V(X_MIN,Y_MIN:Y_MAX)  =   0.0	
RHO(X_MIN,Y_MIN:Y_MAX)=   RHO(X_MIN+1,Y_MIN:Y_MAX)
U(X_MAX,Y_MIN:Y_MAX)  =   U0_IN*IMAGE(X_Max,Y_MIN:Y_MAX)
V(X_MAX,Y_MIN:Y_MAX)  =   0.0	
RHO(X_MAX,Y_MIN:Y_MAX)=   RHO(X_MAX-1,Y_MIN:Y_MAX)
DO K=1,Q_S
DO J=Y_MIN,Y_MAX
    !IF(IMAGE(X_MIN,J).EQ.1)THEN
    F(X_MIN,J,K)    =	FEQ(X_MIN,J,K)
    !ENDIF
    !IF(IMAGE(X_MAX,J).EQ.1)THEN
    F(X_MAX,J,K)    =	FEQ(X_Max,J,K)!F(X_MAX-1,J,K)
    !ENDIF
ENDDO
ENDDO

!NORTH&SOUTH
DO K=1,Q_S
DO I=X_MIN,X_MAX
    F(I,Y_MAX,K)=	F(I,Y_MAX-1,K)
DO J=Y_MIN,Y_MAX
    IF(IS_WALL(I,J).EQ.1)THEN
    F(I,J,K)	=	F(I,J,OPP(K))
    ENDIF
ENDDO
ENDDO 
ENDDO
!DENSITIS AT FOUR CORNERS Y_BOUND
!RHO(X_MIN,Y_BOUND)  =   2.0*RHO(X_MIN+1,Y_BOUND+1)-RHO(X_MIN+2,Y_BOUND+2)
!RHO(X_MAX,Y_BOUND)  =   2.0*RHO(X_MAX-1,Y_BOUND+1)-RHO(X_MAX-2,Y_BOUND+2)
!RHO(X_MIN,Y_MAX)  =   2.0*RHO(X_MIN+1,Y_MAX-1)-RHO(X_MIN+2,Y_MAX-2)
!RHO(X_MAX,Y_MAX)  =   2.0*RHO(X_MAX-1,Y_MAX-1)-RHO(X_MAX-2,Y_MAX-2)
!!CONDITION: INLET FLUX OR VELOCITY IS KNOWN
	!WEST:: INLET VELOCITY (U0_IN,0)
!    DO J=Y_BOUND+1,Y_MAX
!    RHO(X_MIN,J)= (F(X_MIN,J,9)+F(X_MIN,J,2)+F(X_MIN,J,4)+2.0*(F(X_MIN,J,3)+F(X_MIN,J,6)+F(X_MIN,J,7)))/(1-U0_IN)    
!   F(X_MIN,J,1)= F(X_MIN,J,3)+2.0*RHO(X_MIN,J)*U0_IN/3.0
!    F(X_MIN,J,5)= F(X_MIN,J,7)-0.5*(F(X_MIN,J,2)-F(X_MIN,J,4))+RHO(X_MIN,J)*U0_IN/6.0
!	F(X_MIN,J,8)= F(X_MIN,J,6)+0.5*(F(X_MIN,J,2)-F(X_MIN,J,4))+RHO(X_MIN,J)*U0_IN/6.0
    !EAST:: OUTLET VELOCITY U(xmax-1)=U(xmax-1)
!    RHO(X_MAX,J)=F(X_MAX,J,9)+F(X_MAX,J,2)+F(X_MAX,J,4)+2.0*(F(X_MAX,J,1)+F(X_MAX,J,5)+F(X_MAX,J,8))
!    F(X_MAX,J,3)=F(X_MAX,J,1)
!    F(X_MAX,J,6)=F(X_MAX,J,8)-0.5*(F(X_MAX,J,2)-F(X_MAX,J,4))
!    F(X_MAX,J,7)=F(X_MAX,J,5)+0.5*(F(X_MAX,J,2)-F(X_MAX,J,4))
!    ENDDO
    !NORTH: 
!    DO I=X_MIN,X_MAX
!    !RHO(I,Y_MAX)=F(I,Y_MAX,1)+F(I,Y_MAX,3)+F(I,Y_MAX,9)+2.0*(F(I,Y_MAX,2)+F(I,Y_MAX,5)+F(I,Y_MAX,6))
!	F(I,Y_MAX,4)=F(I,Y_MAX-1,4)!F(I,Y_MAX,2)
!    F(I,Y_MAX,7)=F(I,Y_MAX-1,7)!F(I,Y_MAX,5)+0.5*(F(I,Y_MAX,1)-F(I,Y_MAX,3))
!    F(I,Y_MAX,8)=F(I,Y_MAX-1,8)!F(I,Y_MAX,6)-0.5*(F(I,Y_MAX,1)-F(I,Y_MAX,3)
!    ENDDO
    !SOUTH VELOCITY (0,0)
!    DO K=1,Q_S
!    DO I=X_MIN,X_MAX
!    DO J=Y_MIN,Y_MAX
!    IF(IS_WALL(I,J).EQ.1)THEN
!	F(I,J,K)=F(I,J,OPP(K))
!    ENDIF
!    ENDDO
!    ENDDO
!    ENDDO
RETURN
END SUBROUTINE VELOCITY_BOUND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
!!!!!!!!!	     C: BOUNDARY CONDITIONS	    !!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BOUND_C
USE LB_mass_transfer
IMPLICIT NONE
!WEST: INLET CONCENTRATION C=1.0 / a fixed amount of gas
DO K=1,Q_S
DO J=Y_MIN,Y_MAX
!IF(IMAGE(X_MIN,J).EQ.1)THEN
    G(X_MIN,J,K)    =   W(K)*C0+W(OPP(K))*C0-G(X_MIN,J,OPP(K))
!ENDIF
ENDDO
ENDDO
!EAST: DC/DX=0
DO K=1,Q_S
DO J=Y_MIN,Y_MAX
!IF(IMAGE(X_MAX,J).NE.0)THEN
    G(X_MAX,J,K)    =   G(X_MAX-1,J,K)
!ENDIF
ENDDO
ENDDO
!NORTH: DC/DX=0
DO K=1,Q_S
DO I=X_MIN,X_MAX
    G(I,Y_MAX,K)    =   G(I,Y_MAX-1,K)
ENDDO
!SOUTH! ADIABATICS dC/dy=0
DO J =  Y_MIN, Y_MAX-1
DO I =  X_MIN, X_MAX
!IF(IS_WALL(I,J).EQ.1)THEN
    SS(I)   =   0.0  !Kr*(1.0-N_SORP(I,J))*C(I,Y_MAX)
    G(I,J,K)=   G(I,J+1,K)  !-SS(I)*DX/D_COEF
!ENDIF
!IF(IMAGE(I,J).EQ.0)THEN
!NORTH(N),SOUTH(S),WEST(W),EAST(E)
!	IF(IMAGE(I,J+1).EQ.1)THEN 		!N
!		G(I,J,K)  	=	G(I,J+1,K)    
!    ELSEIF(IMAGE(I,J-1).EQ.1)THEN	!S
!		G(I,J,K)  	=	G(I,J-1,K)
!	ELSEIF(IMAGE(I-1,J).EQ.1)THEN	!W
!		G(I,J,K)  	=	G(I-1,J,K)
!	ELSEIF(IMAGE(I+1,J).EQ.1)THEN	!E
!		G(I,J,K)  	=	G(I+1,J,K)
	!ELSEIF(IMAGE(I+1,J+1).EQ.1)THEN	!NE
	!	G(I,J,K)  	=	G(I+1,J+1,K)
	!ELSEIF(IMAGE(I-1,J+1).EQ.1)THEN	!NW
	!	G(I,J,K)  	=	G(I-1,J+1,K)
	!ELSEIF(IMAGE(I-1,J-1).EQ.1)THEN	!SW
	!	G(I,J,K)  	=	G(I-1,J-1,K)
	!ELSEIF(IMAGE(I+1,J-1).EQ.1)THEN	!SE
	!	G(I,J,K)  	=	G(I+1,J-1,K)
!	ENDIF
!ENDIF
ENDDO
ENDDO
ENDDO
RETURN
END SUBROUTINE BOUND_C
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!        OUTPUT  TO TECPLOT  FILES      !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE OUTPUT
USE LB_mass_transfer
IMPLICIT NONE
INTEGER I1,I2,J1,J2,M,N,I0,I00,I3,I4,I5,ZZ
REAL(KIND=8) UV(X_MIN:X_MAX,Y_MIN:Y_MAX)
CHARACTER*8 NAME

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!        NUMBER GENERATOR     !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ZZ=ICHAR('0')
I00=NSTEP/1000000
I0=NSTEP/100000-I00*10
I1=NSTEP/10000-I00*100-I0*10
I2=NSTEP/1000-I00*1000-I0*100-I1*10
I3=NSTEP/100-I00*10000-I0*1000-I1*100-I2*10
I4=NSTEP/10-I00*100000-I0*10000-I1*1000-I2*100-I3*10
I5=NSTEP/1-I00*1000000-I0*100000-I1*10000-I2*1000-I3*100-I4*10

I00=ZZ+I00
I0=ZZ+I0
I1=ZZ+I1
I2=ZZ+I2
I3=ZZ+I3
I4=ZZ+I4
I5=ZZ+I5
DO I=X_MIN,X_MAX
DO J=Y_MIN,Y_MAX
    UV(I,J)=SQRT(U(I,J)**2+V(I,J)**2)
ENDDO
ENDDO

NAME=CHAR(I00)//CHAR(I0)//CHAR(I1)//CHAR(I2)//CHAR(I3)//CHAR(I4)//CHAR(I5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!        MID PLANE PROFILE    !!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(1,FILE='OUTPUT//XPROF'//NAME//'.PLT')
WRITE(1,*)'VARIABLES=X,C,DEN'
DO I=X_MIN,X_MAX
    WRITE(1,*) X(I,Y_MAX),C(I,Y_MAX),RHO(I,Y_MAX)
ENDDO
CLOSE(1)
OPEN(1,FILE='OUTPUT//YPROF'//NAME//'.PLT')
WRITE(1,*)'VARIABLES=X,U,C'
DO J=Y_MIN,Y_MAX
    WRITE(1,*) Y(X_MAX/2,J),U(X_MAX/2,J),C(X_MAX/2,J)
ENDDO
CLOSE(1)
OPEN(1,FILE='OUTPUT//FIELD'//NAME//'.PLT')
WRITE(1,*)'VARIABLES='
WRITE(1,*)'"X" "Y" "C" "U" "V" "D" "N/N_sat"'
WRITE(1,*) 'ZONE T="ZONE 1"'
WRITE(1,*) 'I=',X_MAX, ' J=',Y_MAX
WRITE(1,*) 'F=POINT'
DO J=Y_MIN,Y_MAX
DO I=X_MIN,X_MAX
	WRITE(1,*) X(I,J),Y(I,J),C(I,J),U(I,J),V(I,J),RHO(I,J),N_SORP(I,J)
ENDDO
ENDDO
CLOSE(1)
RETURN
END SUBROUTINE   