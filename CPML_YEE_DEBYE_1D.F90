PROGRAM ADE_FDTD_ON_YEE_CPML                                                                                   
!-----------------------------------------------------------                                                         
! CPML + Debye @1D TM mode Maxwell's equations with ADE-FDTD on Yee's grid system                                                                                                       
! Yu-Wei Chang    
!-----------------------------------------------------------                                                         
    IMPLICIT NONE                                                                                                      
    INTEGER, PARAMETER::	DP = SELECTED_REAL_KIND(P=15)                                          
    
    ! Computational Doman Parameters
    INTEGER::	I, Isource, Nx, N, Ier, T, numTimeSteps  
    INTEGER::	numGridElectric, numGridMagnetic
    REAL(KIND=DP):: Xstart, Xend, Tstart, Tend, DeltaX, DeltaT, Cr
    REAL(KIND=DP)::	D_a, D_b    !coefficient parameter in the field updates
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::  C_a, C_b, Beta_d, k_d	!coefficient array in the field updates			
  
    !Physical Constants
    REAL(KIND=DP):: Mu_0, Epsilon_0, Tn, PI, C, Eta_0 
  
    !Free Space Parameters
    REAL(KIND=DP)::	Mu_freespace, Epsilon_freespace
  
    !Debye Parameters
    REAL(KIND=DP)::	DeltaEpsilon_Debye, Tau_Debye, Epsilon_infinite_Debye, SigmaS_Debye
  
    !Physical Variables: Electric Field, Magnetic Field, Polarization Current
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::	Ez, Hy, Ez_X, Hy_X, J_d, Ez_asterisk	
  
    !Material Variables
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::	Epsilon_r, Sigma, Tau, DeltaEpsilon
  
    !CPML Parameters
    INTEGER::   m, ma, numCPML
    REAL(KIND=DP)::	Sigma_xmax_freeSpace, Sigma_xmax_debye, Kappa_xmax, a_xmax
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::	Psi_Hyx_NegativeDir, Psi_Hyx_PositiveDir, 	&   ! for the magnetic field updates
                                                bhx_negativeDir, bhx_positiveDir, chx_negativeDir, chx_positiveDir, &	
                                                Sigma_xh_negativeDir, Sigma_xh_positiveDir, &   
  												Kappa_xh_negativeDir, Kappa_xh_positiveDir, &   
  												a_xh_negativeDir, a_xh_positiveDir, &   
  												Psi_Ezx_NegativeDir, Psi_Ezx_PositiveDir, 	& ! for the electric field updates
  												bex_negativeDir, bex_positiveDir, cex_negativeDir, cex_positiveDir, &	
  												Sigma_xe_negativeDir, Sigma_xe_positiveDir, &   
  												Kappa_xe_negativeDir, Kappa_xe_positiveDir, &   
  												a_xe_negativeDir, a_xe_positiveDir  
  																							
    REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::  X                                                                      
    CHARACTER::   PARTNAME*10, STEPCOUNTER*6                                                           
  
    !--------------------------------------------
    !Set General Parameters
    PI=   4.0_DP * ATAN(1.0_DP)                   
    C=    299792458.0_DP                           
    Mu_0= 4.0_DP*PI*1.0E-7_DP                  
    Epsilon_0=    1.0_DP / (Mu_0*C**2)   
    Eta_0=    SQRT(Mu_0/Epsilon_0) 
  
    !--------------------------------------------
    !Set Vacuum Parameters
    Epsilon_freespace=  1.0_DP  
    Mu_freespace= 1.0_DP                     
  
    !--------------------------------------------
    !Set Debye Parameters
    DeltaEpsilon_Debye= 3.0_DP
    Tau_Debye=  7.0_DP*1.0E-12_DP
    Epsilon_infinite_Debye= 7.0_DP
    SigmaS_Debye=   0.0_DP
  
    !--------------------------------------------
    !Set Computational Environment
    DeltaX= 75.0_DP*1.0E-6_DP
    DeltaT= 0.125_DP*1.0E-12_DP
    WRITE(*,*) " Set the mesh number of X direction(500) : "
    READ(*,*) numGridElectric
    WRITE(*,*) " Set the total time step(2000) : "
    READ(*,*) numTimeSteps			
    numGridMagnetic=    numGridElectric-1
    numCPML=    11
    Isource=    235
	
	!--------------------------------------------
	!Allocate Arrays                                                                                                     
    ALLOCATE(	Ez(numGridElectric), Hy(numGridMagnetic), &
                Ez_X(numGridElectric-1), Hy_X(numGridMagnetic-1), &
                J_d(numGridElectric-2), Ez_asterisk(numGridElectric), STAT= Ier)
    IF( Ier /= 0 ) STOP 'ALLOCATION ERROR 1' 
  
    ALLOCATE(	Epsilon_r(numGridElectric-2), Sigma(numGridElectric-2), &
                Tau(numGridElectric-2), DeltaEpsilon(numGridElectric-2), STAT= Ier)
    IF( Ier /= 0 ) STOP 'ALLOCATION ERROR 2' 
  
    ALLOCATE(	C_a(numGridElectric-2), C_b(numGridElectric-2), Beta_d(numGridElectric-2), k_d(numGridElectric-2),	STAT= Ier) 
    IF( Ier /= 0 ) STOP 'ALLOCATION ERROR 3'    
  
    ALLOCATE(	Psi_Hyx_NegativeDir(numCPML), Psi_Hyx_PositiveDir(numCPML), &						
                bhx_negativeDir(numCPML), bhx_positiveDir(numCPML), chx_negativeDir(numCPML), chx_positiveDir(numCPML), &	
                Sigma_xh_negativeDir(numCPML), Sigma_xh_positiveDir(numCPML), &					
                Kappa_xh_negativeDir(numCPML), Kappa_xh_positiveDir(numCPML), &					
                a_xh_negativeDir(numCPML), a_xh_positiveDir(numCPML), &									
                Psi_Ezx_NegativeDir(numCPML), Psi_Ezx_PositiveDir(numCPML), 	&
                bex_negativeDir(numCPML), bex_positiveDir(numCPML), cex_negativeDir(numCPML), cex_positiveDir(numCPML), & 
                Sigma_xe_negativeDir(numCPML), Sigma_xe_positiveDir(numCPML), &					
                Kappa_xe_negativeDir(numCPML), Kappa_xe_positiveDir(numCPML), &					
                a_xe_negativeDir(numCPML), a_xe_positiveDir(numCPML), STAT= Ier) 				
    IF( Ier /= 0 ) STOP 'ALLOCATION ERROR 4'
  
    ALLOCATE(X(numGridElectric), STAT= Ier) 
    IF( Ier /= 0 ) STOP 'ALLOCATION ERROR 5'
  
    !--------------------------------------------
    !Set CPML Environment
    m=	3
    ma=	1
    Sigma_xmax_freeSpace=   0.8_DP*(DBLE(m+1)) / (eta_0*DeltaX*SQRT(Epsilon_freespace*Mu_freespace))
    Sigma_xmax_debye=   0.8_DP*(DBLE(m+1)) / (eta_0*DeltaX*SQRT(Epsilon_infinite_Debye*Mu_freespace))
    Kappa_xmax= 45.0_DP
    a_xmax= 0.0_DP
    Psi_Hyx_NegativeDir(:)=	0.0_DP
    Psi_Hyx_PositiveDir(:)=	0.0_DP
    Psi_Ezx_NegativeDir(:)=	0.0_DP
    Psi_Ezx_PositiveDir(:)=	0.0_DP
    
    DO I=	1,	numCPML !Set b,c update coefficents in CPML
  	    !Negative Direction
  	    Sigma_xh_negativeDir(I)=	Sigma_xmax_freeSpace*((DBLE(numCPML)-DBLE(I)+0.5_DP)/DBLE(numCPML-0.5_DP))**m  
        Sigma_xe_negativeDir(I)=	Sigma_xmax_freeSpace*((DBLE(numCPML)-DBLE(I))/DBLE(numCPML-0.5_DP))**m                                                                       
  	    Kappa_xh_negativeDir(I)=	1.0_DP + (Kappa_xmax-1.0_DP)*((DBLE(numCPML)-DBLE(I)+0.5_DP)/DBLE(numCPML-0.5_DP))**m                                                                         
  	    Kappa_xe_negativeDir(I)=	1.0_DP + (Kappa_xmax-1.0_DP)*((DBLE(numCPML)-DBLE(I))/DBLE(numCPML-0.5_DP))**m
  	    a_xh_negativeDir(I)=    a_xmax*((DBLE(I)-1.0_DP)/DBLE(numCPML-0.5_DP))**ma
  	    a_xe_negativeDir(I)=    a_xmax*((DBLE(I)-0.5_DP)/DBLE(numCPML-0.5_DP))**ma
  	
  	    !Positive Direction
  	    Sigma_xh_positiveDir(I)=	Sigma_xmax_debye*((DBLE(I)-0.5_DP)/DBLE(numCPML-0.5_DP))**m    
        Sigma_xe_positiveDir(I)=	Sigma_xmax_debye*((DBLE(I)-1.0_DP)/DBLE(numCPML-0.5_DP))**m                                                                     
  	    Kappa_xh_positiveDir(I)=	1.0_DP + (Kappa_xmax-1.0_DP)*(((DBLE(I)+0.5_DP)-1.0_DP)/DBLE(numCPML-0.5_DP))**m                                                                  
  	    Kappa_xe_positiveDir(I)=	1.0_DP + (Kappa_xmax-1.0_DP)*((DBLE(I)-1.0_DP)/DBLE(numCPML-0.5_DP))**m
  	    a_xh_positiveDir(I)=    a_xmax*((DBLE(numCPML)-(DBLE(I)))/DBLE(numCPML-0.5_DP))**ma
  	    a_xe_positiveDir(I)=    a_xmax*((DBLE(numCPML)-DBLE(I)+0.5_DP)/DBLE(numCPML-0.5_DP))**ma   
    END DO
    !Negative Direction
    bhx_negativeDir(:)=	EXP(-(Sigma_xh_negativeDir(:)/Kappa_xh_negativeDir(:)+a_xh_negativeDir(:))*DeltaT/Epsilon_0)
    bhx_positiveDir(:)=	EXP(-(Sigma_xh_positiveDir(:)/Kappa_xh_positiveDir(:)+a_xh_positiveDir(:))*DeltaT/Epsilon_0)
    chx_negativeDir(:)=	Sigma_xh_negativeDir(:)/(Sigma_xh_negativeDir(:)+Kappa_xh_negativeDir(:)*a_xh_negativeDir(:))/Kappa_xh_negativeDir(:)*(bhx_negativeDir-1.0_DP)
    chx_positiveDir(:)=	Sigma_xh_positiveDir(:)/(Sigma_xh_positiveDir(:)+Kappa_xh_positiveDir(:)*a_xh_positiveDir(:))/Kappa_xh_positiveDir(:)*(bhx_positiveDir-1.0_DP)
    !Positive Direction
    bex_negativeDir(:)=	EXP(-(Sigma_xe_negativeDir(:)/Kappa_xe_negativeDir(:)+a_xe_negativeDir(:))*DeltaT/Epsilon_0)
    bex_positiveDir(:)=	EXP(-(Sigma_xe_positiveDir(:)/Kappa_xe_positiveDir(:)+a_xe_positiveDir(:))*DeltaT/Epsilon_0)
    cex_negativeDir(:)=	Sigma_xe_negativeDir(:)/(Sigma_xe_negativeDir(:)+Kappa_xe_negativeDir(:)*a_xe_negativeDir(:))/Kappa_xe_negativeDir(:)*(bex_negativeDir-1.0_DP)
    cex_positiveDir(:)=	Sigma_xe_positiveDir(:)/(Sigma_xe_positiveDir(:)+Kappa_xe_positiveDir(:)*a_xe_positiveDir(:))/Kappa_xe_positiveDir(:)*(bex_positiveDir-1.0_DP)
    
    CALL ReplaceNaNWithZero(cex_negativeDir)
    CALL ReplaceNaNWithZero(cex_positiveDir)
    !Output parameters for validation
    CALL Outputfile1D('Sigma_xh_negativeDir', DeltaX, Sigma_xh_negativeDir)
    CALL Outputfile1D('Sigma_xh_positiveDir', DeltaX, Sigma_xh_positiveDir(SIZE(Sigma_xh_positiveDir):1:-1))
    CALL Outputfile1D('Sigma_xe_negativeDir', DeltaX, Sigma_xe_negativeDir)
    CALL Outputfile1D('Sigma_xe_positiveDir', DeltaX, Sigma_xe_positiveDir(SIZE(Sigma_xe_positiveDir):1:-1))
    CALL Outputfile1D('Kappa_xh_negativeDir', DeltaX, Kappa_xh_negativeDir)
    CALL Outputfile1D('Kappa_xh_positiveDir', DeltaX, Kappa_xh_positiveDir(SIZE(Kappa_xh_positiveDir):1:-1))
    CALL Outputfile1D('Kappa_xe_negativeDir', DeltaX, Kappa_xe_negativeDir)
    CALL Outputfile1D('Kappa_xe_positiveDir', DeltaX, Kappa_xe_positiveDir(SIZE(Kappa_xe_positiveDir):1:-1))
    CALL Outputfile1D('a_xh_negativeDir', DeltaX, a_xh_negativeDir)
    CALL Outputfile1D('a_xh_positiveDir', DeltaX, a_xh_positiveDir(SIZE(a_xh_positiveDir):1:-1))
    CALL Outputfile1D('a_xe_negativeDir', DeltaX, a_xe_negativeDir)
    CALL Outputfile1D('a_xe_positiveDir', DeltaX, a_xe_positiveDir(SIZE(a_xe_positiveDir):1:-1))
    CALL Outputfile1D('ch_x_negativeDir', DeltaX, chx_negativeDir)
    CALL Outputfile1D('ch_x_positiveDir', DeltaX, chx_positiveDir(SIZE(chx_positiveDir):1:-1))
    CALL Outputfile1D('bhx_negativeDir', DeltaX, bhx_negativeDir)
    CALL Outputfile1D('bhx_positiveDir', DeltaX, bhx_positiveDir(SIZE(bhx_positiveDir):1:-1))
    CALL Outputfile1D('ce_x_negativeDir', DeltaX, cex_negativeDir)
    CALL Outputfile1D('ce_x_positiveDir', DeltaX, cex_positiveDir(SIZE(cex_positiveDir):1:-1))
    CALL Outputfile1D('bex_negativeDir', DeltaX, bex_negativeDir)
    CALL Outputfile1D('bex_positiveDir', DeltaX, bex_positiveDir(SIZE(bex_positiveDir):1:-1))
  
    !--------------------------------------------
    !Initialize Field Values
    Ez(:)=	0.0_DP
    Hy(:)=	0.0_DP
    Ez_X(:)=0.0_DP
    Hy_X(:)=0.0_DP
    J_d(:)= 0.0_DP
    Ez_asterisk(:)=	0.0_DP
  
    !--------------------------------------------
    !Initialize Material Parameter Throughout Whole Computational Area
    Epsilon_r(1:SIZE(Epsilon_r)/2)= Epsilon_freespace ! vacuum permeability filled with half of computational region
    Epsilon_r(SIZE(Epsilon_r)/2+1:SIZE(Epsilon_r))= Epsilon_infinite_Debye ! debye permeability filled with half of computational region
    Sigma(:)= 0.0_DP
    Tau(1:SIZE(Tau)/2)=	0.0_DP
    Tau(SIZE(Tau)/2+1:SIZE(Tau))=	Tau_Debye
    DeltaEpsilon(1:SIZE(DeltaEpsilon)/2)= 0.0_DP
    DeltaEpsilon(SIZE(DeltaEpsilon)/2+1:SIZE(DeltaEpsilon))= DeltaEpsilon_Debye
  
    !--------------------------------------------
    !Update Coefficients
    k_d(:)=	(2.0_DP*Tau(:)-DeltaT)/(2.0_DP*Tau(:)+DeltaT)
    Beta_d(:)=	(2.0_DP*DeltaT*Epsilon_0*DeltaEpsilon(:))/(DeltaT+2.0_DP*Tau(:))
    C_a(:)=	(2.0_DP*Epsilon_0*Epsilon_r(:)-Sigma(:)*DeltaT+2.0_DP*Beta_d(:))/(2.0_DP*Epsilon_0*Epsilon_r(:)-Sigma(:)*DeltaT+2.0_DP*Beta_d(:))
    C_b(:)=	(2.0_DP*DeltaT)/(2.0_DP*Epsilon_0*Epsilon_r(:)+Sigma(:)*DeltaT+2.0_DP*Beta_d(:))
    D_a=	1.0_DP
    D_b=	DeltaT/(Mu_freespace*Mu_0)
    
    OPEN( UNIT = 88, FILE = ' Sample_point.PLT ', STATUS = 'REPLACE')
    
    !--------------------------------------------
    !Time Step Initialization, iteration begins
    WRITE(*,*) "Time Step Initiate, Enter it."
    Pause
    DO N=	1,	numTimeSteps
  	    !CPML terms update at time step (n+1/2) 
        CALL FIRST_DERIVATIVE_OF(Hy, DeltaX, Hy_X)
  	    Psi_Ezx_NegativeDir(:)=	bex_negativeDir(:)*Psi_Ezx_NegativeDir(:)	+	cex_negativeDir(:)*Hy_X(1:numCPML) 
  	    Psi_Ezx_PositiveDir(:)=	bex_positiveDir(:)*Psi_Ezx_PositiveDir(:)	+	cex_positiveDir(:)*Hy_X(SIZE(Hy_X)-numCPML+1:SIZE(Hy_X)) 
        
        !Electic Fields Update at time step (n+1/2) 
        Ez_asterisk(:)=	Ez(:)
        CALL FIRST_DERIVATIVE_WITH_KAPPA_OF(Hy_X, Kappa_xe_negativeDir, Kappa_xe_positiveDir)
        Ez(2:numGridElectric-1)=	C_a(:)*Ez(2:numGridElectric-1)	+	C_b(:)*(Hy_X(:)-0.5_DP*(1.0_DP+k_d(:))*J_d(:)) ! update Ez in computational domain
        Ez(2:numCPML+1)=	Ez(2:numCPML+1)	+	C_b(1:numCPML)*Psi_Ezx_NegativeDir(:) ! update Ez in CPML field at negative direction
  	    Ez(numGridElectric-numCPML:numGridElectric-1)=	Ez(numGridElectric-numCPML:numGridElectric-1)	+	C_b(SIZE(C_b)-numCPML+1:SIZE(C_b))*Psi_Ezx_PositiveDir(:) ! update Ez in CPML field at positive direction
        
        !Polarize Current Update at time step (n+1/2) 
        J_d(:)=	k_d(:)*J_d(:) + Beta_d(:)*(Ez(2:numGridElectric-1)-Ez_asterisk(2:numGridElectric-1))/DeltaT
        
        !Generate Source
  	    CALL GENERATE_SOURCE(Ez(Isource)) 
        
        !CPML term update at Time step (n+1) 
        CALL FIRST_DERIVATIVE_OF(Ez, DeltaX, Ez_X)
        Psi_Hyx_NegativeDir(:)=	bhx_negativeDir(:)*Psi_Hyx_NegativeDir(:)	+	chx_negativeDir(:)*Ez_X(1:numCPML) 
  	    Psi_Hyx_PositiveDir(:)=	bhx_positiveDir(:)*Psi_Hyx_PositiveDir(:)	+	chx_positiveDir(:)*Ez_X(SIZE(Ez_X)-numCPML+1:SIZE(Ez_X))
        
        !Magnetic Fields Update at Time step (n+1)
        CALL FIRST_DERIVATIVE_WITH_KAPPA_OF(Ez_X, Kappa_xh_negativeDir, Kappa_xh_positiveDir)
  	    Hy(:)=	D_a*Hy(:) + D_b*(Ez_X(:)) ! update Hy in computational domain
  	    Hy(1:numCPML)=	Hy(1:numCPML)	+	D_b*Psi_Hyx_NegativeDir(:) ! update Hy in CPML field at negative direction
  	    Hy(numGridMagnetic-numCPML+1:numGridMagnetic)=	Hy(numGridMagnetic-numCPML+1:numGridMagnetic)	+	D_b*Psi_Hyx_PositiveDir(:) ! update Hy in CPML field at positive direction

        ! output results
        IF( MOD(N,20)==0 ) THEN
            WRITE(STEPCOUNTER,'(I6.6)') N
            PARTNAME=	'Ez_'//STEPCOUNTER
            CALL Outputfile1D(PARTNAME, DeltaX, Ez)
        END IF
        WRITE(88,*) (N*DeltaT)/(1E-09), Ez(445)
    END DO
CONTAINS
	SUBROUTINE	GENERATE_SOURCE(Amplitude)
		REAL(KIND=	DP)	:: d, t0
		REAL(KIND=	DP), INTENT(OUT)	:: Amplitude
		d=	40.0_DP*DeltaT
		t0=	3.0_DP*d
		Amplitude=	Amplitude + EXP(-((DBLE(N-1)*DeltaT-t0)/d)**2)
        
	END	SUBROUTINE
	
    SUBROUTINE	FIRST_DERIVATIVE_OF(Phi, sizeMesh, Phi_derivative)
  	    INTEGER							::	I
  	    REAL(KIND=	DP),	INTENT(IN)	:: sizeMesh, Phi(:)
  	    REAL(KIND=	DP),	INTENT(INOUT)	:: Phi_derivative(:)
  	
  	    DO I=	1, SIZE(Phi_derivative)
  		    Phi_derivative(I)=	(Phi(I+1) - Phi(I))/sizeMesh
  	    END DO
    END SUBROUTINE	FIRST_DERIVATIVE_OF
    
    SUBROUTINE	FIRST_DERIVATIVE_WITH_KAPPA_OF(Phi_derivative, kappa_negDir, kappa_posDir)
  	    INTEGER							::	I
  	    REAL(KIND=	DP),	INTENT(IN)	:: kappa_negDir(:), kappa_posDir(:)
  	    REAL(KIND=	DP),	INTENT(INOUT)	:: Phi_derivative(:)
        REAL(KIND=  DP),    DIMENSION(:):: kappa_reciprocal(SIZE(kappa_negDir))
  	
        kappa_reciprocal(:) = 1.0_DP/kappa_negDir(:)
        Phi_derivative(1:SIZE(kappa_negDir))=	Phi_derivative(1:SIZE(kappa_negDir))*kappa_reciprocal(:)
        kappa_reciprocal(:) = 1.0_DP/kappa_posDir(:)
        Phi_derivative((SIZE(Phi_derivative)-SIZE(kappa_posDir)+1):SIZE(Phi_derivative))=	Phi_derivative((SIZE(Phi_derivative)-SIZE(kappa_posDir)+1):SIZE(Phi_derivative))*kappa_reciprocal(:)
    END SUBROUTINE	FIRST_DERIVATIVE_WITH_KAPPA_OF
                                                                                                                     
    FUNCTION L2_Norm( Q, Q_ext )                                                                                                                                                                                                        
        REAL(KIND=DP)  :: L2_Norm                                                                                        
        REAL(KIND=DP), INTENT(IN) :: Q(:), Q_ext(:)                                                                      
        !--------------------------------------------------------                                                                                                                                                                        
        L2_Norm = SQRT( SUM( (Q-Q_ext)**2 )/DBLE(SIZE(Q)) )                                                                                                                                                                                  
    END FUNCTION L2_Norm     
    
    SUBROUTINE ReplaceNaNWithZero(x)
        INTEGER::   I
        REAL(KIND=DP), INTENT(INOUT):: x(:)
        DO I = 1, SIZE(x)
            IF(x(I) /= x(I)) THEN 
                x(I) = 0
            END IF
        END DO
        
    END SUBROUTINE ReplaceNaNWithZero
                                                                                                                     
    SUBROUTINE Outputfile1D(FILE_NAME, sizeMesh, Phi)                                                                               
        IMPLICIT NONE                                                                                                                                                                                
        REAL(KIND=DP), INTENT(IN):: sizeMesh, Phi(:)                                                                                                                                         
        CHARACTER(LEN=*), INTENT(IN):: FILE_NAME
        CHARACTER:: FILENAME*30                                                                                                                                                                                          
        FILENAME=	FILE_NAME//'.PLT'                                                                          
        OPEN( UNIT= 10, FILE= FILENAME, STATUS= 'REPLACE' )                                                              
        DO I = 1, size(Phi)                                                                                                     
    	    WRITE(10,*) DBLE(I-1)*sizeMesh, Phi(I)                                                                                         
        END DO                                                                                                           
        CLOSE( UNIT= 10 )                                                                                                                                                                                                            
    END SUBROUTINE Outputfile1D                                                                                                                                                                                                                                                                                                                                
END PROGRAM ADE_FDTD_ON_YEE_CPML                                                                           