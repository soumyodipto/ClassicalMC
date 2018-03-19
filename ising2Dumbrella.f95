PROGRAM ising2D
IMPLICIT NONE
INTEGER ::i,j,m,n
INTEGER,ALLOCATABLE ::A(:,:) !Matrix containing spins
INTEGER ::nrows,ncols
REAL ::temp,beta
INTEGER ::Configtype
INTEGER ::npass
INTEGER ::ipass


INTEGER ::trial_spin

INTEGER ::win

REAL ::deltaU
REAL ::log_eta

Integer ::magnetization

INTEGER :: visit(401)

OPEN(UNIT=11,FILE='ising_umbrella.txt',STATUS='OLD')
READ(11,*)nrows
READ(11,*)ncols
READ(11,*)npass
READ(11,*)temp
READ(11,*)Configtype
CLOSE(UNIT=11)
ALLOCATE(A(nrows+2,ncols+2)) !Dimensions of spin array i.e. system size

OPEN(UNIT=20,FILE='umbrella.txt',STATUS='REPLACE')
WRITE(20,*)"magnetization   visits"



           
           ! initialize variables
           beta=1.0/temp
           
           ! set up intial spin configuration
           select case (Configtype)
             case(1) ! checkerboard configuration
             A(1,1)=1
             DO i=1, nrows+1
               A(i+1,1)=-A(i,1)
               ENDDO
             DO j=1 , ncols+1
               A(:,j+1)=-A(:,j)
               ENDDO
               case(2) !interface configuration
               DO i=1, nrows+2
                 DO j=1 , (ncols+2)/2
                 A(i,j)=1
                 ENDDO
                 DO j= (ncols+2)/2 +1 ,ncols+2
                   A(i,j)=-1
                   ENDDO
                   ENDDO
                case(3) !unequal interface
                DO i=1 , nrows+2
                  DO j=1 , (ncols+2)/4
                    A(i,j)=1
                    ENDDO
                    DO j=(ncols+2)/4 + 1 , ncols+2
                      A(i,j)=-1
                      ENDDO
                      ENDDO
                      case default
                      PRINT *,"Error"
                      STOP
                      END select
  ! Main loop containing MC program
  DO i=1,401
    visit(i)=0
    END DO
  window_loop: DO win=1,2
  PRINT *,"Running for window",win
  MC_passes: DO ipass=0 , npass
                  if(mod(ipass,100) == 0) print *, "Interval"
                    magnetization = sum(A(2:nrows+1,2:nrows+1))
                  
                        
					!Randomly choose a spin
                        IF(((win-1)*40.le.magnetization).and.(magnetization.le.win*40))THEN 
                          visit((magnetization + 1))=visit((magnetization + 1))+1
                          END IF
                        DO WHILE(((win-1)*40.le.magnetization).and.(magnetization.le.win*40))
                        m=nint((nrows-1)*ran1(5) + 2) !choose a random row
                        n=nint((ncols-1)*ran1(5) + 2) !choose a random coloumn
                        trial_spin = -A(m,n)          ! trial spin value i.e.flip
                        deltaU = -trial_spin*(A(m-1,n)+ A(m+1,n)+ A(m,n-1)+ A(m,n+1))*2
                        log_eta = DLOG(ran1(5)+1.0d-10) ! random number 0-1 with tiny offset
                        IF(-beta*deltaU > log_eta)then
                          A(m,n)=trial_spin
                          IF(m==2)A(nrows+2,n)=trial_spin
                          IF(m==nrows+1) A(1,n)=trial_spin
                            IF(n==2) A(m,ncols+2)=trial_spin
                              IF(n==ncols+1) A(m,1)=trial_spin
                                ENDIF
                                END DO 
								ENDDO MC_passes
                                
								ENDDO window_loop
 
       do i = 0, 400
         write(20,*) i, visit(i+1)
       enddo
         
         CLOSE(UNIT=20)
         
		CONTAINS
         DOUBLE PRECISION FUNCTION ran1(idum)
         IMPLICIT NONE
         DOUBLE PRECISION :: r(97)
         INTEGER ,INTENT(IN) :: idum
         save
         INTEGER , PARAMETER :: M1=259200, IA1=7141, IC1=54773
         REAL , PARAMETER ::RM1=1.0d0/M1
         INTEGER ,PARAMETER ::M2=134456, IA2=8121, IC2=28411
         REAL , PARAMETER :: RM2=1.0d0/M2
         INTEGER , PARAMETER ::M3=243000,IA3=4561,IC3=51349
         INTEGER::IX1,IX2,IX3,jjj
         INTEGER::iff=0
         IF(idum<0.or.iff==0)then
           iff=1
           IX1=mod(IC1-idum,M1)
           IX1=mod(IA1*IX1+IC1,M1)
           IX2=mod(IX1,M2)
           IX1=mod(IA1*IX1+IC1,M1)
           IX3=mod(IX1,M3)
           DO jjj=1,97
             IX1=mod(IA1*IX1+IC1,M1)
             IX2=mod(IA2*IX2+IC2,M2)
             r(jjj)= (dfloat(IX1)+dfloat(IX2)*RM2)*RM1
             ENDDO
             ENDIF
             IX1=mod(IA1*IX1+IC1,M1)
             IX2=mod(IA2*IX2+IC2,M2)
             IX3=mod(IA3*IX3+IC3,M3)
             jjj=1+(97*IX3)/M3
             IF(jjj>97 .or. jjj<1) PAUSE
               ran1=r(jjj)
               r(jjj)=(dfloat(IX1)+dfloat(IX2)*RM2)*RM1

               END FUNCTION ran1
               END PROGRAM ising2D
                            
                         
                              
                   
                        
                    
                                      
                
            

