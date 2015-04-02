      PROGRAM POLLU
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
C ... Author:  John Ernsthausen.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Adapted this file from its corresponding DASSL driver at
C ... http://www.dm.uniba.it/~testset/
C ... Documentation is available at the mentioned website.
C ...
C ...    Pollution problem
C ...    ODE of dimension 20
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
      EXTERNAL          G,DG,IC,SOLOUT
      INTEGER           NEQ,IDID,LRW,LIW
      INTEGER           IW(50),INFO(15),IPAR(10)
      DOUBLE PRECISION  ATOL,RTOL,T
      DOUBLE PRECISION  RW(1000),Y(20),YP(20),RPAR(10)
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Local variables.
C
      CHARACTER*5       TASK
      INTEGER           LISPACE,LRSPACE
      INTEGER           MAXORD,MU,ML
      INTEGER           LUN,IWRT,NUMOUT
      INTEGER           I,J,IER
      REAL              TIMTIM,GETTIM,CPUTIM
      DOUBLE PRECISION  TNEXT,T0,T1,H,POSNEG
      DOUBLE PRECISION  XMAXT,XMAX
      DOUBLE PRECISION ZER,ONE
      PARAMETER(ZER=0.0D0,ONE=1.0D0)
      DOUBLE PRECISION TFACT
      PARAMETER( TFACT=1.01D0 )
C
C ... Space and dimensions.
C
      LRW = 1000         ! Real work space.
      LIW = 50           ! Integer work space.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Setup the integration.
C
      CALL IC(NEQ,T1,T0,Y,YP,IW,RW,RTOL,ATOL,IPAR,RPAR,INFO,IER)
      T = T0
C
C ... Validate there is enough workspace.
C
      ML     = IW(1)
      MU     = IW(2)
      MAXORD = IW(3)
      LISPACE = 20+NEQ
      IF( INFO(6).EQ.0 ) THEN
         LRSPACE = 40+(MAXORD+4)*NEQ+NEQ**2
      ENDIF
      IF( INFO(6).EQ.1 ) THEN
         IF( INFO(5).EQ.1 ) THEN
           LRSPACE = 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
         ENDIF
         IF( INFO(5).EQ.0 ) THEN
           LRSPACE = 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
     &                 +2*(NEQ/(ML+MU+1)+1)
         ENDIF
      ENDIF
      IF ( LIW.LT.LISPACE .OR. LRW.LT.LRSPACE) THEN
        WRITE(6,10) LRSPACE,LRW,LISPACE,LIW
      ENDIF
 10   FORMAT(1X,'DDASSL requires',/
     1       1X,' Require real work space    =',I5/
     2       1X,' Alloted real work space    =',I5/
     3       1X,' Require integer work space =',I5/
     4       1X,' Alotted integer work space =',I5/)
C
C ... Integration direction.
C
      POSNEG = ONE
      IF( T1.LT.T0 ) POSNEG = -POSNEG
C
C ... Print out the solution at T0.
C
      TASK = 'START'
      CALL SOLOUT(TASK,NEQ,T,Y,YP,TNEXT,IW,RW,IPAR,RPAR,IER)
      IF( IER.EQ.1 ) THEN
         WRITE(6,20) IER
 20      FORMAT(///' SOLOUT returned IER = ',I5)
         GOTO 1000
      ENDIF
C
C ... Integration loop.
C
 30   CONTINUE
      CALL DDASSL(G,NEQ,T,Y,YP,TNEXT,INFO,RTOL,ATOL,IDID,
     &            RW,LRW,IW,LIW,RPAR,IPAR,DG)
      IF( IDID.EQ.1 .OR. IDID.EQ.2 .OR. IDID.EQ.3 ) THEN
         IF( T.GE.TNEXT ) THEN
C
C ......... Print out the solution at T and update TNEXT.
C
            IF( T.GE.T1 ) THEN
               TASK = 'FINAL'
            ELSE
               TASK = 'INTER'
            ENDIF
            CALL SOLOUT(TASK,NEQ,T,Y,YP,TNEXT,IW,RW,IPAR,RPAR,IER)
            IF( IER.EQ.1 ) THEN
               WRITE(6,20) IER
               GOTO 1000
            ENDIF
            IF( TASK.EQ.'FINAL') GOTO 40
C
C ......... Check if we are close to the terminal value T1.
C
            H = TNEXT-T
            IF( (T + TFACT*H - T1)*POSNEG .GT. ZER ) THEN
               TNEXT = T1
            ENDIF
         ENDIF
         INFO(1) = 1
         GOTO 30
      ELSEIF( IDID.LT.0 ) THEN 
         WRITE(6,*) 'DASSLD ERROR: DASSL RETURNED IDID = ', IDID 
         GOTO 1000
      ENDIF
  40  CONTINUE
C
C ... Setup the integration.
C
      CALL IC(NEQ,T1,T0,Y,YP,IW,RW,RTOL,ATOL,IPAR,RPAR,INFO,IER)
      T = T0
C
C ... Time the integration.
C
C ... Time the timer call.
C
      TIMTIM = GETTIM() 
      TIMTIM = GETTIM() - TIMTIM
C
C ... Start the timer for this run.
C
      CPUTIM = GETTIM() 
C
C ... Integrate.
C
 50   CONTINUE
      CALL DDASSL(G,NEQ,T,Y,YP,T1,INFO,RTOL,ATOL,IDID,
     &            RW,LRW,IW,LIW,RPAR,IPAR,DG)
      IF( IDID.EQ.1 .OR. IDID.EQ.2 .OR. IDID.EQ.3 ) THEN
         IF( T.GE.T1 ) GOTO 60
         INFO(1) = 1
         GOTO 50
      ELSEIF( IDID.LT.0 ) THEN 
         WRITE(6,*) 'DASSLD ERROR: DASSL RETURNED IDID = ', IDID 
         GOTO 1000
      ENDIF
  60  CONTINUE
C
C ... Stop the timer.
C
      CPUTIM = GETTIM() - CPUTIM - TIMTIM
      RPAR(1) = CPUTIM
C
C ... Print out the run stats.
C
      TASK = 'STATS'
      CALL SOLOUT(TASK,NEQ,T,Y,YP,T1,IW,RW,IPAR,RPAR,IER)
      IF( IER.EQ.1 ) THEN
         WRITE(6,20) IER
         GOTO 1000
      ENDIF
C
C ... Clean up.
C
 1000 CONTINUE
      CLOSE(UNIT= IWRT)
      CLOSE(UNIT= NUMOUT)
      STOP
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
C ... Functions and Subroutines.
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE IC(NEQ,T1,T0,Y,YP,IW,RW,RTOL,ATOL,IPAR,RPAR,INFO,IER)
      IMPLICIT NONE
C
C ... Declare ICs global variables.
C
      INTEGER          IER
      INTEGER          NEQ
      INTEGER          INFO(*)
      INTEGER          IW(*),IPAR(*)
      DOUBLE PRECISION RW(*),RPAR(*)
      DOUBLE PRECISION T0,T1,RTOL,ATOL
      DOUBLE PRECISION Y(*),YP(*)
C
C ... Declare ICs local variables.
C
      INTEGER   MAXORD,MU,ML
      INTEGER   LUN,IWRT,NUMOUT
      INTEGER   I
      DOUBLE PRECISION WK(20)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Space and dimensions.
C
      NEQ = 20           ! Dimension of problem.
      MU  = 20           ! The upper bandwidth of the Jacobian.
      ML  = 20           ! The lower bandwidth of the Jacobian.
C
C ... Initial condition.
C
      T0 = 0.0D0
      T1 = 60.0D0
      DO 10 I=1,20
         Y(I) = 0.0D0
         WK(I)  = 0.0D0
 10   CONTINUE
      Y( 2) = 0.2D0
      Y( 4) = 0.04D0
      Y( 7) = 0.1D0
      Y( 8) = 0.3D0
      Y( 9) = 0.01D0
      Y(17) = 0.007D0
      CALL G(T0,Y,WK,YP,IER,RPAR,IPAR)
C
C ... Define output devices.
C
      LUN    = 6         ! My output device.
      IWRT   = 11        ! My output file.
      NUMOUT = 12        ! My output file for numbers only.
      IPAR(1)= LUN
      IPAR(2)= IWRT
      IPAR(3)= NUMOUT
C
C ... Initialize tolerances RTOL, ATOL.
C
      IF( FIRST ) THEN
         FIRST = .FALSE.
         OPEN(UNIT= IWRT,FILE= 'pollu-das.dat',STATUS= 'UNKNOWN')
         OPEN(UNIT= NUMOUT,FILE= 'pollu-das-n.dat',STATUS= 'UNKNOWN')
         RTOL  = 1.0D-7
         ATOL  = 1.7D-9
         WRITE(*,'('//''' '''//',A'//'   '//')')
     &           'Give relative error tolerance: '
         READ *, RTOL
         WRITE(*,'('//''' '''//',A'//'   '//')')
     &           'Give absolute error tolerance: '
         READ *, ATOL
         WRITE(LUN,*) "RTOL = ",RTOL
         WRITE(LUN,*) "ATOL = ",ATOL
         WRITE(IWRT,*) "RTOL = ",RTOL
         WRITE(IWRT,*) "ATOL = ",ATOL
      ENDIF
C
C ... Set the INFO parameters.
C
C ... A new step.
C
      INFO(1) = 0
C
C ... Scalar RTOL and ATOL.
C
      INFO(2) = 0
C 
C ... Use DASSL intermediate-output mode (restart every step) 
C ... instead of stopping/restarting after DASSL IDID=-1 
C ... (500 STEPS TAKEN ON THIS CALL BEFORE REACHING T1).
C 
      INFO(3) = 1 
C 
C ... For an accurate solution at the endpoint we do not want the
C ... solution to be interpolated. Set RW(1) in integration loop.
C 
      INFO(4) = 1
      RW(1) = T1
C
C ... Jacobian:
C ... Full (dense) user-supplied: INFO(5)=1, INFO(6)=0.
C ... Full (dense) finite-difference-generated: INFO(5)=0, INFO(6)=0.
C ... Banded user-supplied: INFO(5)=1, INFO(6)=1.
C ... Banded finite-difference-generated: INFO(5)=0, INFO(6)=1.
C
      INFO(5) = 1
      INFO(6) = 0
      IW(1) = ML       ! Used whenever INFO(6) = 1
      IW(2) = MU       ! Used whenever INFO(6) = 1
C
C ... Specify a maximum stepsize.
C
      INFO(7) = 0
      RW(2) = 1.0D0    ! Used whenever INFO(7) = 1
C
C ... Set INFO(8)=1 to manually choose the initial stepsize.
C ... Set RW(3) to be your chosen stepsize.
C
      INFO(8) = 0
      RW(3) = 1.0D-05  ! Used whenever INFO(8) = 1
C
C ... Limit the BDF order. MAXORDER is 5.
C
      INFO(9) = 0
      MAXORD  = 5
      IW(3) = MAXORD   ! Used whenever INFO(9) = 1
C
C ... The solution is always positive.
C
      INFO(10)= 0
C
C ... Make the initial condition consistent.
C
      INFO(11)= 0
C
C ... Unused.
C
      INFO(13)= 0
      INFO(14)= 0
      INFO(15)= 0
C ...
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE G(T,Y,YP,GG,IER,RPAR,IPAR)
      IMPLICIT NONE
C
C ... Declare Gs global variables.
C
      INTEGER           IER
      INTEGER           IPAR(*)
      DOUBLE PRECISION  T
      DOUBLE PRECISION  RPAR(*),Y(*),YP(*),GG(*)
C
C ... G calculates GG= G(t,y,yp).
C ...
C ... Variables in the calling sequence:
C ... ----------------------------------
C ...
C ... Name   Type I/O  Description
C ... ----   ---- ---  -----------
C ...
C ... T       D   IN   Current time.
C ... Y       D   IN   Current point.
C ... YP      D   IN   Current velocity.
C ... GG      D   OUT  Current residual.
C ... IER     I   I/O  Error flag:
C ...                  IER = 0 -- on input.
C ...                  IER =-1 -- illegal input, DASSL will correct.
C ...                  IER =-2 -- return control to driver, IDID= -11.
C ... RPAR    D   IN   Array for transmitting parameters to the user 
C ...                  routines. Dimension set by the calling program.
C ... IPAR    I   IN   Array for transmitting parameters to the user
C ...                  routines. Dimension set by the calling program.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Declare Gs local variables.
C
      INTEGER I
      DOUBLE PRECISION R(25)
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,
     &                 K15,K16,K17,K18,K19,K20,K21,K22,K23,K24,K25
      PARAMETER (K1=.35D0,   K2=.266D2,
     &           K3=.123D5,  K4=.86D-3,
     &           K5=.82D-3,  K6=.15D5,
     &           K7=.13D-3,  K8=.24D5,
     &           K9=.165D5,  K10=.9D4,
     &           K11=.22D-1, K12=.12D5,
     &           K13=.188D1, K14=.163D5,
     &           K15=.48D7,  K16=.35D-3,
     &           K17=.175D-1,K18=.1D9,
     &           K19=.444D12,K20=.124D4,
     &           K21=.21D1,  K22=.578D1,
     &           K23=.474D-1,K24=.178D4,
     &           K25=.312D1)
C ...
      R( 1) = K1 *Y( 1)
      R( 2) = K2 *Y( 2)*Y(4)
      R( 3) = K3 *Y( 5)*Y(2)
      R( 4) = K4 *Y( 7)
      R( 5) = K5 *Y( 7)
      R( 6) = K6 *Y( 7)*Y(6)
      R( 7) = K7 *Y( 9)
      R( 8) = K8 *Y( 9)*Y(6)
      R( 9) = K9 *Y(11)*Y(2)
      R(10) = K10*Y(11)*Y(1)
      R(11) = K11*Y(13)
      R(12) = K12*Y(10)*Y(2)
      R(13) = K13*Y(14)
      R(14) = K14*Y( 1)*Y(6)
      R(15) = K15*Y( 3)
      R(16) = K16*Y( 4)
      R(17) = K17*Y( 4)
      R(18) = K18*Y(16)
      R(19) = K19*Y(16)
      R(20) = K20*Y(17)*Y(6)
      R(21) = K21*Y(19)
      R(22) = K22*Y(19)
      R(23) = K23*Y( 1)*Y(4)
      R(24) = K24*Y(19)*Y(1)
      R(25) = K25*Y(20)
C ...
      GG( 1) = -R(1)-R(10)-R(14)-R(23)-R(24)+
     &          R(2)+R(3)+R(9)+R(11)+R(12)+R(22)+R(25)
      GG( 2) = -R(2)-R(3)-R(9)-R(12)+R(1)+R(21)
      GG( 3) = -R(15)+R(1)+R(17)+R(19)+R(22)
      GG( 4) = -R(2)-R(16)-R(17)-R(23)+R(15)
      GG( 5) = -R(3)+R(4)+R(4)+R(6)+R(7)+R(13)+R(20)
      GG( 6) = -R(6)-R(8)-R(14)-R(20)+R(3)+R(18)+R(18)
      GG( 7) = -R(4)-R(5)-R(6)+R(13)
      GG( 8) =  R(4)+R(5)+R(6)+R(7)
      GG( 9) = -R(7)-R(8)
      GG(10) = -R(12)+R(7)+R(9)
      GG(11) = -R(9)-R(10)+R(8)+R(11)
      GG(12) =  R(9)
      GG(13) = -R(11)+R(10)
      GG(14) = -R(13)+R(12)
      GG(15) =  R(14)
      GG(16) = -R(18)-R(19)+R(16)
      GG(17) = -R(20)
      GG(18) =  R(20)
      GG(19) = -R(21)-R(22)-R(24)+R(23)+R(25)
      GG(20) = -R(25)+R(24)
C ...
      DO 10 I=1,20
         GG(I) = GG(I) - YP(I)
 10   CONTINUE
C ...
      IER = 0
C ...
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE DG(T,Y,YP,DGG,CJ,RPAR,IPAR)
      IMPLICIT NONE
C
C ... Declare DGs global variables.
C
      INTEGER           IPAR(*)
      DOUBLE PRECISION  T,CJ
      DOUBLE PRECISION  RPAR(*),Y(*),YP(*),DGG(20,*)
C
C ... DG calculates DGG= DG(t,y,yp)|y + cj*DG(t,y,yp)|yp.
C ...
C ... Variables in the calling sequence:
C ... ----------------------------------
C ...
C ... Name   Type I/O  Description
C ... ----   ---- ---  -----------
C ...
C ... T       D   IN   Current time.
C ... Y       D   IN   Current point.
C ... YP      D   IN   Current velocity.
C ... DGG     D   OUT  Current weighted sum of derivatives.
C ... CJ      D   IN   A weighting factor.
C ... RPAR    D   IN   Array for transmitting parameters to the user 
C ...                  routines. Dimension set by the calling program.
C ... IPAR    I   IN   Array for transmitting parameters to the user
C ...                  routines. Dimension set by the calling program.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Declare DGs local variables.
C
      INTEGER I,J
      DOUBLE PRECISION K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,
     &                 K15,K16,K17,K18,K19,K20,K21,K22,K23,K24,K25
      PARAMETER (K1=.35D0,   K2=.266D2,
     &           K3=.123D5,  K4=.86D-3,
     &           K5=.82D-3,  K6=.15D5,
     &           K7=.13D-3,  K8=.24D5,
     &           K9=.165D5,  K10=.9D4,
     &           K11=.22D-1, K12=.12D5,
     &           K13=.188D1, K14=.163D5,
     &           K15=.48D7,  K16=.35D-3,
     &           K17=.175D-1,K18=.1D9,
     &           K19=.444D12,K20=.124D4,
     &           K21=.21D1,  K22=.578D1,
     &           K23=.474D-1,K24=.178D4,
     &           K25=.312D1)
C ...
      DO 20 I=1,20 
         DO 10 J=1,20 
            DGG(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
C ...
      DGG(1,1)   = -K1-K10*Y(11)-K14*Y(6)-K23*Y(4)-K24*Y(19)
      DGG(1,11)  = -K10*Y(1)+K9*Y(2)
      DGG(1,6)   = -K14*Y(1)
      DGG(1,4)   = -K23*Y(1)+K2*Y(2)
      DGG(1,19)  = -K24*Y(1)+K22
      DGG(1,2)   = K2*Y(4)+K9*Y(11)+K3*Y(5)+K12*Y(10)
      DGG(1,13)  = K11
      DGG(1,20)  = K25
      DGG(1,5)   = K3*Y(2)
      DGG(1,10)  = K12*Y(2)
C ...
      DGG(2,4)   = -K2*Y(2)
      DGG(2,5)   = -K3*Y(2)
      DGG(2,11)  = -K9*Y(2)
      DGG(2,10)  = -K12*Y(2)
      DGG(2,19)  = K21
      DGG(2,1)   = K1
      DGG(2,2)   = -K2*Y(4)-K3*Y(5)-K9*Y(11)-K12*Y(10)
C ...
      DGG(3,1)   = K1
      DGG(3,4)   = K17
      DGG(3,16)  = K19
      DGG(3,19)  = K22
      DGG(3,3)   = -K15
C ...
      DGG(4,4)   = -K2*Y(2)-K16-K17-K23*Y(1)
      DGG(4,2)   = -K2*Y(4)
      DGG(4,1)   = -K23*Y(4)
      DGG(4,3)   = K15
C ...
      DGG(5,5)   = -K3*Y(2)
      DGG(5,2)   = -K3*Y(5)
      DGG(5,7)   = 2D0*K4+K6*Y(6)
      DGG(5,6)   = K6*Y(7)+K20*Y(17)
      DGG(5,9)   = K7
      DGG(5,14)  = K13
      DGG(5,17)  = K20*Y(6)
C ...
      DGG(6,6)   = -K6*Y(7)-K8*Y(9)-K14*Y(1)-K20*Y(17)
      DGG(6,7)   = -K6*Y(6)
      DGG(6,9)   = -K8*Y(6)
      DGG(6,1)   = -K14*Y(6)
      DGG(6,17)  = -K20*Y(6)
      DGG(6,2)   = K3*Y(5)
      DGG(6,5)   = K3*Y(2)
      DGG(6,16)  = 2D0*K18
C ...
      DGG(7,7)   = -K4-K5-K6*Y(6)
      DGG(7,6)   = -K6*Y(7)
      DGG(7,14)  = K13
C ...
      DGG(8,7)   = K4+K5+K6*Y(6)
      DGG(8,6)   = K6*Y(7)
      DGG(8,9)   = K7
C ...
      DGG(9,9)   = -K7-K8*Y(6)
      DGG(9,6)   = -K8*Y(9)
C ...
      DGG(10,10) = -K12*Y(2)
      DGG(10,2)  = -K12*Y(10)+K9*Y(11)
      DGG(10,9)  = K7
      DGG(10,11) = K9*Y(2)
C ...
      DGG(11,11) = -K9*Y(2)-K10*Y(1)
      DGG(11,2)  = -K9*Y(11)
      DGG(11,1)  = -K10*Y(11)
      DGG(11,9)  = K8*Y(6)
      DGG(11,6)  = K8*Y(9)
      DGG(11,13) = K11
C ...
      DGG(12,11) = K9*Y(2)
      DGG(12,2)  = K9*Y(11)
C ...
      DGG(13,13) = -K11
      DGG(13,11) = K10*Y(1)
      DGG(13,1)  = K10*Y(11)
C ...
      DGG(14,14) = -K13
      DGG(14,10) = K12*Y(2)
      DGG(14,2)  = K12*Y(10)
C ...
      DGG(15,1)  = K14*Y(6)
      DGG(15,6)  = K14*Y(1)
C ...
      DGG(16,16) = -K18-K19
      DGG(16,4)  = K16
C ...
      DGG(17,17) = -K20*Y(6)
      DGG(17,6)  = -K20*Y(17)
C ...
      DGG(18,17) = K20*Y(6)
      DGG(18,6)  = K20*Y(17)
C ...
      DGG(19,19) = -K21-K22-K24*Y(1)
      DGG(19,1)  = -K24*Y(19)+K23*Y(4)
      DGG(19,4)  = K23*Y(1)
      DGG(19,20) = K25
C ...
      DGG(20,20) = -K25
      DGG(20,1)  = K24*Y(19)
      DGG(20,19) = K24*Y(1)
C ...
      DO 30 I=1,20
         DGG(I,I) = DGG(I,I) - CJ
 30   CONTINUE
C ...
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE SOLOUT(TASK,NEQ,T,Y,YP,TOUT,IOUT,ROUT,IPAR,RPAR,IER)
      IMPLICIT NONE
C
C ... Declare SOLOUTs global variables.
C
      CHARACTER*5        TASK
      INTEGER            NEQ
      INTEGER            IER
      INTEGER            IOUT(*)
      INTEGER            IPAR(*)
      DOUBLE PRECISION   T,TOUT
      DOUBLE PRECISION   RPAR(*),ROUT(*),Y(*),YP(*)
C
C ... SOLOUT is a subroutine to write the solution at selected points.
C ...
C ... Variables in the calling sequence:
C ... ----------------------------------
C ...
C ... Name Type I/O  Description
C ... ---- ---- ---  -----------
C ...
C ... TASK   C  IN   Control variable.
C ... NEQ    I  IN   Dimension of Y and YP.
C ... T      D  IN   The current time.
C ... Y      D  IN   Array of dimension NX, the current point.
C ... YP     D  IN   Array of dimension NX, the current velocity.
C ... TOUT   D  OUT  The next out time.
C ... IOUT   I  IN   Integer work array, see solver.
C ... ROUT   D  IN   Real work array, see solver.
C ... RPAR   D  IN   Array for transmitting parameters to the user 
C ...                routines. Dimension set by the calling program.
C ... IPAR   I  IN   Array for transmitting parameters to the user
C ...                routines. Dimension set by the calling program.
C ... IER    I  OUT  Return indicator:
C ...                IER = 0 -- code will continue.
C ...                IER = 1 -- code is to stop.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Declare SOLOUTs local variables.
C
      INTEGER   I
      INTEGER   LUN,IWRT,NUMOUT
      INTEGER   NST,RESEVAL,JACEVAL,NNI,NETF,NCF,NGE,KUSED,HUSED,CTIME
      DOUBLE PRECISION EXACT(20),PREC(20),MAXVAL
      DOUBLE PRECISION DT
      DATA DT/0.1D0/
      DATA NST/11/,RESEVAL/12/,JACEVAL/13/,NETF/14/,NCF/15/
      DATA KUSED/8/,HUSED/7/,CTIME/4/
C
C ... Output devices.
C
      LUN    = IPAR(1)
      IWRT   = IPAR(2)
      NUMOUT = IPAR(3)
C
C ... The FINAL output.
C
      IF( TASK.EQ.'FINAL' ) THEN
         WRITE(LUN,10) IOUT(NST),T
         WRITE(LUN,40) (Y(I),I= 1,NEQ,1)
      ENDIF
C
C ... This run stats.
C
      IF( TASK.EQ.'STATS' ) THEN
         WRITE(LUN,70) IOUT(NST)+IOUT(NETF)+IOUT(NCF), IOUT(NST),
     &                 IOUT(RESEVAL), IOUT(JACEVAL),
     &                 -1, IOUT(NETF), IOUT(NCF)
         WRITE(IWRT,70) IOUT(NST)+IOUT(NETF)+IOUT(NCF), IOUT(NST),
     &                  IOUT(RESEVAL), IOUT(JACEVAL),
     &                  -1, IOUT(NETF), IOUT(NCF)
 70      FORMAT(/'Final Run Statistics:', //,
     &            'Number of integration steps        ', I9 , /,
     &            'Number of accepted steps           ', I9 , /,
     &            'Number of f evaluations            ', I9 , /,
     &            'Number of Jacobian evaluations     ', I9 , /,
     &            'Number of nonlinear iterations     ', I9 , /,
     &            'Number of error test failures      ', I9 , /,
     &            'Number of nonlinear conv. failures ', I9 )
C
C ...... Compute precision.
C
         CALL SOLUT(NEQ,T,EXACT)
         MAXVAL = LOG10(ABS(Y(1)-EXACT(1)))
         DO I=1,NEQ
            PREC(I) = LOG10(ABS(Y(I)-EXACT(I)))
            IF(PREC(I) .GT. MAXVAL) THEN
               MAXVAL = PREC(I)
            ENDIF
         ENDDO
         WRITE(LUN,80) (PREC(I), I = 1,NEQ)
         WRITE(LUN,90) ' SCD of norm(Y-Exact, infty) ',ABS(MAXVAL)
         WRITE(IWRT,80) (PREC(I), I = 1,NEQ)
         WRITE(IWRT,90) ' SCD of norm(Y-Exact, infty) ',ABS(MAXVAL)
 80      FORMAT(6(1X,G10.4))
 90      FORMAT(A,1X,F10.2,1X,F10.2,1X,F10.2)
         WRITE(LUN,100)' CPU-TIME sec (microsecond accuracy) ',RPAR(1)
         WRITE(IWRT,100)' CPU-TIME sec (microsecond accuracy) ',RPAR(1)
100      FORMAT(A,F8.4)
C
C ...... Only the stats are expected. Then return.
C
         RETURN
      ENDIF
C
C ... The START header.
C
      IF( TASK.EQ.'START' ) THEN
         IOUT(NST)   = 0
         IOUT(KUSED) = 0
         ROUT(HUSED) = 0.0D0
         WRITE(LUN,5)
         WRITE(IWRT,5)
 5       FORMAT(/'Pollution problem for DASSL')
         WRITE(LUN,10) IOUT(NST),T
         WRITE(LUN,40) (Y(I),I= 1,NEQ,1)
      ENDIF
C
C ... Write the continuation point.
C
      WRITE(IWRT,10) IOUT(NST)+IOUT(NETF)+IOUT(NCF),T
      WRITE(IWRT,40) (Y(I),I= 1,NEQ,1)
      WRITE(IWRT,50) ROUT(CTIME),ROUT(HUSED),IOUT(KUSED)
      WRITE(NUMOUT,60) T,(Y(I),I=1,NEQ,1)
 10   FORMAT(1X,' nstep= ',I9,' Interpolation Time= ',G16.9 )
 40   FORMAT(1X,' Solution:'/(1X,4G16.9))
 50   FORMAT(1X,' Step Time= ',G16.9,' Step= ',G16.9,' Order= ',I9/)
 60   FORMAT(9(1X,G16.9))
C
C ... Increment time.
C
      TOUT = T + DT
C
      IER = 0
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      SUBROUTINE SOLUT(NEQ,T,Y)
      IMPLICIT NONE
C
C ... Declare SOLUTs global variables.
C
      INTEGER            NEQ
      DOUBLE PRECISION   T
      DOUBLE PRECISION   Y(*)
C
C ... SOLOUT is a subroutine to write the solution at selected points.
C ...
C ... Variables in the calling sequence:
C ... ----------------------------------
C ...
C ... Name Type I/O  Description
C ... ---- ---- ---  -----------
C ...
C ... NEQ    I  IN   Dimension of Y and YP.
C ... T      D  IN   The current time.
C ... Y      D  IN   Array of dimension NX, the current point.
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Computed using true double precision RADAU5 on Cray C90
C ...      uround = work(1) = 1.01d-19
C ...      rtol = atol = h0 = 1.1d-18
C
      y( 1) = 0.5646255480022769d-01
      y( 2) = 0.1342484130422339d+00
      y( 3) = 0.4139734331099427d-08
      y( 4) = 0.5523140207484359d-02
      y( 5) = 0.2018977262302196d-06
      y( 6) = 0.1464541863493966d-06
      y( 7) = 0.7784249118997964d-01
      y( 8) = 0.3245075353396018d+00
      y( 9) = 0.7494013383880406d-02
      y(10) = 0.1622293157301561d-07
      y(11) = 0.1135863833257075d-07
      y(12) = 0.2230505975721359d-02
      y(13) = 0.2087162882798630d-03
      y(14) = 0.1396921016840158d-04
      y(15) = 0.8964884856898295d-02
      y(16) = 0.4352846369330103d-17
      y(17) = 0.6899219696263405d-02
      y(18) = 0.1007803037365946d-03
      y(19) = 0.1772146513969984d-05
      y(20) = 0.5682943292316392d-04
C ...
      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==
C
      REAL FUNCTION GETTIM()
      IMPLICIT NONE
C
C ... GETTIM IS SUPPOSED TO DELIVER THE CPU TIME (USER TIME)
C ... BEWARE OF THE CLOCK RESOLUTION: TYPICALLY 1/50 - 1/100 S
C ...
C ... IN FORTRAN 77 TIMING IS HIGHLY MACHINE/COMPILER SPECIFIC
C ... ON SOME MACHINES YOU MIGHT EVEN HAVE TO USE A C FUNCTION INSTEAD
C

C ... DUMMY; RETURNS 0
C     WRITE(*,*) 'CURRENTLY NO TIMING; ACTIVATE TIMING BY CONFIGURING',
C    &           'THE FUNCTION GETTIM'
C     GETTIM = 0
C
C ... IF THE COMPILER SUPPORTS THE INTRINSIC CPU_TIME
C
C      CALL CPU_TIME(GETTIM)
C
C ... BEST BET FOR UNIX
C
      REAL*4 ETIME
      REAL*4 TOTAL, TARRAY(2)
      TOTAL = ETIME(TARRAY)
      GETTIM = TARRAY(1)
C
C ... CRAY
C
C     REAL SECOND
C     GETTIM = SECOND()

      RETURN
      END
C
C234567==1=========2=========3=========4=========5=========6=========7==

