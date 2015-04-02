      PROGRAM MEDAKZO
      IMPLICIT NONE
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
C ...    Medical Akzo Nobel problem
C ...    ODE of dimension 400
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
      EXTERNAL          G,DG,IC,SOLOUT
      INTEGER           NEQ,IDID,LRW,LIW
      INTEGER           IW(500),INFO(15),IPAR(10)
      DOUBLE PRECISION  ATOL,RTOL,T
      DOUBLE PRECISION  RW(500000),Y(400),YP(400),RPAR(10)
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
      LRW = 500000     ! Real work space.
      LIW = 500        ! Integer work space.
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
      INTEGER   J
      DOUBLE PRECISION WK(400)
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Space and dimensions.
C
      NEQ = 400            ! Dimension of problem.
      MU  = 2              ! The upper bandwidth of the Jacobian.
      ML  = 2              ! The lower bandwidth of the Jacobian.
C
C ... Initial condition.
C
      T0 = 0.0D0
      T1 = 20.0D0
      DO 10 J=1,NEQ/2
         Y(2*J-1) = 0.0D0
         Y(2*J)   = 1.0D0
         WK(2*J-1)= 0.0D0
         WK(2*J)  = 0.0D0
  10  CONTINUE
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
         OPEN(UNIT= IWRT,FILE= 'medakzo-das.dat',STATUS= 'UNKNOWN')
         OPEN(UNIT= NUMOUT,FILE= 'medakzo-das-n.dat',STATUS= 'UNKNOWN')
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
      INFO(6) = 1
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
      INTEGER N,NEQ,I,J,L
      DOUBLE PRECISION ZETA,DZETA,DZETA2,KK,CC,PHI,ALPHA,BETA,DUM
      PARAMETER(KK=100.0D0,CC=4.0D0)
C ...
      NEQ    = 400
      N      = NEQ/2
      DZETA  = 1.0D0/DBLE(N)
      DZETA2 = DZETA*DZETA
      DUM    = (DZETA-1.0D0)*((DZETA-1.0D0)/CC)
      ALPHA  = 2.0D0*(DZETA-1.0D0)*(DUM/CC)
      BETA   = DUM*DUM
C ...
      IF( T.LE.5.0D0 ) THEN
         PHI = 2.0D0
      ELSE
         PHI = 0.0D0
      ENDIF
C ...
      GG(1) = (PHI-2.0D0*Y(1)+Y(3))*(BETA/DZETA2)
     &        +ALPHA*(Y(3)-PHI)/(2.0D0*DZETA)-KK*Y(1)*Y(2)
      GG(2) = -KK*Y(2)*Y(1)
C ...
      DO 10 J=2,N-1
         ZETA  = DBLE(J)*DZETA
         DUM   = (ZETA-1.0D0)*((ZETA-1.0D0)/CC)
         ALPHA = 2.0D0*(ZETA-1D0)*(DUM/CC)
         BETA  = DUM*DUM
         I     = 2*J-1
         GG(I) = (Y(I-2)-2.0D0*Y(I)+Y(I+2))*(BETA/DZETA2)
     &           +ALPHA*(Y(I+2)-Y(I-2))/(2.0D0*DZETA)
     &           -KK*Y(I)*Y(I+1)
       GG(I+1) = -KK*Y(I+1)*Y(I)
   10 CONTINUE
C ...
      I     = 2*N-1
      GG(I) = -KK*Y(I+1)*Y(I)
      GG(I+1) = -KK*Y(I+1)*Y(I)
C ...
      DO 20 I=1,NEQ
         GG(I) = GG(I) - YP(I)
 20   CONTINUE
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
      DOUBLE PRECISION  RPAR(*),Y(*),YP(*),DGG(7,*)
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
      INTEGER N,NEQ,ML,MU,M,K,I,J
      DOUBLE PRECISION ZETA,DZETA,DZETA2,ALPHA,BETA,KK,CC,DUM,BZ
      PARAMETER(KK=100.0D0,CC=4.0D0)
C ...
      NEQ= 400
      ML = 2
      MU = 2
      M = ML + MU + 1
C ... Initialize DGG to zero
      DO 20 J=1,NEQ
         DO 10 I=1,7         ! 2*ML+MU+1
            DGG(I,J) = 0D0
   10    CONTINUE
   20 CONTINUE
C ... I believe there is no overwriting!
      N      = NEQ/2
      DZETA  = 1D0/DBLE(N)
      DZETA2 = DZETA*DZETA
      DUM    = (DZETA-1D0)*(DZETA-1D0)/CC
      ALPHA  = 2D0*(DZETA-1D0)*DUM/CC
      BETA   = DUM*DUM
C ... I=1 and J=1,2,3 => K=5,4,3
      DGG(5,1) = -BETA*2.0D0/DZETA2-KK*Y(2)
      DGG(4,2) = -KK*Y(1)
      DGG(3,3) = BETA/DZETA2+ALPHA/(2.0D0*DZETA)
C ... I=2 and J=1,2 => K=6,5
      DGG(6,1) = -KK*Y(2)
      DGG(5,2) = -KK*Y(1)
C ...
      DO 30 J=2,N-1
         I          = 2*J-1
         ZETA       = J*DZETA
         DUM        = (ZETA-1.0D0)*((ZETA-1.0D0)/CC)
         ALPHA      = 2.0D0*(ZETA-1.0D0)*(DUM/CC)
         BETA       = DUM*DUM
         BZ         = BETA/DZETA2
C ...... I=2*J-1 and J=I-2,I,I+1,I+2 => K=7,5,4,3
         DGG(7,I-2) = BZ-ALPHA/(2.0D0*DZETA)
         DGG(5,I)   = -2.0D0*BZ-KK*Y(I+1)
         DGG(4,I+1) = -KK*Y(I)
         DGG(3,I+2) = BZ+ALPHA/(2.0D0*DZETA)
C ...... I=2*J and J=2*J-1,2*J => K=6,5
         DGG(6,I)   = -KK*Y(I+1)
         DGG(5,I+1) = -KK*Y(I)
   30 CONTINUE
C ... I=2*N-1 and J=2*N-1,2*N => K=5,4
      DGG(5,2*N-1) = -KK*Y(2*N)
      DGG(4,2*N)   = -KK*Y(2*N-1)
C ... I=2*N and J=2*N-1,2*N => K=6,5
      DGG(6,2*N-1) = -KK*Y(2*N)
      DGG(5,2*N)   = -KK*Y(2*N-1)
C ... I=1:NEQ and J=I => K=5
      DO 40 I=1,NEQ
         DGG(5,I) = DGG(5,I) - CJ
 40   CONTINUE
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
      DOUBLE PRECISION EXACT(400),PREC(400),MAXVAL
      DOUBLE PRECISION DT
      DATA DT/0.01D0/
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
 5       FORMAT(/'Medical Akzo Nobel problem for DASSL')
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
C ... Computed at Cray C90 using Cray double precision
C ... Solving Medical Akzo Nobel problem using PSIDE
C ...
C ... User input:
C ...
C ... Give relative error tolerance: 1d-10
C ... Give absolute error tolerance: 1d-10
C ...
C ...
C ... Integration characteristics:
C ...
C ...    number of integration steps         551
C ...    number of accepted steps            537
C ...    number of f evaluations            8914
C ...    number of Jacobian evaluations       21
C ...    number of LU decompositions         620
C ...
C ... CPU-time used:                          78.24 sec
C
      y(  1) =  0.5113983840919909d-005
      y(  2) =  0.1925112884312553d-143
      y(  3) =  0.1027858770570419d-004
      y(  4) =  0.1890518289312031d-142
      y(  5) =  0.1549349862635799d-004
      y(  6) =  0.1774199325357386d-142
      y(  7) =  0.2075835344757462d-004
      y(  8) =  0.5897341137981092d-143
      y(  9) =  0.2607273610116854d-004
      y( 10) =  0.1093527900908030d-143
      y( 11) =  0.3143617475695002d-004
      y( 12) =  0.1188834841626416d-144
      y( 13) =  0.3684813884509626d-004
      y( 14) =  0.9968323236025642d-147
      y( 15) =  0.4230803594492533d-004
      y( 16) = -0.2801994001528093d-146
      y( 17) =  0.4781520853483223d-004
      y( 18) = -0.7337417669341249d-147
      y( 19) =  0.5336893059800053d-004
      y( 20) = -0.1209033101530330d-147
      y( 21) =  0.5896840407836044d-004
      y( 22) = -0.1430357497530360d-148
      y( 23) =  0.6461275518112516d-004
      y( 24) = -0.1063952641824646d-149
      y( 25) =  0.7030103051210320d-004
      y( 26) =  0.7939969136126717d-152
      y( 27) =  0.7603219304985662d-004
      y( 28) =  0.1568246940545520d-150
      y( 29) =  0.8180511794465543d-004
      y( 30) =  0.4074950357924872d-150
      y( 31) =  0.8761858813806752d-004
      y( 32) =  0.5592746648679992d-150
      y( 33) =  0.9347128979692480d-004
      y( 34) = -0.5510388943414421d-151
      y( 35) =  0.9936180755532036d-004
      y( 36) = -0.2724738349250769d-149
      y( 37) =  0.1052886195582220d-003
      y( 38) = -0.9327772452398718d-149
      y( 39) =  0.1112500923002360d-003
      y( 40) = -0.2182885200987554d-148
      y( 41) =  0.1172444752530255d-003
      y( 42) = -0.4041450806475518d-148
      y( 43) =  0.1232698952748828d-003
      y( 44) = -0.5608157478395261d-148
      y( 45) =  0.1293243507959787d-003
      y( 46) = -0.2639662630908699d-148
      y( 47) =  0.1354057057728661d-003
      y( 48) =  0.1801866277537073d-147
      y( 49) =  0.1415116834059119d-003
      y( 50) =  0.8464449882759417d-147
      y( 51) =  0.1476398596134615d-003
      y( 52) =  0.2245234937355967d-146
      y( 53) =  0.1537876562567258d-003
      y( 54) =  0.3359213489153582d-146
      y( 55) =  0.1599523341096154d-003
      y( 56) = -0.3085721171916412d-146
      y( 57) =  0.1661309855680449d-003
      y( 58) = -0.4465322607423735d-145
      y( 59) =  0.1723205270935920d-003
      y( 60) = -0.1970925996866384d-144
      y( 61) =  0.1785176913868402d-003
      y( 62) = -0.6070953121563027d-144
      y( 63) =  0.1847190192862588d-003
      y( 64) = -0.1412011918930335d-143
      y( 65) =  0.1909208513890961d-003
      y( 66) = -0.2378861987352203d-143
      y( 67) =  0.1971193193914910d-003
      y( 68) = -0.2380432473186974d-143
      y( 69) =  0.2033103371458565d-003
      y( 70) = -0.6522557638254663d-145
      y( 71) =  0.2094895914345677d-003
      y( 72) =  0.1784305601809064d-143
      y( 73) =  0.2156525324601176d-003
      y( 74) = -0.1007474781780816d-142
      y( 75) =  0.2217943640531935d-003
      y( 76) = -0.5281511349479423d-142
      y( 77) =  0.2279100336016016d-003
      y( 78) = -0.1117525482975987d-141
      y( 79) =  0.2339942217046434d-003
      y( 80) = -0.1127916494884468d-141
      y( 81) =  0.2400413315594459d-003
      y( 82) = -0.1633306916231411d-142
      y( 83) =  0.2460454780878912d-003
      y( 84) =  0.2708874035585891d-143
      y( 85) =  0.2520004768152150d-003
      y( 86) = -0.2501941069702609d-142
      y( 87) =  0.2578998325140575d-003
      y( 88) = -0.2642308070750020d-141
      y( 89) =  0.2637367276308081d-003
      y( 90) = -0.3684887530751217d-139
      y( 91) =  0.2695040105145025d-003
      y( 92) = -0.3647274179805887d-138
      y( 93) =  0.2751941834723564d-003
      y( 94) = -0.1255641406397419d-137
      y( 95) =  0.2807993906802854d-003
      y( 96) = -0.1694257216823904d-138
      y( 97) =  0.2863114059815211d-003
      y( 98) = -0.1785516142939602d-136
      y( 99) =  0.2917216206117258d-003
      y(100) = -0.3935939757647002d-135
      y(101) =  0.2970210308948898d-003
      y(102) = -0.2514765666933440d-134
      y(103) =  0.3022002259608294d-003
      y(104) = -0.7200873856605984d-134
      y(105) =  0.3072493755423352d-003
      y(106) = -0.7539683247227422d-134
      y(107) =  0.3121582179180383d-003
      y(108) =  0.3738577086039426d-135
      y(109) =  0.3169160480759169d-003
      y(110) = -0.2493582962172335d-131
      y(111) =  0.3215117061821543d-003
      y(112) =  0.3039632438293726d-130
      y(113) =  0.3259335664508512d-003
      y(114) =  0.5321044068586611d-128
      y(115) =  0.3301695265219917d-003
      y(116) = -0.1918129324351378d-126
      y(117) =  0.3342069974681551d-003
      y(118) = -0.1336929159252586d-124
      y(119) =  0.3380328945648600d-003
      y(120) =  0.9521748754010357d-123
      y(121) =  0.3416336289752354d-003
      y(122) =  0.1001197393324181d-120
      y(123) =  0.3449951005170561d-003
      y(124) =  0.2703860993866771d-119
      y(125) =  0.3481026916991771d-003
      y(126) =  0.4365133580297076d-119
      y(127) =  0.3509412632351946d-003
      y(128) =  0.4898111237855383d-115
      y(129) =  0.3534951512648823d-003
      y(130) =  0.1621439381962246d-112
      y(131) =  0.3557481665387581d-003
      y(132) =  0.3003220203772183d-110
      y(133) =  0.3576835958481664d-003
      y(134) =  0.5931668289615909d-108
      y(135) =  0.3592842060126915d-003
      y(136) =  0.2235590472383775d-105
      y(137) =  0.3605322507686931d-003
      y(138) =  0.1025457293602057d-102
      y(139) =  0.3614094809374544d-003
      y(140) =  0.3496613568296336d-100
      y(141) =  0.3618971582890092d-003
      y(142) =  0.4767073568395508d-098
      y(143) =  0.3619760735583436d-003
      y(144) = -0.2410784286794997d-095
      y(145) =  0.3616265691144918d-003
      y(146) = -0.9188398110576038d-093
      y(147) =  0.3608285668302233d-003
      y(148) =  0.1146623087995081d-089
      y(149) =  0.3595616017506735d-003
      y(150) =  0.1649638439865233d-086
      y(151) =  0.3578048622135169d-003
      y(152) =  0.1215140240350217d-083
      y(153) =  0.3555372371311931d-003
      y(154) =  0.7134490346394154d-081
      y(155) =  0.3527373712073181d-003
      y(156) =  0.4502515392738464d-078
      y(157) =  0.3493837289247301d-003
      y(158) =  0.7138395988310312d-075
      y(159) =  0.3454546682115489d-003
      y(160) =  0.9941693919247076d-071
      y(161) =  0.3409285247640208d-003
      y(162) =  0.2012859826753015d-066
      y(163) =  0.3357837080804970d-003
      y(164) =  0.3598261520662423d-062
      y(165) =  0.3299988103392750d-003
      y(166) =  0.5466580008990664d-058
      y(167) =  0.3235527293336597d-003
      y(168) =  0.6945384844951550d-054
      y(169) =  0.3164248067597393d-003
      y(170) =  0.7275415527806026d-050
      y(171) =  0.3085949832350532d-003
      y(172) =  0.6193143746524996d-046
      y(173) =  0.3000439715082906d-003
      y(174) =  0.4219255556214135d-042
      y(175) =  0.2907534493998412d-003
      y(176) =  0.2263678154715720d-038
      y(177) =  0.2807062740884081d-003
      y(178) =  0.9401607967545219d-035
      y(179) =  0.2698867194275612d-003
      y(180) =  0.2968231730793053d-031
      y(181) =  0.2582807380350103d-003
      y(182) =  0.6987463944434805d-028
      y(183) =  0.2458762499428408d-003
      y(184) =  0.1201641789884051d-024
      y(185) =  0.2326634596245027d-003
      y(186) =  0.1477169946829840d-021
      y(187) =  0.2186352032185982d-003
      y(188) =  0.1268462422099779d-018
      y(189) =  0.2037873277440060d-003
      y(190) =  0.7425015664001834d-016
      y(191) =  0.1881191040379240d-003
      y(192) =  0.2886826929895103d-013
      y(193) =  0.1716336750388461d-003
      y(194) =  0.7252477041900172d-011
      y(195) =  0.1543385408702044d-003
      y(196) =  0.1143390654212691d-008
      y(197) =  0.1362460820444338d-003
      y(198) =  0.1096625145716966d-006
      y(199) =  0.1173741304462833d-003
      y(200) =  0.6190822732534586d-005
      y(201) =  0.9774701310627047d-004
      y(202) =  0.1986273404756002d-003
      y(203) =  0.7740788649977313d-004
      y(204) =  0.3489773624098464d-002
      y(205) =  0.5657119003189305d-004
      y(206) =  0.3234526094359604d-001
      y(207) =  0.3643334879766658d-004
      y(208) =  0.1548747348410801d+000
      y(209) =  0.2003152841880950d-004
      y(210) =  0.4026980529594953d+000
      y(211) =  0.9608297851720770d-005
      y(212) =  0.6649744834198490d+000
      y(213) =  0.4215537698495267d-005
      y(214) =  0.8409284546320647d+000
      y(215) =  0.1753504402754791d-005
      y(216) =  0.9314946676956936d+000
      y(217) =  0.7048158429518009d-006
      y(218) =  0.9720896201631835d+000
      y(219) =  0.2760943506466737d-006
      y(220) =  0.9890204872799944d+000
      y(221) =  0.1057554501281432d-006
      y(222) =  0.9957930123519514d+000
      y(223) =  0.3965142250779033d-007
      y(224) =  0.9984246531478463d+000
      y(225) =  0.1455273204279008d-007
      y(226) =  0.9994229325942358d+000
      y(227) =  0.5226348147846279d-008
      y(228) =  0.9997932125999319d+000
      y(229) =  0.1835610545325733d-008
      y(230) =  0.9999275409325039d+000
      y(231) =  0.6301078589385454d-009
      y(232) =  0.9999751869380269d+000
      y(233) =  0.2112538351365564d-009
      y(234) =  0.9999917015131560d+000
      y(235) =  0.6912550453447044d-010
      y(236) =  0.9999972914302640d+000
      y(237) =  0.2205932132514696d-010
      y(238) =  0.9999991378543379d+000
      y(239) =  0.6860095639285670d-011
      y(240) =  0.9999997325855174d+000
      y(241) =  0.2077324462852526d-011
      y(242) =  0.9999999192384585d+000
      y(243) =  0.6120038908594393d-012
      y(244) =  0.9999999762710279d+000
      y(245) =  0.1752695518797070d-012
      y(246) =  0.9999999932230490d+000
      y(247) =  0.4875001992978682d-013
      y(248) =  0.9999999981203191d+000
      y(249) =  0.1315706848908981d-013
      y(250) =  0.9999999994941428d+000
      y(251) =  0.3442274192104633d-014
      y(252) =  0.9999999998680372d+000
      y(253) =  0.8721783456154470d-015
      y(254) =  0.9999999999666630d+000
      y(255) =  0.2137938962858872d-015
      y(256) =  0.9999999999918528d+000
      y(257) =  0.5064735930780995d-016
      y(258) =  0.9999999999980759d+000
      y(259) =  0.1158284928109727d-016
      y(260) =  0.9999999999995613d+000
      y(261) =  0.2554350586347124d-017
      y(262) =  0.9999999999999036d+000
      y(263) =  0.5425563935887811d-018
      y(264) =  0.9999999999999796d+000
      y(265) =  0.1108623976460997d-018
      y(266) =  0.9999999999999958d+000
      y(267) =  0.2176490922739810d-019
      y(268) =  0.9999999999999992d+000
      y(269) =  0.4100180074816888d-020
      y(270) =  0.9999999999999998d+000
      y(271) =  0.7401919443964595d-021
      y(272) =  0.1000000000000000d+001
      y(273) =  0.1278745657114596d-021
      y(274) =  0.1000000000000000d+001
      y(275) =  0.2111087049605767d-022
      y(276) =  0.1000000000000000d+001
      y(277) =  0.3325632734364699d-023
      y(278) =  0.1000000000000000d+001
      y(279) =  0.4991515592566292d-024
      y(280) =  0.1000000000000000d+001
      y(281) =  0.7126950428617158d-025
      y(282) =  0.1000000000000000d+001
      y(283) =  0.9664740804131475d-026
      y(284) =  0.1000000000000000d+001
      y(285) =  0.1242716896959521d-026
      y(286) =  0.1000000000000000d+001
      y(287) =  0.1512543532243458d-027
      y(288) =  0.1000000000000000d+001
      y(289) =  0.1739533019752215d-028
      y(290) =  0.1000000000000000d+001
      y(291) =  0.1886942537979667d-029
      y(292) =  0.1000000000000000d+001
      y(293) =  0.1926965705022792d-030
      y(294) =  0.1000000000000000d+001
      y(295) =  0.1849021812823421d-031
      y(296) =  0.1000000000000000d+001
      y(297) =  0.1663798767415642d-032
      y(298) =  0.1000000000000000d+001
      y(299) =  0.1401076830818626d-033
      y(300) =  0.1000000000000000d+001
      y(301) =  0.1101818149402153d-034
      y(302) =  0.1000000000000000d+001
      y(303) =  0.8074224739509168d-036
      y(304) =  0.1000000000000000d+001
      y(305) =  0.5501249196662931d-037
      y(306) =  0.1000000000000000d+001
      y(307) =  0.3476859813132770d-038
      y(308) =  0.1000000000000000d+001
      y(309) =  0.2033489290876775d-039
      y(310) =  0.1000000000000000d+001
      y(311) =  0.1097880013869247d-040
      y(312) =  0.1000000000000000d+001
      y(313) =  0.5457825200381417d-042
      y(314) =  0.1000000000000000d+001
      y(315) =  0.2491675366427318d-043
      y(316) =  0.1000000000000000d+001
      y(317) =  0.1041801880291617d-044
      y(318) =  0.1000000000000000d+001
      y(319) =  0.3978066491064419d-046
      y(320) =  0.1000000000000000d+001
      y(321) =  0.1383174699098532d-047
      y(322) =  0.1000000000000000d+001
      y(323) =  0.4365911791079500d-049
      y(324) =  0.1000000000000000d+001
      y(325) =  0.1247057764661705d-050
      y(326) =  0.1000000000000000d+001
      y(327) =  0.3212728839963712d-052
      y(328) =  0.1000000000000000d+001
      y(329) =  0.7439366703571565d-054
      y(330) =  0.1000000000000000d+001
      y(331) =  0.1542770387822259d-055
      y(332) =  0.1000000000000000d+001
      y(333) =  0.2854454245592573d-057
      y(334) =  0.1000000000000000d+001
      y(335) =  0.4693220411250150d-059
      y(336) =  0.1000000000000000d+001
      y(337) =  0.6828458274546624d-061
      y(338) =  0.1000000000000000d+001
      y(339) =  0.8752952529541412d-063
      y(340) =  0.1000000000000000d+001
      y(341) =  0.9838541433761416d-065
      y(342) =  0.1000000000000000d+001
      y(343) =  0.9649177728609193d-067
      y(344) =  0.1000000000000000d+001
      y(345) =  0.8213596936190817d-069
      y(346) =  0.1000000000000000d+001
      y(347) =  0.6033986647865674d-071
      y(348) =  0.1000000000000000d+001
      y(349) =  0.3802531117966294d-073
      y(350) =  0.1000000000000000d+001
      y(351) =  0.2042261117698575d-075
      y(352) =  0.1000000000000000d+001
      y(353) =  0.9282595096128614d-078
      y(354) =  0.1000000000000000d+001
      y(355) =  0.3543587864454877d-080
      y(356) =  0.1000000000000000d+001
      y(357) =  0.1126779423370979d-082
      y(358) =  0.1000000000000000d+001
      y(359) =  0.2957534367766753d-085
      y(360) =  0.1000000000000000d+001
      y(361) =  0.6344600529877694d-088
      y(362) =  0.1000000000000000d+001
      y(363) =  0.1100279075462365d-090
      y(364) =  0.1000000000000000d+001
      y(365) =  0.1523845293461783d-093
      y(366) =  0.1000000000000000d+001
      y(367) =  0.1662696161555950d-096
      y(368) =  0.1000000000000000d+001
      y(369) =  0.1407578290673998d-099
      y(370) =  0.1000000000000000d+001
      y(371) =  0.9086150803567186d-103
      y(372) =  0.1000000000000000d+001
      y(373) =  0.4384339596163745d-106
      y(374) =  0.1000000000000000d+001
      y(375) =  0.1545482064392824d-109
      y(376) =  0.1000000000000000d+001
      y(377) =  0.3874172613928345d-113
      y(378) =  0.1000000000000000d+001
      y(379) =  0.6689452219441953d-117
      y(380) =  0.1000000000000000d+001
      y(381) =  0.7655680935317283d-121
      y(382) =  0.1000000000000000d+001
      y(383) =  0.5538543899545850d-125
      y(384) =  0.1000000000000000d+001
      y(385) =  0.2386173886563501d-129
      y(386) =  0.1000000000000000d+001
      y(387) =  0.5664887497790931d-134
      y(388) =  0.1000000000000000d+001
      y(389) =  0.6671124967149171d-139
      y(390) =  0.1000000000000000d+001
      y(391) =  0.3351973480286951d-144
      y(392) =  0.1000000000000000d+001
      y(393) =  0.5684315818559200d-150
      y(394) =  0.1000000000000000d+001
      y(395) =  0.2142121793294590d-156
      y(396) =  0.1000000000000000d+001
      y(397) =  0.6727117900187205d-164
      y(398) =  0.1000000000000000d+001
      y(399) =  0.0000000000000000d+000
      y(400) =  0.1000000000000000d+001
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

