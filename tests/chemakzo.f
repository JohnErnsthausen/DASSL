      PROGRAM CHEMAKZO
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
C ...    Chemical Akzo Nobel problem
C ...    index 1 DAE of dimension 6
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
      EXTERNAL          G,DG,IC,SOLOUT
      INTEGER           NEQ,IDID,LRW,LIW
      INTEGER           IW(50),INFO(15),IPAR(10)
      DOUBLE PRECISION  ATOL,RTOL,T
      DOUBLE PRECISION  RW(500),Y(6),YP(6),RPAR(10)
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
      LRW = 500          ! Real work space.
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
      DOUBLE PRECISION K1,K2,K3,K4,KBIG,KLA,PO2,HEN,KS
      PARAMETER (
     &   K1   = 18.7D0,
     &   K2   = 0.58D0,
     &   K3   = 0.09D0,
     &   K4   = 0.42D0,
     &   KBIG = 34.4D0,
     &   KLA  = 3.3D0,
     &   KS   = 115.83D0,
     &   PO2  = 0.9D0,
     &   HEN  = 737D0
     &)
      DOUBLE PRECISION R1,R2,R3,R4,R5,FIN
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Space and dimensions.
C
      NEQ = 6            ! Dimension of problem.
      MU  = 6            ! The upper bandwidth of the Jacobian.
      ML  = 6            ! The lower bandwidth of the Jacobian.
C
C ... Initial condition.
C
      T0   = 0.0D0
      T1   = 180.0D0
      Y(1) = 0.444D0
      Y(2) = 0.00123D0
      Y(3) = 0.0D0
      Y(4) = 0.007D0
      Y(5) = 0.0D0
      Y(6) = KS*Y(1)*Y(4)
C ...
      R1  = K1*(Y(1)**4)*SQRT(Y(2))
      R2  = K2*Y(3)*Y(4)
      R3  = (K2/KBIG)*Y(1)*Y(5)
      R4  = K3*Y(1)*(Y(4)**2)
      R5  = K4*(Y(6)**2)*SQRT(Y(2))
      FIN = KLA*(PO2/HEN-Y(2))
C ...
      YP(1) = -2.0D0*R1 +R2 -R3       -R4
      YP(2) = -0.5D0*R1               -R4     -0.5D0*R5 + FIN
      YP(3) =        R1 -R2 +R3
      YP(4) =           -R2 +R3 -2.0D0*R4
      YP(5) =            R2 -R3                     +R5
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
         OPEN(UNIT= IWRT,FILE= 'chemakzo-das.dat',STATUS= 'UNKNOWN')
         OPEN(UNIT= NUMOUT,FILE= 'chemakzo-das-n.dat',STATUS= 'UNKNOWN')
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
      DOUBLE PRECISION K1,K2,K3,K4,KBIG,KLA,PO2,HEN,KS
      PARAMETER (
     &   K1   = 18.7D0,
     &   K2   = 0.58D0,
     &   K3   = 0.09D0,
     &   K4   = 0.42D0,
     &   KBIG = 34.4D0,
     &   KLA  = 3.3D0,
     &   KS   = 115.83D0,
     &   PO2  = 0.9D0,
     &   HEN  = 737D0
     &)
      DOUBLE PRECISION R1,R2,R3,R4,R5,FIN
C ...
      R1  = K1*(Y(1)**4)*SQRT(Y(2))
      R2  = K2*Y(3)*Y(4)
      R3  = (K2/KBIG)*Y(1)*Y(5)
      R4  = K3*Y(1)*(Y(4)**2)
      R5  = K4*(Y(6)**2)*SQRT(Y(2))
      FIN = KLA*(PO2/HEN-Y(2))
C ...
      GG(1) =-YP(1) -2.0D0*R1 +R2 -R3       -R4
      GG(2) =-YP(2) -0.5D0*R1               -R4     -0.5D0*R5 + FIN
      GG(3) =-YP(3)       +R1 -R2 +R3
      GG(4) =-YP(4)           -R2 +R3 -2.0D0*R4
      GG(5) =-YP(5)           +R2 -R3                     +R5
      GG(6) = KS*Y(1)*Y(4)-Y(6)
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
      DOUBLE PRECISION  RPAR(*),Y(*),YP(*),DGG(6,*)
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
      DOUBLE PRECISION K1,K2,K3,K4,KBIG,KLA,PO2,HEN,KS
      PARAMETER (
     &   K1   = 18.7D0,
     &   K2   = 0.58D0,
     &   K3   = 0.09D0,
     &   K4   = 0.42D0,
     &   KBIG = 34.4D0,
     &   KLA  = 3.3D0,
     &   KS   = 115.83D0,
     &   PO2  = 0.9D0,
     &   HEN  = 737D0
     &)
      DOUBLE PRECISION R11,R12,R23,R24,R31,R35,R41,R44,R52,R56,FIN2
C ...
      IF( Y(2).LT.0.0D0 ) THEN
         RETURN
      ENDIF
C ...
      DO 20 I=1,6 
         DO 10 J=1,6 
            DGG(I,J) = 0.0D0
   10    CONTINUE
   20 CONTINUE
C ...
      R11  = 4D0*K1*(Y(1)**3)*SQRT(Y(2))
      R12  = 0.5D0*K1*(Y(1)**4)/SQRT(Y(2))
      R23  = K2*Y(4)
      R24  = K2*Y(3)
      R31  = (K2/KBIG)*Y(5)
      R35  = (K2/KBIG)*Y(1)
      R41  = K3*Y(4)**2
      R44  = 2D0*K3*Y(1)*Y(4)
      R52  = 0.5D0*K4*(Y(6)**2)/SQRT(Y(2))
      R56  = 2D0*K4*Y(6)*SQRT(Y(2))
      FIN2 = -KLA
C ...
      DGG(1,1) = -CJ-2D0*R11-R31-R41
      DGG(1,2) = -2D0*R12
      DGG(1,3) = R23
      DGG(1,4) = R24-R44
      DGG(1,5) = -R35
C ...
      DGG(2,1) = -0.5D0*R11-R41
      DGG(2,2) = -CJ-0.5D0*R12-0.5D0*R52+FIN2
      DGG(2,4) = -R44
      DGG(2,6) = -0.5D0*R56
C ...
      DGG(3,1) = R11+R31
      DGG(3,2) = R12
      DGG(3,3) = -CJ-R23
      DGG(3,4) = -R24
      DGG(3,5) = R35
C ...
      DGG(4,1) = R31-2D0*R41
      DGG(4,3) = -R23
      DGG(4,4) = -CJ-R24-2D0*R44
      DGG(4,5) = R35
C ...
      DGG(5,1) = -R31
      DGG(5,2) = R52
      DGG(5,3) = R23
      DGG(5,4) = R24
      DGG(5,5) = -CJ-R35
      DGG(5,6) = R56
C ...
      DGG(6,1) = KS*Y(4)
      DGG(6,4) = KS*Y(1)
      DGG(6,6) = -1D0
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
      DOUBLE PRECISION EXACT(6),PREC(6),MAXVAL
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
 5       FORMAT(/'Chemical Akzo Nobel problem for DASSL')
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
C ... Computed using true double precision on a Cray C90
C ... 
C ... Solving Chemical Akzo Nobel problem using PSIDE
C ... 
C ... User input:
C ... 
C ... Give relative error tolerance: 1d-19
C ... Give absolute error tolerance: 1d-19
C ... 
C ... Integration characteristics:
C ... 
C ...    number of integration steps        5535
C ...    number of accepted steps           5534
C ...    number of f evaluations          100411
C ...    number of Jacobian evaluations        7
C ...    number of LU decompositions         368
C ...
C ... CPU-time used:                          30.77 sec
C ...
      Y(  1) =  0.1150794920661702d0
      Y(  2) =  0.1203831471567715d-2
      Y(  3) =  0.1611562887407974d0
      Y(  4) =  0.3656156421249283d-3
      Y(  5) =  0.1708010885264404d-1
      Y(  6) =  0.4873531310307455d-2
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

