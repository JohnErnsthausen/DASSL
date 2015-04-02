      PROGRAM TBA
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
C ...    Two bit adding unit
C ...    index 1 DAE of dimension 350
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
      EXTERNAL          G,DG,IC,SOLOUT
      INTEGER           NEQ,IDID,LRW,LIW
      INTEGER           IW(500),INFO(15),IPAR(10)
      DOUBLE PRECISION  ATOL,RTOL,T
      DOUBLE PRECISION  RW(200000),Y(350),YP(350),RPAR(10)
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
      LRW = 200000       ! Real work space.
      LIW = 500          ! Integer work space.
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
      INTEGER   I,J
      DOUBLE PRECISION U(175)
      DOUBLE PRECISION RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
       COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
      EXTERNAL GCN
      LOGICAL FIRST
      DATA FIRST/.TRUE./
      SAVE FIRST
C
C234567--1---------2---------3---------4---------5---------6---------7--
C
C ... Space and dimensions.
C
      NEQ = 350          ! Dimension of problem.
      MU  = 350          ! The upper bandwidth of the Jacobian.
      ML  = 350          ! The lower bandwidth of the Jacobian.
C
C ... Initial condition.
C
C ... y{0}=(g(U), U), where U is the vector of voltage variables
C ... g(u) is computed by GCN
      t0 = 0.0D0
      t1 = 320.0D0
C ... Scaling of time axis.
      CTIME=1.D4
C ... Stiffness factor for the voltage courses.
      STIFF=5.D0
C ... The MOS-Parameters.
      RGS=.4D2/(CTIME*STIFF)
      RGD=.4D2/(CTIME*STIFF)
      RBS=.1D3/(CTIME*STIFF)
      RBD=.1D3/(CTIME*STIFF)
      CGS=.6D-4*CTIME
      CGD=.6D-4*CTIME
      CBD=2.4D-5*CTIME
      CBS=2.4D-5*CTIME
      DELTA=0.2D-1
      CURIS=1.D-15*CTIME*STIFF
      VTH=25.85D0
      VDD=5.D0
      VBB=-2.5D0
C ... Load capacitance for the logical subciruits.
C ... CLOAD=.5D-4*CTIME
C ... No artificial load capacitances.
      CLOAD=0.0D0
C ... Load capacitance for the output nodes.
      COUT=2.D-4*CTIME-CLOAD
      U(  1) =    4.999999999996544D0
      U(  2) =    4.999999999999970D0
      U(  3) =   -2.499999999999975D0
      U(  4) =   -2.499999999999975D0
      U(  5) =    4.999999999996514D0
      U(  6) =    0.000000000000000D0
      U(  7) =    4.999999999996514D0
      U(  8) =   -2.499999999999991D0
      U(  9) =   -2.499999999999975D0
      U( 10) =    0.000000000000000D0
      U( 11) =    4.999999999996514D0
      U( 12) =   -2.499999999999991D0
      U( 13) =   -2.499999999999975D0
      U( 14) =    0.215858486765796D0
      U( 15) =    4.988182208251953D0
      U( 16) =   -2.499999999999990D0
      U( 17) =   -2.499999999999975D0
      U( 18) =    0.204040695017748D0
      U( 19) =    0.011817791748026D0
      U( 20) =    0.192222903269723D0
      U( 21) =   -2.499999999999991D0
      U( 22) =   -2.499999999999990D0
      U( 23) =   -0.228160951881239D0
      U( 24) =    0.204040695017748D0
      U( 25) =   -2.499999999999992D0
      U( 26) =   -2.499999999999990D0
      U( 27) =   -0.228160951881241D0
      U( 28) =    0.000000000000000D0
      U( 29) =   -0.228160951881239D0
      U( 30) =   -2.499999999999991D0
      U( 31) =   -2.499999999999992D0
      U( 32) =    4.999999999996547D0
      U( 33) =    4.999999999999970D0
      U( 34) =   -2.499999999999975D0
      U( 35) =   -2.499999999999975D0
      U( 36) =    4.999999999996517D0
      U( 37) =    0.000000000000000D0
      U( 38) =    4.999999999996517D0
      U( 39) =   -2.499999999999991D0
      U( 40) =   -2.499999999999975D0
      U( 41) =    0.000000000000000D0
      U( 42) =    4.999999999996517D0
      U( 43) =   -2.499999999999991D0
      U( 44) =   -2.499999999999975D0
      U( 45) =    0.215858484247529D0
      U( 46) =    4.988182208251953D0
      U( 47) =   -2.499999999999990D0
      U( 48) =   -2.499999999999975D0
      U( 49) =    0.204040692499482D0
      U( 50) =    0.011817791748035D0
      U( 51) =    0.192222900751447D0
      U( 52) =   -2.499999999999991D0
      U( 53) =   -2.499999999999990D0
      U( 54) =   -0.026041071738432D0
      U( 55) =    0.204040692499482D0
      U( 56) =   -2.499999999999992D0
      U( 57) =   -2.499999999999990D0
      U( 58) =   -0.026041071738434D0
      U( 59) =    0.000000000000000D0
      U( 60) =   -0.026041071738432D0
      U( 61) =   -2.499999999999991D0
      U( 62) =   -2.499999999999992D0
      U( 63) =    0.215858484880918D0
      U( 64) =    4.988182208251953D0
      U( 65) =   -2.499999999999990D0
      U( 66) =   -2.499999999999975D0
      U( 67) =    0.204040693132870D0
      U( 68) =    0.011817791748026D0
      U( 69) =    0.192222901384845D0
      U( 70) =   -2.499999999999991D0
      U( 71) =   -2.499999999999990D0
      U( 72) =   -0.026041071737961D0
      U( 73) =    0.204040693132870D0
      U( 74) =   -2.499999999999992D0
      U( 75) =   -2.499999999999990D0
      U( 76) =   -0.026041071737963D0
      U( 77) =    0.000000000000000D0
      U( 78) =   -0.026041071737961D0
      U( 79) =   -2.499999999999991D0
      U( 80) =   -2.499999999999992D0
      U( 81) =    4.999999999996546D0
      U( 82) =    4.999999999999970D0
      U( 83) =   -2.499999999999975D0
      U( 84) =   -2.499999999999975D0
      U( 85) =    4.999999999996516D0
      U( 86) =    0.000000000000000D0
      U( 87) =    4.999999999996516D0
      U( 88) =   -2.499999999999991D0
      U( 89) =   -2.499999999999975D0
      U( 90) =    0.000000000000000D0
      U( 91) =    4.999999999996516D0
      U( 92) =   -2.499999999999991D0
      U( 93) =   -2.499999999999975D0
      U( 94) =    0.215858481060569D0
      U( 95) =    4.988182208251953D0
      U( 96) =   -2.499999999999990D0
      U( 97) =   -2.499999999999975D0
      U( 98) =    0.204040689312522D0
      U( 99) =    0.011817791748023D0
      U(100) =    0.192222897564498D0
      U(101) =   -2.499999999999991D0
      U(102) =   -2.499999999999990D0
      U(103) =    4.734672533390068D0
      U(104) =    0.204040689312522D0
      U(105) =   -2.499999999999977D0
      U(106) =   -2.499999999999990D0
      U(107) =    4.734672533390062D0
      U(108) =    0.000000000000000D0
      U(109) =    4.734672533390068D0
      U(110) =   -2.499999999999991D0
      U(111) =   -2.499999999999977D0
      U(112) =    4.999999999996870D0
      U(113) =    4.999999999999972D0
      U(114) =   -2.499999999999975D0
      U(115) =   -2.499999999999975D0
      U(116) =    4.999999999996843D0
      U(117) =   -0.025968303070038D0
      U(118) =    4.999999999996843D0
      U(119) =   -2.499999999999992D0
      U(120) =   -2.499999999999975D0
      U(121) =   -0.025968303070040D0
      U(122) =    0.000000000000000D0
      U(123) =   -0.025968303070038D0
      U(124) =   -2.499999999999991D0
      U(125) =   -2.499999999999992D0
      U(126) =    4.999999999997699D0
      U(127) =    4.999999999999980D0
      U(128) =   -2.499999999999975D0
      U(129) =   -2.499999999999975D0
      U(130) =    4.999999999997678D0
      U(131) =    4.744923533081106D0
      U(132) =    4.999999999997678D0
      U(133) =   -2.499999999999977D0
      U(134) =   -2.499999999999975D0
      U(135) =    4.744923533081098D0
      U(136) =    0.000000000000000D0
      U(137) =    4.744923533081106D0
      U(138) =   -2.499999999999991D0
      U(139) =   -2.499999999999977D0
      U(140) =    0.000000000000000D0
      U(141) =    4.744923533081106D0
      U(142) =   -2.499999999999991D0
      U(143) =   -2.499999999999977D0
      U(144) =    0.215858484844162D0
      U(145) =    4.988182208251953D0
      U(146) =   -2.499999999999990D0
      U(147) =   -2.499999999999975D0
      U(148) =    0.204040693096114D0
      U(149) =    0.011817791748023D0
      U(150) =    0.192222901348091D0
      U(151) =   -2.499999999999991D0
      U(152) =   -2.499999999999990D0
      U(153) =    0.204040693096045D0
      U(154) =    0.204040693096107D0
      U(155) =   -2.499999999999990D0
      U(156) =   -2.499999999999990D0
      U(157) =    0.204040693096037D0
      U(158) =    0.000000000000000D0
      U(159) =    0.204040693096037D0
      U(160) =   -2.499999999999991D0
      U(161) =   -2.499999999999990D0
      U(162) =   -0.026017361873565D0
      U(163) =    0.204040693096114D0
      U(164) =   -2.499999999999992D0
      U(165) =   -2.499999999999990D0
      U(166) =   -0.026017361873568D0
      U(167) =   -0.026017590106916D0
      U(168) =   -0.026017361873565D0
      U(169) =   -2.499999999999992D0
      U(170) =   -2.499999999999992D0
      U(171) =   -0.026017590106918D0
      U(172) =    0.000000000000000D0
      U(173) =   -0.026017590106916D0
      U(174) =   -2.499999999999991D0
      U(175) =   -2.499999999999992D0
C ...
      CALL GCN(175,U,Y)
C ...
      DO 20 I=1,175
         Y(I+175)=U(I)
   20 CONTINUE
      DO 30 J=1,NEQ
         YP(J)=0.0D0
   30 CONTINUE
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
         OPEN(UNIT= IWRT,FILE= 'tba-das.dat',STATUS= 'UNKNOWN')
         OPEN(UNIT= NUMOUT,FILE= 'tba-das-n.dat',STATUS= 'UNKNOWN')
         RTOL  = 1.0D-5
         ATOL  = 1.7D-6
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
      INFO(5) = 0
      INFO(6) = 0
      IW(1) = ML       ! Used whenever INFO(6) = 1
      IW(2) = MU       ! Used whenever INFO(6) = 1
C
C ... Specify a maximum stepsize.
C
      INFO(7) = 0
      RW(2) = 10.0D0    ! Used whenever INFO(7) = 1
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
      DOUBLE PRECISION X(175),WRK(175)
      EXTERNAL FCN,GCN
C ...
      DO 10 I=1,175
  10  X(I)=Y(I+175)
C ...
      CALL FCN(175,T,X,WRK,IER)
      IF( IER.EQ.-1 ) RETURN
C ...
      DO 20 I=1,175
  20  GG(I)=-YP(I) + WRK(I)
C ...
      CALL GCN(175,X,WRK)
C ...
      DO 30 I=1,175
  30  GG(I+175)=Y(I)-WRK(I)
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
      DOUBLE PRECISION  RPAR(*),Y(*),YP(*),DGG(350,*)
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
      DOUBLE PRECISION EXACT(350),PREC(350),MAXVAL
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
 5       FORMAT(/'Two bit adding machine example problem for DASSL')
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
C ... COMPUTED ON CRAY C90 USING CRAY DOUBLE PRECISION
C ... SOLVING TRANSISTOR AMPLIFIER USING PSIDE
C ...
C ... USER INPUT:
C ...
C ... GIVE RELATIVE ERROR TOLERANCE: 1D-14
C ... GIVE ABSOLUTE ERROR TOLERANCE: 1D-14
C ...
C ...
C ... INTEGRATION CHARACTERISTICS:
C ...
C ...    NUMBER OF INTEGRATION STEPS       16061
C ...    NUMBER OF ACCEPTED STEPS          15824
C ...    NUMBER OF F EVALUATIONS          401944
C ...    NUMBER OF JACOBIAN EVALUATIONS      458
C ...    NUMBER OF LU DECOMPOSITIONS        4884
C ...
C ... CPU-TIME USED:                         182.44 SEC
C
C ... The reference solution was computed with RADAU5,
C ... atol=rtol=1d-5, h0=4d-5
C
C ... # steps                2518
C ... # steps accepted       2044
C ... # f-eval              21762
C ... # Jac-eval             1742
C ... # LU-decomp            2508
C
C ... The scd values reported in the section (Numerical solution of the
C ... problem) were computed by considering only the output signals.
C ... They are given by:
C
C ...      x(49)  = y(49+175)  = y(224)
C ...      x(130) = y(130+175) = y(305)
C ...      x(148) = y(148+175) = y(323)
C
      Y(  1) =   0.5329676348348893D-05
      Y(  2) =   0.6083199215822746D-03
      Y(  3) =  -0.5802788414936773D+00
      Y(  4) =  -0.1934408303298687D+00
      Y(  5) =   0.1924683208304517D+02
      Y(  6) =  -0.5143798198402260D-61
      Y(  7) =   0.2999386350402069D+01
      Y(  8) =  -0.3048568278060051D+00
      Y(  9) =  -0.5802788414936773D+00
      Y( 10) =  -0.5143798198402260D-61
      Y( 11) =   0.2999386350402069D+01
      Y( 12) =  -0.3048568278060051D+00
      Y( 13) =  -0.5794434527225972D+00
      Y( 14) =   0.7090602915275927D-02
      Y( 15) =   0.2870386839346849D+01
      Y( 16) =  -0.3201995162025645D+00
      Y( 17) =  -0.1934408303298687D+00
      Y( 18) =  -0.3831797456449475D+01
      Y( 19) =  -0.2992291455828521D+01
      Y( 20) =  -0.2883958687237743D+01
      Y( 21) =  -0.3048568278060051D+00
      Y( 22) =  -0.3201995162025645D+00
      Y( 23) =  -0.1408045046739925D+00
      Y( 24) =   0.1225221993879181D+00
      Y( 25) =  -0.2863915327854759D+00
      Y( 26) =  -0.3201995162025645D+00
      Y( 27) =   0.5727830655709517D+00
      Y( 28) =  -0.3568994249538271D-06
      Y( 29) =  -0.1408045061245250D+00
      Y( 30) =  -0.3048568278060051D+00
      Y( 31) =  -0.2863915327854759D+00
      Y( 32) =   0.1046493633015554D-06
      Y( 33) =   0.1194233059399546D-04
      Y( 34) =  -0.5803216340894781D+00
      Y( 35) =  -0.1934408303298687D+00
      Y( 36) =   0.7618503612049611D+01
      Y( 37) =  -0.1225225577378757D+00
      Y( 38) =   0.2877465395282166D+01
      Y( 39) =  -0.3048568278060051D+00
      Y( 40) =  -0.5803216340894781D+00
      Y( 41) =  -0.5144043562440774D-61
      Y( 42) =   0.2999987953020043D+01
      Y( 43) =  -0.3048568278060051D+00
      Y( 44) =  -0.5794087680278276D+00
      Y( 45) =   0.7090674536833933D-02
      Y( 46) =   0.2870484175935613D+01
      Y( 47) =  -0.3201875688249403D+00
      Y( 48) =  -0.1934408303298687D+00
      Y( 49) =  -0.1508928314544720D+01
      Y( 50) =  -0.2992897263001197D+01
      Y( 51) =  -0.2884653493511335D+01
      Y( 52) =  -0.3048568278060051D+00
      Y( 53) =  -0.3201875688249403D+00
      Y( 54) =  -0.2207815941707912D-01
      Y( 55) =   0.1224251495275525D+00
      Y( 56) =  -0.3020231631441768D+00
      Y( 57) =  -0.3201875688249403D+00
      Y( 58) =   0.6040463262883536D+00
      Y( 59) =  -0.1225235639587947D+00
      Y( 60) =  -0.1445997109340358D+00
      Y( 61) =  -0.3048568278060051D+00
      Y( 62) =  -0.3020231631441768D+00
      Y( 63) =   0.7090617704627438D-02
      Y( 64) =   0.2870406940405611D+01
      Y( 65) =  -0.3201970487035549D+00
      Y( 66) =  -0.1934408303298687D+00
      Y( 67) =  -0.7839517558869167D+01
      Y( 68) =  -0.2992292602180288D+01
      Y( 69) =  -0.2883977656734086D+01
      Y( 70) =  -0.3048568278060051D+00
      Y( 71) =  -0.3201970487035549D+00
      Y( 72) =  -0.2207815941706669D-01
      Y( 73) =   0.1225024418897568D+00
      Y( 74) =  -0.3020231631441784D+00
      Y( 75) =  -0.3201970487035548D+00
      Y( 76) =   0.6040463262883567D+00
      Y( 77) =  -0.1225235639587947D+00
      Y( 78) =  -0.1445997109340232D+00
      Y( 79) =  -0.3048568278060051D+00
      Y( 80) =  -0.3020231631441784D+00
      Y( 81) =   0.2827557455215537D-06
      Y( 82) =   0.3226083665770413D-04
      Y( 83) =  -0.5803201765692529D+00
      Y( 84) =  -0.1934408303298687D+00
      Y( 85) =   0.1349596894124011D+02
      Y( 86) =  -0.5144029989944615D-61
      Y( 87) =   0.2999967456407596D+01
      Y( 88) =  -0.3048568278060051D+00
      Y( 89) =  -0.5803201765692529D+00
      Y( 90) =  -0.5144029989944615D-61
      Y( 91) =   0.2999967456407596D+01
      Y( 92) =  -0.3048568278060051D+00
      Y( 93) =  -0.5794099530682210D+00
      Y( 94) =   0.7090627069218701D-02
      Y( 95) =   0.2870419659132725D+01
      Y( 96) =  -0.3201954889813194D+00
      Y( 97) =  -0.1934408303298687D+00
      Y( 98) =  -0.4564991969572258D+01
      Y( 99) =  -0.2992873132356891D+01
      Y(100) =  -0.2884572066660236D+01
      Y(101) =  -0.3048568278060051D+00
      Y(102) =  -0.3201954889813193D+00
      Y(103) =  -0.1408098389641259D+00
      Y(104) =   0.1224893544588554D+00
      Y(105) =  -0.2863908148025885D+00
      Y(106) =  -0.3201954889813194D+00
      Y(107) =   0.5727816296051769D+00
      Y(108) =  -0.3578850567906977D-06
      Y(109) =  -0.1408098404182732D+00
      Y(110) =  -0.3048568278060051D+00
      Y(111) =  -0.2863908148025885D+00
      Y(112) =  -0.5030971301117276D-05
      Y(113) =  -0.5741341370739308D-03
      Y(114) =  -0.5803636949068757D+00
      Y(115) =  -0.1934408303298687D+00
      Y(116) =   0.1342161027531811D+01
      Y(117) =   0.4675775943965132D+00
      Y(118) =   0.2878076723218615D+01
      Y(119) =  -0.3737359865057522D+00
      Y(120) =  -0.5803636949068773D+00
      Y(121) =   0.7474719730115046D+00
      Y(122) =  -0.1224897137980597D+00
      Y(123) =   0.4675903224882113D+00
      Y(124) =  -0.3048568278060051D+00
      Y(125) =  -0.3737359865057511D+00
      Y(126) =   0.1438859667693932D-04
      Y(127) =   0.1643466830238340D-02
      Y(128) =  -0.5802045234551009D+00
      Y(129) =  -0.1934408303298687D+00
      Y(130) =   0.1115322810290722D+02
      Y(131) =  -0.1786085072486350D+00
      Y(132) =  -0.2245965361045817D-02
      Y(133) =  -0.5675131482652329D+00
      Y(134) =  -0.5802045234550925D+00
      Y(135) =   0.1702539444795697D+01
      Y(136) =  -0.1225024418897565D+00
      Y(137) =   0.2699459271144151D+01
      Y(138) =  -0.3048568278060051D+00
      Y(139) =  -0.5675131482652329D+00
      Y(140) =  -0.1224897137980597D+00
      Y(141) =   0.2699471999235856D+01
      Y(142) =  -0.3048568278060051D+00
      Y(143) =  -0.5675131482652517D+00
      Y(144) =   0.7090737749608171D-02
      Y(145) =   0.2870570104712901D+01
      Y(146) =  -0.3201770176390597D+00
      Y(147) =  -0.1934408303298687D+00
      Y(148) =  -0.1189155590887148D+01
      Y(149) =  -0.2992881712346612D+01
      Y(150) =  -0.2884714042931089D+01
      Y(151) =  -0.3048568278060051D+00
      Y(152) =  -0.3201770176390597D+00
      Y(153) =  -0.2877042169326476D+01
      Y(154) =  -0.2877046623845863D+01
      Y(155) =  -0.3201777033645129D+00
      Y(156) =  -0.3201770176390597D+00
      Y(157) =   0.6403554067290257D+00
      Y(158) =  -0.1224897137980597D+00
      Y(159) =  -0.1449637037525297D-03
      Y(160) =  -0.3048568278060051D+00
      Y(161) =  -0.3201777033645128D+00
      Y(162) =  -0.1512233366530073D+00
      Y(163) =  -0.1567296590745326D-03
      Y(164) =  -0.3011638803836413D+00
      Y(165) =  -0.3201770176390597D+00
      Y(166) =   0.6023277607672893D+00
      Y(167) =  -0.1545119766856435D+00
      Y(168) =  -0.1512657482779617D+00
      Y(169) =  -0.3007435275083740D+00
      Y(170) =  -0.3011638803836447D+00
      Y(171) =   0.6014870550167479D+00
      Y(172) =  -0.5144056888603891D-61
      Y(173) =  -0.3199281323439955D-01
      Y(174) =  -0.3048568278060051D+00
      Y(175) =  -0.3007435275083740D+00
      Y(176) =   0.4998987079433206D+01
      Y(177) =   0.4999992063175263D+01
      Y(178) =  -0.2499999831773528D+01
      Y(179) =  -0.2499999999999975D+01
      Y(180) =   0.4998978196639293D+01
      Y(181) =  -0.8572996997335921D-61
      Y(182) =   0.4998977250670116D+01
      Y(183) =  -0.2499999999999991D+01
      Y(184) =  -0.2499999831773528D+01
      Y(185) =  -0.8572996997335921D-61
      Y(186) =   0.4998977250670116D+01
      Y(187) =  -0.2499999999999991D+01
      Y(188) =  -0.2500000136347459D+01
      Y(189) =   0.2160218554281297D+00
      Y(190) =   0.4988182249480750D+01
      Y(191) =  -0.2500000024216612D+01
      Y(192) =  -0.2499999999999975D+01
      Y(193) =   0.2042041839026698D+00
      Y(194) =   0.1182577025842485D-01
      Y(195) =   0.1923803845763885D+00
      Y(196) =  -0.2499999999999991D+01
      Y(197) =  -0.2500000024216612D+01
      Y(198) =  -0.2346741744566541D+00
      Y(199) =   0.2042036656465303D+00
      Y(200) =  -0.2499999823257662D+01
      Y(201) =  -0.2500000024216613D+01
      Y(202) =  -0.2346742463623664D+00
      Y(203) =  -0.5948323749230453D-06
      Y(204) =  -0.2346741768742084D+00
      Y(205) =  -0.2499999999999991D+01
      Y(206) =  -0.2499999823257662D+01
      Y(207) =   0.4999980114197910D+01
      Y(208) =   0.4999999843666627D+01
      Y(209) =  -0.2499999996784775D+01
      Y(210) =  -0.2499999999999975D+01
      Y(211) =   0.4999979939782303D+01
      Y(212) =  -0.7899378970256425D-07
      Y(213) =   0.4999979842706281D+01
      Y(214) =  -0.2499999999999991D+01
      Y(215) =  -0.2499999996784775D+01
      Y(216) =  -0.8573405937400700D-61
      Y(217) =   0.4999979921700072D+01
      Y(218) =  -0.2499999999999991D+01
      Y(219) =  -0.2500000002605499D+01
      Y(220) =   0.2158597056211766D+00
      Y(221) =   0.4988182207952476D+01
      Y(222) =  -0.2500000000353345D+01
      Y(223) =  -0.2499999999999975D+01
      Y(224) =   0.2040419147264534D+00
      Y(225) =   0.1181783478030806D-01
      Y(226) =   0.1922241172634104D+00
      Y(227) =  -0.2499999999999991D+01
      Y(228) =  -0.2500000000353345D+01
      Y(229) =  -0.3679693236179853D-01
      Y(230) =   0.2040419158792541D+00
      Y(231) =  -0.2499999772012851D+01
      Y(232) =  -0.2500000000353345D+01
      Y(233) =  -0.3679622453612054D-01
      Y(234) =  -0.1756028654736632D-05
      Y(235) =  -0.3679533432072315D-01
      Y(236) =  -0.2499999999999991D+01
      Y(237) =  -0.2499999772012851D+01
      Y(238) =   0.2159883644910842D+00
      Y(239) =   0.4988182235659396D+01
      Y(240) =  -0.2500000020901020D+01
      Y(241) =  -0.2499999999999975D+01
      Y(242) =   0.2041706683167043D+00
      Y(243) =   0.1182385967214655D-01
      Y(244) =   0.1923487687491322D+00
      Y(245) =  -0.2499999999999991D+01
      Y(246) =  -0.2500000020901020D+01
      Y(247) =  -0.3679693236177788D-01
      Y(248) =   0.2041707364829278D+00
      Y(249) =  -0.2499999772012851D+01
      Y(250) =  -0.2500000020901019D+01
      Y(251) =  -0.3679622453609997D-01
      Y(252) =  -0.1756028654731466D-05
      Y(253) =  -0.3679533432070254D-01
      Y(254) =  -0.2499999999999991D+01
      Y(255) =  -0.2499999772012850D+01
      Y(256) =   0.4999946292289309D+01
      Y(257) =   0.4999999589090829D+01
      Y(258) =  -0.2499999989270124D+01
      Y(259) =  -0.2499999999999975D+01
      Y(260) =   0.4999945821029733D+01
      Y(261) =  -0.8573383316572590D-61
      Y(262) =   0.4999945760679328D+01
      Y(263) =  -0.2499999999999991D+01
      Y(264) =  -0.2499999989270124D+01
      Y(265) =  -0.8573383316572590D-61
      Y(266) =   0.4999945760679329D+01
      Y(267) =  -0.2499999999999991D+01
      Y(268) =  -0.2500000008700431D+01
      Y(269) =   0.2159672042326096D+00
      Y(270) =   0.4988182257671780D+01
      Y(271) =  -0.2500000009368004D+01
      Y(272) =  -0.2499999999999975D+01
      Y(273) =   0.2041494924505785D+00
      Y(274) =   0.1182393376824584D-01
      Y(275) =   0.1923257099293294D+00
      Y(276) =  -0.2499999999999991D+01
      Y(277) =  -0.2500000009368003D+01
      Y(278) =  -0.2346830649402098D+00
      Y(279) =   0.2041489240980933D+00
      Y(280) =  -0.2499999822769364D+01
      Y(281) =  -0.2500000009368002D+01
      Y(282) =  -0.2346831370442538D+00
      Y(283) =  -0.5964750946511629D-06
      Y(284) =  -0.2346830673637887D+00
      Y(285) =  -0.2499999999999991D+01
      Y(286) =  -0.2499999822769363D+01
      Y(287) =   0.5000956289101704D+01
      Y(288) =   0.5000007783825420D+01
      Y(289) =  -0.2500000106879903D+01
      Y(290) =  -0.2499999999999975D+01
      Y(291) =   0.5000964674053872D+01
      Y(292) =   0.9834666589775638D+00
      Y(293) =   0.5000965207014408D+01
      Y(294) =  -0.2500000011950726D+01
      Y(295) =  -0.2500000106879903D+01
      Y(296) =   0.9834666825678713D+00
      Y(297) =  -0.3054618845645094D-07
      Y(298) =   0.9834666965975977D+00
      Y(299) =  -0.2499999999999991D+01
      Y(300) =  -0.2500000011950726D+01
      Y(301) =   0.4997262436706446D+01
      Y(302) =   0.4999977567095746D+01
      Y(303) =  -0.2499999724698131D+01
      Y(304) =  -0.2499999999999975D+01
      Y(305) =   0.4997238455712048D+01
      Y(306) =   0.4703283828639512D+01
      Y(307) =   0.4997221398452132D+01
      Y(308) =  -0.2499999197373328D+01
      Y(309) =  -0.2499999724698136D+01
      Y(310) =   0.4703273936740555D+01
      Y(311) =  -0.6816622292221013D-07
      Y(312) =   0.4703269453557049D+01
      Y(313) =  -0.2499999999999991D+01
      Y(314) =  -0.2499999197373328D+01
      Y(315) =  -0.3054618845625812D-07
      Y(316) =   0.4703269491177070D+01
      Y(317) =  -0.2499999999999991D+01
      Y(318) =  -0.2499999197373338D+01
      Y(319) =   0.2157164867589083D+00
      Y(320) =   0.4988182098364394D+01
      Y(321) =  -0.2500000001646564D+01
      Y(322) =  -0.2499999999999975D+01
      Y(323) =   0.2038985905095614D+00
      Y(324) =   0.1180963378537931D-01
      Y(325) =   0.1920890828112512D+00
      Y(326) =  -0.2499999999999991D+01
      Y(327) =  -0.2500000001646564D+01
      Y(328) =   0.2039079144284985D+00
      Y(329) =   0.2039004902295213D+00
      Y(330) =  -0.2500000004493590D+01
      Y(331) =  -0.2500000001646564D+01
      Y(332) =   0.2039079021505135D+00
      Y(333) =  -0.3054618845623932D-07
      Y(334) =   0.2039078862776569D+00
      Y(335) =  -0.2499999999999991D+01
      Y(336) =  -0.2500000004493590D+01
      Y(337) =  -0.4788940197110148D-01
      Y(338) =   0.2038882763521206D+00
      Y(339) =  -0.2499999353490320D+01
      Y(340) =  -0.2500000001646564D+01
      Y(341) =  -0.4789765786972171D-01
      Y(342) =  -0.5331577724007001D-01
      Y(343) =  -0.4790539656059958D-01
      Y(344) =  -0.2499999201728321D+01
      Y(345) =  -0.2499999353490320D+01
      Y(346) =  -0.5331888562403928D-01
      Y(347) =  -0.8573428147674035D-61
      Y(348) =  -0.5332135539066590D-01
      Y(349) =  -0.2499999999999991D+01
      Y(350) =  -0.2499999201728321D+01

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
C    +           'THE FUNCTION GETTIM'
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
C
      SUBROUTINE FCN(N,X,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the two-bit adder
C computing the sum of two two-bit numbers and one carry-in:
C
C   A1*2+A0 + B1*2+B0 + CIN = C*2^2 + S1*2 + S0
C
C Input signals: V1  === A0_INV
C                V2  === B0_INV
C                V3  === A1_INV
C                V4  === B1_INV
C                CIN === CARRY_IN
C
C Output signals: SO === AT NODE 49
C                 S1 === AT NODE 130
C                  C  === AT NODE 148
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   X    Time-point at which FCN evaluated
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      integer*4 N,ierr
      double precision IDS,IBS,IBD
      DIMENSION Y(N),F(N)
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
      EXTERNAL IDS, IBS, IBD

C
C --- Evaluating the input signals at time point X
C
c      print *,X,Y(176)
      CALL PULSE(X,V1,V1D,0.d0,5.d0,0.d0,5.d0,5.d0,5.d0,20.d0)
      CALL PULSE(X,V2,V2D,0.d0,5.d0,10.d0,5.d0,15.d0,5.d0,40.d0)
      CALL PULSE(X,V3,V3D,0.d0,5.d0,30.d0,5.d0,35.d0,5.d0,80.d0)
      CALL PULSE(X,V4,V4D,0.d0,5.d0,70.d0,5.d0,75.d0,5.d0,160.d0)
      CALL PULSE(X,CIN,CIND,0.d0,5.d0,150.d0,5.d0,155.d0,5.d0,320.d0)

C
C ---
C --- Right-hand side of the two-bit adder: The ten logical subcircuits
C ---
C

C
C --- NOR-gate 1: nodes 1 -- 13
C
      CALL NOR(N,1,V1,V2,V1D,V2D,Y,F,ierr)

C
C --- ANDOI-gate 1: nodes 14 -- 31
C

      CALL ANDOI(N,14,Y(5),V2,V1,0.d0,V2D,V1D,Y,F,ierr)

C
C --- NOR-gate 2: nodes 32 -- 44
C

      CALL NOR(N,32,Y(18),CIN,0.d0,CIND,Y,F,ierr)

C
C --- ANDOI-gate 2: nodes 45-- 62
C

      CALL ANDOI(N,45,Y(36),CIN,Y(18),0.d0,CIND,0.d0,Y,F,ierr)

C
C --- ANDOI-gate 3: nodes 63-- 80
C

      CALL ANDOI(N,63,Y(5),CIN,Y(18),0.d0,CIND,0.d0,Y,F,ierr)

C
C --- NOR-gate 3: nodes 81 -- 93
C

      CALL NOR(N,81,V3,V4,V3D,V4D,Y,F,ierr)

C
C --- ANDOI-gate 4: nodes 94 -- 111
C

      CALL ANDOI(N,94,Y(85),V4,V3,0.d0,V4D,V3D,Y,F,ierr)

C
C --- NAND-gate: nodes 112 -- 125
C

      CALL NAND(N,112,Y(67),Y(98),0.d0,0.d0,Y,F,ierr)

C
C --- ORANI-gate 1: nodes 126 -- 143
C

      CALL ORANI(N,126,Y(116),Y(67),Y(98),0.d0,0.d0,0.d0,Y,F,ierr)

C
C --- ANDOI-gate 5 (ANDOI-gate with capacitive coupling
C ---               of result node): nodes 144 -- 161
C

      CALL ANDOIP(N,144,Y(85),Y(5),Y(98),0.d0,0.d0,0.d0,Y,F,ierr)

C
C --- Three additional enhancement transistors in series
C --- First transistor:  Internal nodes 162 -- 165
C ---                    DRAIN = node 148 (First node of ANDOI-gate 5)
C ---                    GATE = node 98 (First node of ANDOI-gate 4)
C ---                    SOURCE = node  166
C --- Second transistor: Internal nodes 167 -- 171
C ---                    DRAIN = node 166
C ---                    GATE = node 18 (First node of ANDOI-gate 1)
C ---                    SOURCE = node 171
C --- Third transistor:  Internal nodes 172 -- 175
C ---                    DRAIN = node 171
C ---                    GATE = CIN (carry_in voltage)
C ---                    Source = MASS
C
      F(162)=-(Y(162)-Y(166))/RGS - IDS(3,Y(163)-Y(162),Y(98)-Y(162),
     *         Y(164)-Y(166),Y(98)-Y(163),Y(165)-Y(148),ierr)
      F(163)=-(Y(163)-Y(148))/RGD + IDS(3,Y(163)-Y(162),Y(98)-Y(162),
     *         Y(164)-Y(166),Y(98)-Y(163),Y(165)-Y(148),ierr)
      F(164)=-(Y(164)-VBB)/RBS +IBS(Y(164)-Y(166))
      F(165)=-(Y(165)-VBB)/RBD +IBD(Y(165)-Y(148))
      F(166)=-IBS(Y(164)-Y(166))-(Y(166)-Y(162))/RGS-
     *        IBD(Y(170)-Y(166))-(Y(166)-Y(168))/RGD
      F(167)=-(Y(167)-Y(171))/RGS-IDS(3,Y(168)-Y(167),Y(18)-Y(167),
     *         Y(169)-Y(171),Y(18)-Y(168),Y(170)-Y(166),ierr)
      F(168)=-(Y(168)-Y(166))/RGD+IDS(3,Y(168)-Y(167),Y(18)-Y(167),
     *         Y(169)-Y(171),Y(18)-Y(168),Y(170)-Y(166),ierr)
      F(169)=-(Y(169)-VBB)/RBS + IBS(Y(169)-Y(171))
      F(170)=-(Y(170)-VBB)/RBD + IBD(Y(170)-Y(166))
      F(171)=-IBS(Y(169)-Y(171))-(Y(171)-Y(167))/RGS-
     *        IBD(Y(175)-Y(171))-(Y(171)-Y(173))/RGD
      F(172)= CGS*CIND-Y(172)/RGS - IDS(3,Y(173)-Y(172),
     *        CIN-Y(172),Y(174),CIN-Y(173),Y(175)-Y(171),ierr)
      F(173)= CGD*CIND-(Y(173)-Y(171))/RGD + IDS(3,Y(173)-Y(172),
     *        CIN-Y(172),Y(174),CIN-Y(173),Y(175)-Y(171),ierr)
      F(174)=-(Y(174)-VBB)/RBS + IBS(Y(174))
      F(175)=-(Y(175)-VBB)/RBD + IBD(Y(175)-Y(171))

       if(ierr.eq.-1)return

       RETURN
       END

      SUBROUTINE PULSE(X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD)
C ---------------------------------------------------------------------------
C
C Evaluating input signal at time point X
C
C Structure of input signal:
C
C                -----------------------                       HIGH
C               /                       \
C              /                         \
C             /                           \
C            /                             \
C           /                               \
C          /                                 \
C         /                                   \
C        /                                     \
C  ------                                       ---------      LOW
C
C |DELAY|   T1  |         T2           |   T3  |
C |          P     E     R     I     O     D            |
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   X                      Time-point at which input signal is evaluated
C   LOW                    Low-level of input signal
C   HIGH                   High-level of input signal
C   DELAY,T1,T2,T3, PERIOD Parameters to specify signal structure
C
C Output parameter:
C   VIN    Voltage of input signal at time point X
C   VIND   Derivative of VIN at time point X
C
C ---------------------------------------------------------------------------

      double precision X,VIN,VIND,LOW,HIGH,DELAY,T1,T2,T3,PERIOD,TIME

      TIME = DMOD(X,PERIOD)

      IF (TIME.GT.(DELAY+T1+T2+T3)) THEN
        VIN = LOW
        VIND= 0.d0
      ELSE IF (TIME.GT.(DELAY+T1+T2)) THEN
        VIN = ((HIGH-LOW)/T3)*(DELAY+T1+T2+T3-TIME) + LOW
        VIND= -((HIGH-LOW)/T3)
      ELSE IF (TIME.GT.(DELAY+T1)) THEN
        VIN = HIGH
        VIND= 0.d0
      ELSE IF (TIME.GT.DELAY) THEN
        VIN = ((HIGH-LOW)/T1)*(TIME-DELAY) + LOW
        VIND= ((HIGH-LOW)/T1)
      ELSE
        VIN = LOW
        VIND=0.d0
      END IF

      RETURN
      END

      SUBROUTINE NOR(N,I,U1,U2,U1D,U2D,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the NOR-gate:
C                   NOT (U1 OR U2)
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   I    Number of first node
C   U1   Voltage of first input signal
C   U2   Voltage of second input signal
C   U1D  Derivative of U1
C   U2D  Derivative of U2
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side of NOR-gate
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I,ierr
        double precision IDS,IBS,IBD
        DIMENSION Y(N),F(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
        EXTERNAL IDS, IBS, IBD

        F(I)=-(Y(I)-Y(I+4))/RGS-IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *        Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+1)=-(Y(I+1)-VDD)/RGD+IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *         Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+2)=-(Y(I+2)-VBB)/RBS + IBS(Y(I+2)-Y(I+4))
        F(I+3)=-(Y(I+3)-VBB)/RBD + IBD(Y(I+3)-VDD)
C
C --- Result node of NOR-gate: Node I+4
C
        F(I+4)=-(Y(I+4)-Y(I))/RGS-IBS(Y(I+2)-Y(I+4))-
     *         (Y(I+4)-Y(I+6))/RGD-IBD(Y(I+8)-Y(I+4))-
     *         (Y(I+4)-Y(I+10))/RGD-IBD(Y(I+12)-Y(I+4))
        F(I+5)=CGS*U1D-Y(I+5)/RGS-IDS(1,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+6)=CGD*U1D-(Y(I+6)-Y(I+4))/RGD+IDS(1,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+7)=-(Y(I+7)-VBB)/RBS + IBS(Y(I+7))
        F(I+8)=-(Y(I+8)-VBB)/RBD + IBD(Y(I+8)-Y(I+4))
        F(I+9)=CGS*U2D-Y(I+9)/RGS-IDS(1,Y(I+10)-Y(I+9),
     *         U2-Y(I+9),Y(I+11),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+10)=CGD*U2D-(Y(I+10)-Y(I+4))/RGD+IDS(1,Y(I+10)-Y(I+9),
     *          U2-Y(I+9),Y(I+11),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+11)=-(Y(I+11)-VBB)/RBS + IBS(Y(I+11))
        F(I+12)=-(Y(I+12)-VBB)/RBD + IBD(Y(I+12)-Y(I+4))

        if(ierr.eq.-1)return

        RETURN
        END

      SUBROUTINE ANDOIP(N,I,U1,U2,U3,U1D,U2D,U3D,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the ANDOI-gate with capacitive
C coupling at result node (used for two-bit adder)
C                   NOT ( U1 OR ( U2 AND U3) )
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   I    Number of first node
C   U1   Voltage of first input signal
C   U2   Voltage of second input signal
C   U3   Voltage of third input signal
C   U1D  Derivative of U1
C   U2D  Derivative of U2
C   U3D  Derivative of U3
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side of ANDOI-gate
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I,ierr
        double precision IDS,IBS,IBD
        DIMENSION Y(N),F(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
        EXTERNAL IDS, IBS, IBD

        F(I)=-(Y(I)-Y(I+4))/RGS-IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *       Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+1)=-(Y(I+1)-VDD)/RGD+IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *         Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+2)=-(Y(I+2)-VBB)/RBS + IBS(Y(I+2)-Y(I+4))
        F(I+3)=-(Y(I+3)-VBB)/RBD + IBD(Y(I+3)-VDD)
C
C --- Result node of ANDOI-gate: Node I+4
C
        F(I+4)=-(Y(I+4)-Y(I))/RGS-IBS(Y(I+2)-Y(I+4))-
     *         (Y(I+4)-Y(I+6))/RGD-IBD(Y(I+8)-Y(I+4))-
     *         (Y(I+4)-Y(I+10))/RGD-IBD(Y(I+12)-Y(I+4))-
     *         (Y(I+4)-Y(163))/RGD-IBD(Y(165)-Y(I+4))
        F(I+5)=CGS*U1D-Y(I+5)/RGS-IDS(1,Y(I+6)-Y(I+5),U1-Y(I+5),
     *         Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+6)=CGD*U1D-(Y(I+6)-Y(I+4))/RGD+IDS(1,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+7)=-(Y(I+7)-VBB)/RBS + IBS(Y(I+7))
        F(I+8)=-(Y(I+8)-VBB)/RBD + IBD(Y(I+8)-Y(I+4))
        F(I+9)=CGS*U2D-(Y(I+9)-Y(I+13))/RGS-IDS(2,Y(I+10)-Y(I+9),
     *         U2-Y(I+9),Y(I+11)-Y(I+13),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+10)=CGD*U2D-(Y(I+10)-Y(I+4))/RGD+IDS(2,Y(I+10)-Y(I+9),
     *         U2-Y(I+9),Y(I+11)-Y(I+13),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+11)=-(Y(I+11)-VBB)/RBS + IBS(Y(I+11)-Y(I+13))
        F(I+12)=-(Y(I+12)-VBB)/RBD + IBD(Y(I+12)-Y(I+4))
C
C --- Coupling node of ANDOI-gate: Node I+13
C
        F(I+13)=-(Y(I+13)-Y(I+9))/RGS-IBS(Y(I+11)-Y(I+13))-
     *           (Y(I+13)-Y(I+15))/RGD-IBD(Y(I+17)-Y(I+13))
        F(I+14)=CGS*U3D-Y(I+14)/RGS-IDS(2,Y(I+15)-Y(I+14),
     *          U3-Y(I+14),Y(I+16),U3-Y(I+15),Y(I+17)-Y(I+13),ierr)
        F(I+15)=CGD*U3D-(Y(I+15)-Y(I+13))/RGD+IDS(2,Y(I+15)-Y(I+14),
     *          U3-Y(I+14),Y(I+16),U3-Y(I+15),Y(I+17)-Y(I+13),ierr)
        F(I+16)=-(Y(I+16)-VBB)/RBS+IBS(Y(I+16))
        F(I+17)=-(Y(I+17)-VBB)/RBD+IBD(Y(I+17)-Y(I+13))

        if(ierr.eq.-1)return

        RETURN
        END

       SUBROUTINE ANDOI(N,I,U1,U2,U3,U1D,U2D,U3D,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the ANDOI-gate
C                   NOT ( U1 OR ( U2 AND U3) )
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   I    Number of first node
C   U1   Voltage of first input signal
C   U2   Voltage of second input signal
C   U3   Voltage of third input signal
C   U1D  Derivative of U1
C   U2D  Derivative of U2
C   U3D  Derivative of U3
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side of ANDOI-gate
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        double precision IDS,IBS,IBD
        integer*4 N,I,ierr
        DIMENSION Y(N),F(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
        EXTERNAL IDS, IBS, IBD

        F(I)=-(Y(I)-Y(I+4))/RGS-IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *       Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+1)=-(Y(I+1)-VDD)/RGD+IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *         Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+2)=-(Y(I+2)-VBB)/RBS + IBS(Y(I+2)-Y(I+4))
        F(I+3)=-(Y(I+3)-VBB)/RBD + IBD(Y(I+3)-VDD)
C
C --- Result node of ANDOI-gate: Node I+4
C
        F(I+4)=-(Y(I+4)-Y(I))/RGS-IBS(Y(I+2)-Y(I+4))-
     *         (Y(I+4)-Y(I+6))/RGD-IBD(Y(I+8)-Y(I+4))-
     *         (Y(I+4)-Y(I+10))/RGD-IBD(Y(I+12)-Y(I+4))
        F(I+5)=CGS*U1D-Y(I+5)/RGS-IDS(1,Y(I+6)-Y(I+5),U1-Y(I+5),
     *         Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+6)=CGD*U1D-(Y(I+6)-Y(I+4))/RGD+IDS(1,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+7)=-(Y(I+7)-VBB)/RBS + IBS(Y(I+7))
        F(I+8)=-(Y(I+8)-VBB)/RBD + IBD(Y(I+8)-Y(I+4))
        F(I+9)=CGS*U2D-(Y(I+9)-Y(I+13))/RGS-IDS(2,Y(I+10)-Y(I+9),
     *         U2-Y(I+9),Y(I+11)-Y(I+13),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+10)=CGD*U2D-(Y(I+10)-Y(I+4))/RGD+IDS(2,Y(I+10)-Y(I+9),
     *         U2-Y(I+9),Y(I+11)-Y(I+13),U2-Y(I+10),Y(I+12)-Y(I+4),ierr)
        F(I+11)=-(Y(I+11)-VBB)/RBS + IBS(Y(I+11)-Y(I+13))
        F(I+12)=-(Y(I+12)-VBB)/RBD + IBD(Y(I+12)-Y(I+4))
C
C --- Coupling node of ANDOI-gate: Node I+13
C
        F(I+13)=-(Y(I+13)-Y(I+9))/RGS-IBS(Y(I+11)-Y(I+13))-
     *           (Y(I+13)-Y(I+15))/RGD-IBD(Y(I+17)-Y(I+13))
        F(I+14)=CGS*U3D-Y(I+14)/RGS-IDS(2,Y(I+15)-Y(I+14),
     *          U3-Y(I+14),Y(I+16),U3-Y(I+15),Y(I+17)-Y(I+13),ierr)
        F(I+15)=CGD*U3D-(Y(I+15)-Y(I+13))/RGD+IDS(2,Y(I+15)-Y(I+14),
     *          U3-Y(I+14),Y(I+16),U3-Y(I+15),Y(I+17)-Y(I+13),ierr)
        F(I+16)=-(Y(I+16)-VBB)/RBS+IBS(Y(I+16))
        F(I+17)=-(Y(I+17)-VBB)/RBD+IBD(Y(I+17)-Y(I+13))

        if(ierr.eq.-1)return

        RETURN
        END

      SUBROUTINE NAND(N,I,U1,U2,U1D,U2D,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the NAND-gate:
C                   NOT (U1 AND U2)
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   I    Number of first node
C   U1   Voltage of first input signal
C   U2   Voltage of second input signal
C   U1D  Derivative of U1
C   U2D  Derivative of U2
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side of NAND-gate
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        double precision IDS,IBS,IBD
        integer*4 N,I,ierr
        DIMENSION Y(N),F(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
        EXTERNAL IDS, IBS, IBD

        F(I)=-(Y(I)-Y(I+4))/RGS-IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *        Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+1)=-(Y(I+1)-VDD)/RGD+IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *         Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+2)=-(Y(I+2)-VBB)/RBS + IBS(Y(I+2)-Y(I+4))
        F(I+3)=-(Y(I+3)-VBB)/RBD + IBD(Y(I+3)-VDD)
C
C --- Result node of NAND-gate: Node I+4
C
        F(I+4)=-(Y(I+4)-Y(I))/RGS-IBS(Y(I+2)-Y(I+4))-
     *         (Y(I+4)-Y(I+6))/RGD-IBD(Y(I+8)-Y(I+4))
        F(I+5)=CGS*U1D-(Y(I+5)-Y(I+9))/RGS-IDS(2,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7)-Y(I+9),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+6)=CGD*U1D-(Y(I+6)-Y(I+4))/RGD+IDS(2,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7)-Y(I+9),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+7)=-(Y(I+7)-VBB)/RBS + IBS(Y(I+7)-Y(I+9))
        F(I+8)=-(Y(I+8)-VBB)/RBD + IBD(Y(I+8)-Y(I+4))
C
C --- Coupling node of NAND-gate: Node I+9
C
        F(I+9)=-(Y(I+9)-Y(I+5))/RGS-IBS(Y(I+7)-Y(I+9))-
     *         (Y(I+9)-Y(I+11))/RGD-IBD(Y(I+13)-Y(I+9))
        F(I+10)=CGS*U2D-Y(I+10)/RGS-IDS(2,Y(I+11)-Y(I+10),
     *          U2-Y(I+10),Y(I+12),U2-Y(I+11),Y(I+13)-Y(I+9),ierr)
        F(I+11)=CGD*U2D-(Y(I+11)-Y(I+9))/RGD+IDS(2,Y(I+11)-Y(I+10),
     *          U2-Y(I+10),Y(I+12),U2-Y(I+11),Y(I+13)-Y(I+9),ierr)
        F(I+12)=-(Y(I+12)-VBB)/RBS + IBS(Y(I+12))
        F(I+13)=-(Y(I+13)-VBB)/RBD + IBD(Y(I+13)-Y(I+9))

        if(ierr.eq.-1)return

        RETURN
        END

      SUBROUTINE ORANI(N,I,U1,U2,U3,U1D,U2D,U3D,Y,F,ierr)
C ---------------------------------------------------------------------------
C
C Right-hand side (static currents) of the ORANI-gate
C                   NOT ( U1 AND ( U2 OR U3) )
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   I    Number of first node
C   U1   Voltage of first input signal
C   U2   Voltage of second input signal
C   U3   Voltage of third input signal
C   U1D  Derivative of U1
C   U2D  Derivative of U2
C   U3D  Derivative of U3
C   Y    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   F    Computed right-hand side of ORANI-gate
C
C External references: Currents of MOS-model due to Shichman and Hodges
C IDS      Function evaluating the drain current
C IBS, IBD Currents of pn-junction
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        double precision IDS,IBS,IBD
        integer*4 N,I,ierr
        DIMENSION Y(N),F(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
        EXTERNAL IDS, IBS, IBD

        F(I)=-(Y(I)-Y(I+4))/RGS-IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *       Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+1)=-(Y(I+1)-VDD)/RGD+IDS(0,Y(I+1)-Y(I),Y(I+4)-Y(I),
     *         Y(I+2)-Y(I+4),Y(I+4)-Y(I+1),Y(I+3)-VDD,ierr)
        F(I+2)=-(Y(I+2)-VBB)/RBS + IBS(Y(I+2)-Y(I+4))
        F(I+3)=-(Y(I+3)-VBB)/RBD + IBD(Y(I+3)-VDD)
C
C --- Result node of ORANI-gate: Node I+4
C
        F(I+4)=-(Y(I+4)-Y(I))/RGS-IBS(Y(I+2)-Y(I+4))-
     *         (Y(I+4)-Y(I+6))/RGD-IBD(Y(I+8)-Y(I+4))
        F(I+5)=CGS*U1D-(Y(I+5)-Y(I+9))/RGS-
     *         IDS(2,Y(I+6)-Y(I+5),U1-Y(I+5),Y(I+7)-Y(I+9),
     *         U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+6)=CGD*U1D-(Y(I+6)-Y(I+4))/RGD+IDS(2,Y(I+6)-Y(I+5),
     *         U1-Y(I+5),Y(I+7)-Y(I+9),U1-Y(I+6),Y(I+8)-Y(I+4),ierr)
        F(I+7)=-(Y(I+7)-VBB)/RBS + IBS(Y(I+7)-Y(I+9))
        F(I+8)=-(Y(I+8)-VBB)/RBD + IBD(Y(I+8)-Y(I+4))
C
C --- Coupling node of ORANI-gate: Node I+9
C
        F(I+9)=-(Y(I+9)-Y(I+5))/RGS-IBS(Y(I+7)-Y(I+9))-
     *         (Y(I+9)-Y(I+11))/RGD-IBD(Y(I+13)-Y(I+9))-
     *         (Y(I+9)-Y(I+15))/RGD-IBD(Y(I+17)-Y(I+9))
        F(I+10)=CGS*U2D-Y(I+10)/RGS-IDS(2,Y(I+11)-Y(I+10),
     *          U2-Y(I+10),Y(I+12),U2-Y(I+11),Y(I+13)-Y(I+9),ierr)
         F(I+11)=CGD*U2D-(Y(I+11)-Y(I+9))/RGD+IDS(2,Y(I+11)-Y(I+10),
     *          U2-Y(I+10),Y(I+12),U2-Y(I+11),Y(I+13)-Y(I+9),ierr)
        F(I+12)=-(Y(I+12)-VBB)/RBS + IBS(Y(I+12))
        F(I+13)=-(Y(I+13)-VBB)/RBD + IBD(Y(I+13)-Y(I+9))
        F(I+14)=CGS*U3D-Y(I+14)/RGS-IDS(2,Y(I+15)-Y(I+14),U3-Y(I+14),
     *          Y(I+16),U3-Y(I+15),Y(I+17)-Y(I+9),ierr)
        F(I+15)=CGD*U3D-(Y(I+15)-Y(I+9))/RGD+
     *        IDS(2,Y(I+15)-Y(I+14),U3-Y(I+14),Y(I+16),U3-Y(I+15),
     *        Y(I+17)-Y(I+9),ierr)
        F(I+16)=-(Y(I+16)-VBB)/RBS + IBS(Y(I+16))
        F(I+17)=-(Y(I+17)-VBB)/RBD + IBD(Y(I+17)-Y(I+9))

        if(ierr.eq.-1)return

        RETURN
        END

        SUBROUTINE GCN(N,U,G)
C ---------------------------------------------------------------------------
C
C Charge-function of two-bit adder
C computing the sum of two two-bit numbers and one carry-in:
C
C   A1*2+A0 + B1*2+B0 + CIN = C*2^2 + S1*2 + S0
C
C Input signals: V1  === A0_INV
C                V2  === B0_INV
C                V3  === A1_INV
C                V4  === B1_INV
C                CIN === CARRY_IN
C
C Output signals: SO === AT NODE 49
C                 S1 === AT NODE 130
C                  C  === AT NODE 148
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   U    State-vector (only node potentials in reduced model)
C
C Output parameter:
C   G    Computed charge function (only function of node potentials!)
C
C External reference:
C   CBDBS Function evaluating the voltage-dependent bulk-capacitance
C         in the MOS-model due to Shichman and Hodges
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I
        DIMENSION U(N),G(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

        EXTERNAL CBDBS

        DO 10 I=1,N
            G(I)=0.d0
 10     CONTINUE

C
C ---
C --- Charge-function of the two-bit adder: The ten logical subcircuits
C ---
C

C
C --- NOR-gate 1: nodes 1 -- 13
C
        CALL DNOR(N,U,1,G)

C
C --- ANDOI-gate 1: nodes 14 -- 31
C
        CALL DANDOI(N,U,14,G)

C
C --- NOR-gate 2: nodes 32 -- 44
C
        CALL DNOR(N,U,32,G)

C
C --- ANDOI-gate 2: nodes 45-- 62
C
        CALL DANDOI(N,U,45,G)

C
C --- ANDOI-gate 3: nodes 63-- 80
C
        CALL DANDOI(N,U,63,G)

C
C --- NOR-gate 3: nodes 81 -- 93
C
        CALL DNOR(N,U,81,G)

C
C --- ANDOI-gate 4: nodes 94 -- 111
C
        CALL DANDOI(N,U,94,G)
C
C --- NAND-gate: nodes 112 -- 125
C
        CALL DNAND(N,U,112,G)

C
C --- ORANI-gate 1: nodes 126 -- 143
C
        CALL DORANI(N,U,126,G)

C
C --- ANDOI-gate 5 (ANDOI-gate with capacitive coupling
C ---               of result node): nodes 144 -- 161
C
        CALL DANDOI(N,U,144,G)

C
C --- Capacitive coupling result node nor-gate 1
C
        G(5)=G(5)+CGS*(U(5)-U(19))+CGD*(U(5)-U(20))+
     *           CGS*(U(5)-U(68))+CGD*(U(5)-U(69))+
     *           CGS*(U(5)-U(153))+CGD*(U(5)-U(154))
        G(19)=G(19)-CGS*U(5)
        G(20)=G(20)-CGD*U(5)
        G(68)=G(68)-CGS*U(5)
        G(69)=G(69)-CGD*U(5)
        G(153)=G(153)-CGS*U(5)
        G(154)=G(154)-CGD*U(5)

C
C --- Capacitive coupling result node andoi-gate 1
C
        G(18)=G(18)+CGS*(U(18)-U(37))+CGD*(U(18)-U(38))+
     *             CGS*(U(18)-U(59))+CGD*(U(18)-U(60))+
     *             CGS*(U(18)-U(77))+CGD*(U(18)-U(78))+
     *             CGS*(U(18)-U(167))+CGD*(U(18)-U(168))
        G(37)=G(37)-CGS*U(18)
        G(38)=G(38)-CGD*U(18)
        G(59)=G(59)-CGS*U(18)
        G(60)=G(60)-CGD*U(18)
        G(77)=G(77)-CGS*U(18)
        G(78)=G(78)-CGD*U(18)

C
C --- Capacitive coupling result node nor-gate 2
C
        G(36)=G(36)+CGS*(U(36)-U(50))+CGD*(U(36)-U(51))
        G(50)=G(50)-CGS*U(36)
        G(51)=G(51)-CGD*U(36)

C
C --- Capacitive coupling result node andoi-gate 2 === s0
C
        G(49)=G(49)+COUT*U(49)

C
C --- Capacitive coupling result node andoi-gate 3
C
        G(67)=G(67)+CGS*(U(67)-U(117))+CGD*(U(67)-U(118))+
     *              CGS*(U(67)-U(136))+CGD*(U(67)-U(137))
        G(117)=G(117)-CGS*U(67)
        G(118)=G(118)-CGD*U(67)
        G(136)=G(136)-CGS*U(67)
        G(137)=G(137)-CGD*U(67)

C
C --- Capacitive coupling result node nor-gate 3
C
        G(85)=G(85)+CGS*(U(85)-U(99))+CGD*(U(85)-U(100))+
     *             CGS*(U(85)-U(149))+CGD*(U(85)-U(150))
        G(99)=G(99)-CGS*U(85)
        G(100)=G(100)-CGD*U(85)
        G(149)=G(149)-CGS*U(85)
        G(150)=G(150)-CGD*U(85)

C
C --- Capacitive coupling result node andoi-gate 4
C
        G(98)=G(98)+CGS*(U(98)-U(122))+CGD*(U(98)-U(123))+
     *              CGS*(U(98)-U(140))+CGD*(U(98)-U(141))+
     *              CGS*(U(98)-U(158))+CGD*(U(98)-U(159))+
     *              CGS*(U(98)-U(162))+CGD*(U(98)-U(163))
        G(122)=G(122)-CGS*U(98)
        G(123)=G(123)-CGD*U(98)
        G(140)=G(140)-CGS*U(98)
        G(141)=G(141)-CGD*U(98)
        G(158)=G(158)-CGS*U(98)
        G(159)=G(159)-CGD*U(98)

C
C --- Capacitive coupling result nand-gate
C
        G(116)=G(116)+CGS*(U(116)-U(131))+CGD*(U(116)-U(132))
        G(131)=G(131)-CGS*U(116)
        G(132)=G(132)-CGD*U(116)

C
C --- Capacitive coupling result node orani-gate === s1
C
        G(130)=G(130)+COUT*U(130)

C
C --- Capacitive coupling result andoi-gate 5 === Cinvers
C
        G(148)=G(148)+CBDBS(U(165)-U(148))*(U(148)-U(165))+COUT*U(148)

C
C --- Charge-function of three additional transistors
C
        G(162)=G(162)+CGS*(U(162)-U(98))
        G(163)=G(163)+CGD*(U(163)-U(98))
        G(164)=G(164)+CBDBS(U(164)-U(166))*(U(164)-U(166))
        G(165)=G(165)+CBDBS(U(165)-U(148))*(U(165)-U(148))
        G(166)=G(166)+CBDBS(U(164)-U(166))*(U(166)-U(164))+
     *                CBDBS(U(170)-U(166))*(U(166)-U(170))+
     *                CLOAD*U(166)
        G(167)=G(167)+CGS*(U(167)-U(18))
        G(168)=G(168)+CGD*(U(168)-U(18))
        G(169)=G(169)+CBDBS(U(169)-U(171))*(U(169)-U(171))
        G(170)=G(170)+CBDBS(U(170)-U(166))*(U(170)-U(166))
        G(171)=G(171)+CBDBS(U(169)-U(171))*(U(171)-U(169))+
     *                CBDBS(U(175)-U(171))*(U(171)-U(175))+
     *                CLOAD*U(171)
        G(172)=G(172)+CGS*U(172)
        G(173)=G(173)+CGD*U(173)
        G(174)=G(174)+CBDBS(U(174))*U(174)
        G(175)=G(175)+CBDBS(U(175)-U(171))*(U(175)-U(171))

        RETURN
        END

        SUBROUTINE DNOR(N,U,I,G)
C ---------------------------------------------------------------------------
C
C Charge-function of the NOR-gate:   NOT (U1 OR U2)
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   U    State-vector (only node potentials in reduced model)
C   I    Number of first node
C
C Output parameter:
C   G    Computed charge-function of theNOR-gate
C External reference:
C   CBDBS Function evaluating the voltage-dependent bulk-capacitance
C         in the MOS-model due to Shichman and Hodges
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I
        DIMENSION U(N),G(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

        EXTERNAL CBDBS

        G(I)=G(I)+CGS*(U(I)-U(I+4))
        G(I+1)=G(I+1)+CGD*(U(I+1)-U(I+4))
        G(I+2)=G(I+2)+CBDBS(U(I+2)-U(I+4))*(U(I+2)-U(I+4))
        G(I+3)=G(I+3)+CBDBS(U(I+3)-VDD)*U(I+3)
        G(I+4)=G(I+4)+CGS*(U(I+4)-U(I))+CGD*(U(I+4)-U(I+1))
     *               + CBDBS(U(I+2)-U(I+4))*(U(I+4)-U(I+2))
     *               + CBDBS(U(I+8)-U(I+4))*(U(I+4)-U(I+8))
     *               + CBDBS(U(I+12)-U(I+4))*(U(I+4)-U(I+12))
     *               + CLOAD*U(I+4)
        G(I+5)=G(I+5)+CGS*U(I+5)
        G(I+6)=G(I+6)+CGD*U(I+6)
        G(I+7)=G(I+7)+CBDBS(U(I+7))*U(I+7)
        G(I+8)=G(I+8)+CBDBS(U(I+8)-U(I+4))*(U(I+8)-U(I+4))
        G(I+9)=G(I+9)+CGS*U(I+9)
        G(I+10)=G(I+10)+CGD*U(I+10)
        G(I+11)=G(I+11)+CBDBS(U(I+11))*U(I+11)
        G(I+12)=G(I+12)+CBDBS(U(I+12)-U(I+4))*(U(I+12)-U(I+14))

        RETURN
        END

        SUBROUTINE DANDOI(N,U,I,G)
C ---------------------------------------------------------------------------
C
C Charge-function of the ANDOI-gate:  NOT ( U1 OR ( U2 AND U3) )
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   U    State-vector (only node potentials in reduced model)
C   I    Number of first node
C
C Output parameter:
C   G    Computed charge-function of the ANDOI-gate
C External reference:
C   CBDBS Function evaluating the voltage-dependent bulk-capacitance
C         in the MOS-model due to Shichman and Hodges
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I
        DIMENSION U(N),G(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

        EXTERNAL CBDBS

        G(I)=G(I)+CGS*(U(I)-U(I+4))
        G(I+1)=G(I+1)+CGD*(U(I+1)-U(I+4))
        G(I+2)=G(I+2)+CBDBS(U(I+2)-U(I+4))*(U(I+2)-U(I+4))
        G(I+3)=G(I+3)+CBDBS(U(I+3)-VDD)*U(I+3)
        G(I+4)=G(I+4)+CGS*(U(I+4)-U(I))+CGD*(U(I+4)-U(I+1))
     *               + CBDBS(U(I+2)-U(I+4))*(U(I+4)-U(I+2))
     *               + CBDBS(U(I+8)-U(I+4))*(U(I+4)-U(I+8))
     *               + CBDBS(U(I+12)-U(I+4))*(U(I+4)-U(I+12))
     *               + CLOAD*U(I+4)
        G(I+5)=G(I+5)+CGS*U(I+5)
        G(I+6)=G(I+6)+CGD*U(I+6)
        G(I+7)=G(I+7)+CBDBS(U(I+7))*U(I+7)
        G(I+8)=G(I+8)+CBDBS(U(I+8)-U(I+4))*(U(I+8)-U(I+4))
        G(I+9)=G(I+9)+CGS*U(I+9)
        G(I+10)=G(I+10)+CGD*U(I+10)
        G(I+11)=G(I+11)+CBDBS(U(I+11)-U(I+13))*(U(I+11)-U(I+13))
        G(I+12)=G(I+12)+CBDBS(U(I+12)-U(I+4))*(U(I+12)-U(I+4))
        G(I+13)=G(I+13)+CBDBS(U(I+11)-U(I+13))*(U(I+13)-U(I+11))
     *                 +CBDBS(U(I+17)-U(I+13))*(U(I+13)-U(I+17))
     *                 +CLOAD*U(I+13)
        G(I+14)=G(I+14)+CGS*U(I+14)
        G(I+15)=G(I+15)+CGD*U(I+15)
        G(I+16)=G(I+16)+CBDBS(U(I+16))*U(I+16)
        G(I+17)=G(I+17)+CBDBS(U(I+17)-U(I+13))*(U(I+17)-U(I+13))

        RETURN
        END

        SUBROUTINE DNAND(N,U,I,G)
C ---------------------------------------------------------------------------
C
C Charge-function of the NAND-gate:  NOT (U1 AND U2)
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   U    State-vector (only node potentials in reduced model)
C   I    Number of first node
C
C Output parameter:
C   G    Computed charge-function of the NAND-gate
C External reference:
C   CBDBS Function evaluating the voltage-dependent bulk-capacitance
C         in the MOS-model due to Shichman and Hodges
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I
        DIMENSION U(N),G(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

        EXTERNAL CBDBS

        G(I)=G(I)+CGS*(U(I)-U(I+4))
        G(I+1)=G(I+1)+CGD*(U(I+1)-U(I+4))
        G(I+2)=G(I+2)+CBDBS(U(I+2)-U(I+4))*(U(I+2)-U(I+4))
        G(I+3)=G(I+3)+CBDBS(U(I+3)-VDD)*U(I+3)
        G(I+4)=G(I+4)+CGS*(U(I+4)-U(I))+CGD*(U(I+4)-U(I+1))
     *               + CBDBS(U(I+2)-U(I+4))*(U(I+4)-U(I+2))
     *               + CBDBS(U(I+8)-U(I+4))*(U(I+4)-U(I+8))
     *               + CLOAD*U(I+4)
        G(I+5)=G(I+5)+CGS*U(I+5)
        G(I+6)=G(I+6)+CGD*U(I+6)
        G(I+7)=G(I+7)+CBDBS(U(I+7)-U(I+9))*(U(I+7)-U(I+9))
        G(I+8)=G(I+8)+CBDBS(U(I+8)-U(I+4))*(U(I+8)-U(I+4))
        G(I+9)=G(I+9)+CBDBS(U(I+7)-U(I+9))*(U(I+9)-U(I+7))
     *               +CBDBS(U(I+13)-U(I+9))*(U(I+9)-U(I+13))
     *               +CLOAD*U(I+9)
        G(I+10)=G(I+10)+CGS*U(I+10)
        G(I+11)=G(I+11)+CGD*U(I+11)
        G(I+12)=G(I+12)+CBDBS(U(I+12))*U(I+12)
        G(I+13)=G(I+13)+CBDBS(U(I+13)-U(I+9))*(U(I+13)-U(I+9))

        RETURN
        END

        SUBROUTINE DORANI(N,U,I,G)
C ---------------------------------------------------------------------------
C
C Charge-function of the ORANI-gate:  NOT ( U1 AND ( U2 OR U3) )
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   N    Dimension of the system
C   U    State-vector (only node potentials in reduced model)
C   I    Number of first node
C
C Output parameter:
C   G    Computed charge-function of the ORANI-gate
C External reference:
C   CBDBS Function evaluating the voltage-dependent bulk-capacitance
C         in the MOS-model due to Shichman and Hodges
C
C ---------------------------------------------------------------------------

        IMPLICIT double precision (A-H,O-Z)
        integer*4 N,I
        DIMENSION U(N),G(N)
        COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

        EXTERNAL CBDBS

        G(I)=G(I)+CGS*(U(I)-U(I+4))
        G(I+1)=G(I+1)+CGD*(U(I+1)-U(I+4))
        G(I+2)=G(I+2)+CBDBS(U(I+2)-U(I+4))*(U(I+2)-U(I+4))
        G(I+3)=G(I+3)+CBDBS(U(I+3)-VDD)*U(I+3)
        G(I+4)=G(I+4)+CGS*(U(I+4)-U(I))+CGD*(U(I+4)-U(I+1))
     *               + CBDBS(U(I+2)-U(I+4))*(U(I+4)-U(I+2))
     *               + CBDBS(U(I+8)-U(I+4))*(U(I+4)-U(I+8))
     *               + CLOAD*U(I+4)
        G(I+5)=G(I+5)+CGS*U(I+5)
        G(I+6)=G(I+6)+CGD*U(I+6)
        G(I+7)=G(I+7)+CBDBS(U(I+7)-U(I+9))*(U(I+7)-U(I+9))
        G(I+8)=G(I+8)+CBDBS(U(I+8)-U(I+4))*(U(I+8)-U(I+4))
        G(I+9)=G(I+9)+CBDBS(U(I+7)-U(I+9))*(U(I+9)-U(I+7)) +
     *                 CBDBS(U(I+13)-U(I+9))*(U(I+9)-U(I+13)) +
     *                 CBDBS(U(I+17)-U(I+9))*(U(I+9)-U(I+17))+
     *                 CLOAD*U(I+9)
        G(I+10)=G(I+10)+CGS*U(I+10)
        G(I+11)=G(I+11)+CGD*U(I+11)
        G(I+12)=G(I+12)+CBDBS(U(I+12))*U(I+12)
        G(I+13)=G(I+13)+CBDBS(U(I+13)-U(I+9))*(U(I+13)-U(I+9))
        G(I+14)=G(I+14)+CGS*U(I+14)
        G(I+15)=G(I+15)+CGD*U(I+15)
        G(I+16)=G(I+16)+CBDBS(U(I+16))*U(I+16)
        G(I+17)=G(I+17)+CBDBS(U(I+17)-U(I+9))*(U(I+17)-U(I+9))

        RETURN
        END

      double precision FUNCTION IDS (NED,VDS, VGS, VBS, VGD, VBD, ierr)
C ---------------------------------------------------------------------------
C
C Function evaluating the drain-current due to the model of
C Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   NED  Integer parameter for MOSFET-type
C   VDS  Voltage between drain and source
C   VGS  Voltage between gate and source
C   VGD  Voltage between gate and drain
C   VBD  Voltage between bulk and drain
C   I    Number of first node
C
C External reference:
C   GDSP, GDSM Drain function for VDS > 0 resp. VDS < 0
C
C ---------------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      integer*4 NED,ierr
      EXTERNAL GDSP, GDSM

      IF ( VDS .GT. 0.d0 ) THEN
         IDS = GDSP (NED,VDS, VGS, VBS,ierr)
      ELSE IF ( VDS .EQ. 0.d0) THEN
         IDS = 0.d0
      ELSE IF ( VDS .LT. 0.d0) THE N
         IDS = GDSM (NED,VDS, VGD, VBD,ierr)
      END IF

      if(ierr.eq.-1)return

      RETURN
      END

      double precision FUNCTION GDSP (NED,VDS, VGS, VBS,ierr)
      IMPLICIT double precision (A-H,O-Z)
      integer*4 NED,ierr

      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT
C

      IF(NED.EQ.0) THEN
C --- Depletion-type
        VT0=-2.43d0
        CGAMMA=.2d0
        PHI=1.28d0
        BETA=53.5D-6*CTIME*STIFF
      ELSE IF (NED.EQ.1) THEN
C --- Enhancement-type
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=4*43.7D-6*CTIME*STIFF
      ELSE IF (NED.EQ.2) THEN
C --- Two enhancement-type transistors in series
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=8*43.7D-6*CTIME*STIFF
      ELSE
C --- Three enhancement-type transistors in series
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=12*43.7D-6*CTIME*STIFF
      END IF

      if(phi-vbs.lt.0d0.or.phi.lt.0d0)then
         GDSP = 0.d0
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBS) - DSQRT(PHI) )

C      IF ( VGS-VTE .LE. 0.d0) THEN
C         GDSP = 0.d0
C      ELSE IF ( 0.d0 .LT. VGS-VTE .AND. VGS-VTE .LE. VDS ) THEN
C         GDSP = - BETA * (VGS - VTE)**2.d0 * (1.d0 + DELTA*VDS)
C      ELSE IF ( 0.d0 .LT. VDS .AND. VDS .LT. VGS-VTE ) THEN
C         GDSP = - BETA * VDS * (2.d0*(VGS - VTE) - VDS) *
C     *          (1.d0 + DELTA*VDS)
C      END IF
      IF ( VGS-VTE .LE. 0.d0) THEN
          GDSP = 0.d0
      ELSE IF (  VGS-VTE .LE. VDS ) THEN
          GDSP = - BETA * (VGS - VTE)**2.d0 * (1.d0 + DELTA*VDS)
      ELSE 
          GDSP = - BETA * VDS * (2.d0*(VGS - VTE) - VDS) *
     *          (1.d0 + DELTA*VDS)
      END IF

      RETURN
      END

      double precision FUNCTION GDSM (NED,VDS, VGD, VBD, ierr)
      IMPLICIT double precision (A-H,O-Z)
      integer*4 NED,ierr

      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

      IF(NED.EQ.0) THEN
C --- Depletion-type
        VT0=-2.43d0
        CGAMMA=.2d0
        PHI=1.28d0
        BETA=53.5D-6*CTIME*STIFF
      ELSE IF (NED.EQ.1) THEN
C --- Enhancement-type
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=4*43.7D-6*CTIME*STIFF
      ELSE IF (NED.EQ.2) THEN
C --- Two enhancement-type transistors in series
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=8*43.7D-6*CTIME*STIFF
      ELSE
C --- Three enhancement-type transistors in series
        VT0=.2d0
        CGAMMA=0.035d0
        PHI=1.01d0
        BETA=12*43.7D-6*CTIME*STIFF
      END IF

      if(phi-vbd.lt.0d0.or.phi.lt.0d0)then
          GDSM = 0.d0
         ierr=-1
         return
      end if

      VTE = VT0 + CGAMMA * ( DSQRT(PHI-VBD) - DSQRT(PHI) )

C      IF ( VGD-VTE .LE. 0.d0) THEN
C         GDSM = 0.d0
C      ELSE IF ( 0.d0 .LT. VGD-VTE .AND. VGD-VTE .LE. -VDS ) THEN
C         GDSM = BETA * (VGD - VTE)*(VGD - VTE) * (1.d0 - DELTA*VDS)
C      ELSE IF ( 0.d0 .LT. -VDS .AND. -VDS .LT. VGD-VTE ) THEN
C         GDSM = - BETA * VDS * (2.0d0*(VGD - VTE) + VDS) *
C     *          (1.d0 - DELTA*VDS)
C      END IF
      IF ( VGD-VTE .LE. 0.d0) THEN
         GDSM = 0.d0
      ELSE IF ( VGD-VTE .LE. -VDS ) THEN
         GDSM = BETA * (VGD - VTE)*(VGD - VTE) * (1.d0 - DELTA*VDS)
      ELSE 
        GDSM = - BETA * VDS * (2.0d0*(VGD - VTE) + VDS) *
     *          (1.d0 - DELTA*VDS)
      END IF

      RETURN
      END

      double precision FUNCTION IBS (VBS)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C source due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and source
C
C ---------------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

C

C
C     IBS = GBS (VBS)
C
      IF ( VBS .LE. 0.d0 ) THEN
         IBS = - CURIS * ( DEXP( VBS/VTH ) - 1.d0 )
C         IBS = - CURIS * ( DEXP(MAX(8., VBS/VTH )) - 1.d0 )
      ELSE
         IBS = 0.d0
      END IF

      RETURN
      END

      double precision FUNCTION IBD (VBD)
C ---------------------------------------------------------------------------
C
C Function evaluating the current of the pn-junction between bulk and
C drain  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   VBS  Voltage between bulk and drain
C
C ---------------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

C

C
C     IBD = GBD (VBD)
C
      IF ( VBD .LE. 0.d0 ) THEN
         IBD = - CURIS * ( DEXP( VBD/VTH ) - 1.d0 )
      ELSE
         IBD = 0.d0
      END IF
      RETURN
      END

      double precision FUNCTION CBDBS (V)
C ---------------------------------------------------------------------------
C
C Function evaluating the voltage-dependent capacitance between bulk and
C drain resp. source  due to the model of Shichman and Hodges
C
C ---------------------------------------------------------------------------
C
C The input parameters are:
C   V    Voltage between bulk and drain resp. source
C
C ---------------------------------------------------------------------------

      IMPLICIT double precision (A-H,O-Z)
      COMMON /CONST/ RGS, RGD, RBS, RBD, CGS, CGD, CBD, CBS,
     *               DELTA, CTIME, STIFF,
     *               CURIS, VTH, VDD, VBB, CLOAD, COUT

      PHIB=0.87d0

      IF ( V .LE. 0.d0 ) THEN
         CBDBS = CBD/DSQRT(1.d0-V/PHIB)
      ELSE
         CBDBS = CBD*(1.d0+V/(2.d0*PHIB))
      END IF

      RETURN
      END
c
