PROGRAM TEST_INSTALL
! Driver code that tests the installation of DELAUNAYSPARSES and
! DELAUNAYSPARSEP. To do so, a toy interpolation problem is
! computed and the results are compared to the known solution.

! Last Update: February, 2019
! Primary Author: Tyler Chang
USE DELSPARSE_MOD
USE OMP_LIB
IMPLICIT NONE

! Declare data.
INTEGER :: SIMPS(3,6), IERR(6)
REAL(KIND=R8) :: EPS
REAL(KIND=R8) :: INTERP_IN(1,20), INTERP_OUT(1,6), EXPECTED_OUT(1,6), &
  & PTS(2,20), PTS_TMP(2,20), Q(2,6), Q_TMP(2,6), WEIGHTS(3,6)

EPS = SQRT(EPSILON(0.0_R8))
PTS = TRANSPOSE( RESHAPE( (/  &
  0.10877683233208346_R8,     &
  0.65747571677546268_R8,     &
  0.74853271200744009_R8,     &
  0.25853058969031051_R8,     &
  0.38508322804628770_R8,     &
  0.19855613243388937_R8,     &
  0.88590610193360986_R8,     &
  0.73957680789581970_R8,     &
  0.46130107231752082_R8,     &
  0.61044888569019906_R8,     &
  0.88848755836796889_R8,     &
  0.56504950910258156_R8,     &
  0.63374920061262452_R8,     &
  0.47642100637444385_R8,     &
  0.89167673297718886_R8,     &
  0.85575976312324076_R8,     &
  0.36741400280848768_R8,     &
  0.22540743314109113_R8,     &
  0.57887702455276135_R8,     &
  0.33794226559725304_R8,     &
  0.76211800269757757_R8,     &
  0.082963515866522064_R8,    &
  0.016220459783666152_R8,    &
  0.17155847087049503_R8,     &
  0.12930597950925682_R8,     &
  0.91552991190955113_R8,     &
  0.30469899967300274_R8,     &
  0.064234640774060825_R8,    &
  0.67129213095523377_R8,     &
  0.56860397761470494_R8,     &
  0.10547481357911370_R8,     &
  0.59408216854500884_R8,     &
  0.90989152079869851_R8,     &
  0.91232248805035077_R8,     &
  0.13873375923421827_R8,     &
  0.68652421762380056_R8,     &
  0.53775708104383380_R8,     &
  0.63512621583969442_R8,     &
  0.98798019619988187_R8,     &
  0.87480704030477330_R8  /), &
  (/ 20, 2 /) ) )
Q = TRANSPOSE( RESHAPE( (/ &
  0.500000000000000000_R8,    &
  0.250000000000000000_R8,    &
  0.250000000000000000_R8,    &
  0.750000000000000000_R8,    &
  0.750000000000000000_R8,    &
  0.100000000000000000_R8,    &
  0.500000000000000000_R8,    &
  0.250000000000000000_R8,    &
  0.750000000000000000_R8,    &
  0.250000000000000000_R8,    &
  0.750000000000000000_R8,    &
  0.500000000000000000_R8 /), &
  (/6, 2/) ) )
INTERP_IN = RESHAPE( (/ &
  0.87089483502966103_R8,     &
  0.74043923264198475_R8,     &
  0.76475317179110625_R8,     &
  0.43008906056080554_R8,     &
  0.51438920755554451_R8,     &
   1.1140860443434404_R8,     &
   1.1906051016066126_R8,     &
  0.80381144866988052_R8,     &
   1.1325932032727546_R8,     &
   1.1790528633049040_R8,     &
  0.99396237194708259_R8,     &
   1.1591316776475904_R8,     &
   1.5436407214113230_R8,     &
   1.3887434944247947_R8,     &
   1.0304104922114070_R8,     &
   1.5422839807470412_R8,     &
  0.90517108385232148_R8,     &
  0.86053364898078555_R8,     &
   1.5668572207526432_R8,     &
   1.2127493059020265_R8 /),  &
  (/ 1, 20 /) )
EXPECTED_OUT = RESHAPE( (/ &
  1.00000000000000000_R8,    &
  0.50000000000000000_R8,    &
  1.00000000000000000_R8,    &
  1.00000000000000000_R8,    &
  1.50000000000000000_R8,    &
  0.68862615900613189_R8 /), &
  (/ 1, 6/) )

! Test DELAUNAYSPARSES.
PTS_TMP = PTS; Q_TMP = Q
CALL DELAUNAYSPARSES(2, 20, PTS_TMP, 6, Q_TMP, SIMPS, WEIGHTS, IERR, &
  & INTERP_IN=INTERP_IN, INTERP_OUT=INTERP_OUT)
IF(ANY(ABS(INTERP_OUT - EXPECTED_OUT) > EPS)) THEN
   WRITE(*,*) "DELAUNAYSPARSES produced an incorrect result. ", &
     & " The installation is not correct."
   STOP
END IF

! Test DELAUNAYSPARSEP, PMODE=1.
PTS_TMP = PTS; Q_TMP = Q
CALL OMP_SET_NUM_THREADS(4)
CALL DELAUNAYSPARSEP(2, 20, PTS_TMP, 6, Q_TMP, SIMPS, WEIGHTS, IERR, &
  & INTERP_IN=INTERP_IN, INTERP_OUT=INTERP_OUT, PMODE=1)
IF(ANY(ABS(INTERP_OUT - EXPECTED_OUT) > EPS)) THEN
   WRITE(*,*) "DELAUNAYSPARSEP produced an incorrect result. ", &
     & " The installation is not correct."
   STOP
END IF

! Test DELAUNAYSPARSEP, PMODE=2.
PTS_TMP = PTS; Q_TMP = Q
CALL OMP_SET_NUM_THREADS(4)
CALL DELAUNAYSPARSEP(2, 20, PTS_TMP, 6, Q_TMP, SIMPS, WEIGHTS, IERR, &
  & INTERP_IN=INTERP_IN, INTERP_OUT=INTERP_OUT, PMODE=2)
IF(ANY(ABS(INTERP_OUT - EXPECTED_OUT) > EPS)) THEN
   WRITE(*,*) "DELAUNAYSPARSEP produced an incorrect result. ", &
     & " The installation is not correct."
   STOP
END IF

! Test DELAUNAYSPARSEP, PMODE=3.
CALL OMP_SET_NESTED(.TRUE.)
CALL OMP_SET_NUM_THREADS(2)
PTS_TMP = PTS; Q_TMP = Q
CALL DELAUNAYSPARSEP(2, 20, PTS_TMP, 6, Q_TMP, SIMPS, WEIGHTS, IERR, &
  & INTERP_IN=INTERP_IN, INTERP_OUT=INTERP_OUT, PMODE=3)
IF(ANY(ABS(INTERP_OUT - EXPECTED_OUT) > EPS)) THEN
   WRITE(*,*) "DELAUNAYSPARSEP produced an incorrect result. ", &
     & " The installation is not correct."
   STOP
END IF

! If all the tests passed, then the installation is correct.
WRITE(*,*) "The installation of DELAUNAYSPARSE appears correct."

END PROGRAM TEST_INSTALL