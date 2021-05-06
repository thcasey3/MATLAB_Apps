#DESC	1.2 * DESCRIPTOR INFORMATION ***********************
*
*	Dataset Type and Format:
*
DSRC	EXP
BSEQ	BIG
IKKF	REAL
XTYP	IDX
YTYP	IGD
ZTYP	NODATA
*
*	Item Formats:
*
IRFMT	D
YFMT	D
*
*	Data Ranges and Resolutions:
*
XPTS	4500
XMIN	3416.800000
XWID	149.966667
YPTS	100
YMIN	0.000000
YWID	461.400000
*
*	Documentational Text:
*
TITL	'Example2D'
IRNAM	'Intensity'
XNAM	'Field'
YNAM	'Time'
IRUNI	''
XUNI	'G'
YUNI	's'
*
************************************************************
*
#SPL	1.2 * STANDARD PARAMETER LAYER
*
OPER    xuser
DATE    05/31/19
TIME    14:05:29
CMNT    
SAMP    
SFOR    
STAG    C
EXPT    CW
OXS1    IADC
AXS1    B0VL
AXS2    ETIM
AXS3    
A2RS    100
A1CT    0.34918
A1SW    0.015
MWFQ    9.808098e+09
MWPW    0.001002
AVGS    1
SPTP    0.00067
RCAG    30
RCHM    1
B0MA    0.0001
B0MF    100000
RCPH    0.0
RCOF    0.0
A1RS    4500
RCTC    0.00032
*
************************************************************
*
#DSL	1.0 * DEVICE SPECIFIC LAYER
*

.DVC     acqStart, 1.0


.DVC     delay, 1.0

Delay              100 ms
NbPoints           100

.DVC     ezcw, 1.0

MetalIonsDef       False
OrgRadicalsDef     False

.DVC     fieldCtrl, 1.0

AllegroMode        True
CenterField        3491.80 G
Delay              0.0 s
FieldFlyback       On
FieldWait          Wait LED off
GFactor            2.000000
MeasuringHall      False
SetToSampleG       False
StaticFieldMon     3416.800 G
SweepDirection     Up
SweepWidth         150.0 G
WidthTM            200.0 G

.DVC     fieldSweep, 1.0


.DVC     freqCounter, 1.0

FrequencyMon       9.808098 GHz
QMonitBridge       On

.DVC     mwBridge, 1.0

AcqFineTuning      Never
AcqScanFTuning     Off
AcqSliceFTuning    Off
BridgeCalib        83.0
Power              1.002 mW
PowerAtten         23.0 dB
QValue             7900

.DVC     recorder, 1.0

BaselineCorr       On
NbScansAcc         1
NbScansDone        1
NbScansToDo        1
ReplaceMode        Off
SmoothMode         Manual
SmoothPoints       1

.DVC     scanEnd, 1.0


.DVC     signalChannel, 1.0

AFCTrap            True
AllegroMode        True
AllowShortCt       False
Calibrated         True
ConvFact           9.271e-09
ConvTime           0.67 ms
DModDetectSCT      First
DualDetect         OFF
EliDelay           1.0 us
EnableImag         Disable
ExtLockIn          False
ExtTrigger         False
Gain               30 dB
Harmonic           1
HighPass           True
Integrator         False
IsCalibExp         False
ModAmp             1.000 G
ModFreq            100.00 kHz
ModPhase           0.0
ModResol           30.000
Offset             0.0 %
PolyCof            {2;3,9;0} 0,0.99652,0,0,0.00737177,0,0,-0.00559614,0,0,-2.88221e-05,0,0,1.00404e-05,0,0,3.43695e-08,0,0,-5.0404e-09,0,0,-1.4783e-11,0,0,-1.29132e-12,0
QuadMode           False
ResCenter          62.5 mm
ResLength          40 mm
Resolution         4500
Resonator          1
SctNorm            True
SctRevision        Allegro
SpuExtension       False
SweepTime          3.02 s
TimeConst          0.32 ms
TimeExp            False
TuneCaps           42

*
************************************************************
