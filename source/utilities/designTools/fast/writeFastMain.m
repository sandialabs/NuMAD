function writeFastMain(fast,output_file)
%WRITEFASTMAIN  Write a FAST primary input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeFastMain(fast,'file_name') 
%      Created for FAST v7.00
%
%      fast = FAST data structure; view default structure by utilizing 
%           readFastMain()
%      output_file = name of file to be created.

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],'--------------------------------------------------------------------------------')
wip(fid,[],'------- FAST INPUT FILE --------------------------------------------------------')
wip(fid,[],fast.title{1})
wip(fid,[],fast.title{2})
wip(fid,[],'---------------------- SIMULATION CONTROL --------------------------------------')
wip(fid,fast.SimCtrl.Echo,         'Echo        - Echo input data to "echo.out" (flag)')
wip(fid,fast.SimCtrl.ADAMSPrep,    'ADAMSPrep   - ADAMS preprocessor mode {1: Run FAST, 2: use FAST as a preprocessor to create an ADAMS model, 3: do both} (switch)')
wip(fid,fast.SimCtrl.AnalMode,     'AnalMode    - Analysis mode {1: Run a time-marching simulation, 2: create a periodic linearized model} (switch)')
wip(fid,fast.SimCtrl.NumBl,        'NumBl       - Number of blades (-)')
wip(fid,fast.SimCtrl.TMax,         'TMax        - Total run time (s)')
wip(fid,fast.SimCtrl.DT,           'DT          - Integration time step (s)')
wip(fid,[],'---------------------- TURBINE CONTROL -----------------------------------------')
wip(fid,fast.TurbCtrl.YCMode,      'YCMode      - Yaw control mode {0: none, 1: user-defined from routine UserYawCont, 2: user-defined from Simulink} (switch)')
wip(fid,fast.TurbCtrl.TYCOn,       'TYCOn       - Time to enable active yaw control (s) [unused when YCMode=0]')
wip(fid,fast.TurbCtrl.PCMode,      'PCMode      - Pitch control mode {0: none, 1: user-defined from routine PitchCntrl, 2: user-defined from Simulink} (switch)')
wip(fid,fast.TurbCtrl.TPCOn,       'TPCOn       - Time to enable active pitch control (s) [unused when PCMode=0]')
wip(fid,fast.TurbCtrl.VSContrl,    'VSContrl    - Variable-speed control mode {0: none, 1: simple VS, 2: user-defined from routine UserVSCont, 3: user-defined from Simulink} (switch)')
wip(fid,fast.TurbCtrl.VS_RtGnSp,   'VS_RtGnSp   - Rated generator speed for simple variable-speed generator control (HSS side) (rpm) [used only when VSContrl=1]')
wip(fid,fast.TurbCtrl.VS_RtTq,     'VS_RtTq     - Rated generator torque/constant generator torque in Region 3 for simple variable-speed generator control (HSS side) (N-m) [used only when VSContrl=1]')
wip(fid,fast.TurbCtrl.VS_Rgn2K,    'VS_Rgn2K    - Generator torque constant in Region 2 for simple variable-speed generator control (HSS side) (N-m/rpm^2) [used only when VSContrl=1]')
wip(fid,fast.TurbCtrl.VS_SlPc,     'VS_SlPc     - Rated generator slip percentage in Region 2 1/2 for simple variable-speed generator control (%%) [used only when VSContrl=1]')
wip(fid,fast.TurbCtrl.GenModel,    'GenModel    - Generator model {1: simple, 2: Thevenin, 3: user-defined from routine UserGen} (switch) [used only when VSContrl=0]')
wip(fid,fast.TurbCtrl.GenTiStr,    'GenTiStr    - Method to start the generator {T: timed using TimGenOn, F: generator speed using SpdGenOn} (flag)')
wip(fid,fast.TurbCtrl.GenTiStp,    'GenTiStp    - Method to stop the generator {T: timed using TimGenOf, F: when generator power = 0} (flag)')
wip(fid,fast.TurbCtrl.SpdGenOn,    'SpdGenOn    - Generator speed to turn on the generator for a startup (HSS speed) (rpm) [used only when GenTiStr=False]')
wip(fid,fast.TurbCtrl.TimGenOn,    'TimGenOn    - Time to turn on the generator for a startup (s) [used only when GenTiStr=True]')
wip(fid,fast.TurbCtrl.TimGenOf,    'TimGenOf    - Time to turn off the generator (s) [used only when GenTiStp=True]')
wip(fid,fast.TurbCtrl.HSSBrMode,   'HSSBrMode   - HSS brake model {1: simple, 2: user-defined from routine UserHSSBr} (switch)')
wip(fid,fast.TurbCtrl.THSSBrDp,    'THSSBrDp    - Time to initiate deployment of the HSS brake (s)')
wip(fid,fast.TurbCtrl.TiDynBrk,    'TiDynBrk    - Time to initiate deployment of the dynamic generator brake [CURRENTLY IGNORED] (s)')
wip(fid,fast.TurbCtrl.TTpBrDp(1),  'TTpBrDp(1)  - Time to initiate deployment of tip brake 1 (s)')
wip(fid,fast.TurbCtrl.TTpBrDp(2),  'TTpBrDp(2)  - Time to initiate deployment of tip brake 2 (s)')
wip(fid,fast.TurbCtrl.TTpBrDp(3),  'TTpBrDp(3)  - Time to initiate deployment of tip brake 3 (s) [unused for 2 blades]')
wip(fid,fast.TurbCtrl.TBDepISp(1), 'TBDepISp(1) - Deployment-initiation speed for the tip brake on blade 1 (rpm)')
wip(fid,fast.TurbCtrl.TBDepISp(2), 'TBDepISp(2) - Deployment-initiation speed for the tip brake on blade 2 (rpm)')
wip(fid,fast.TurbCtrl.TBDepISp(3), 'TBDepISp(3) - Deployment-initiation speed for the tip brake on blade 3 (rpm) [unused for 2 blades]')
wip(fid,fast.TurbCtrl.TYawManS,    'TYawManS    - Time to start override yaw maneuver and end standard yaw control (s)')
wip(fid,fast.TurbCtrl.TYawManE,    'TYawManE    - Time at which override yaw maneuver reaches final yaw angle (s)')
wip(fid,fast.TurbCtrl.NacYawF,     'NacYawF     - Final yaw angle for yaw maneuvers (degrees)')
wip(fid,fast.TurbCtrl.TPitManS(1), 'TPitManS(1) - Time to start override pitch maneuver for blade 1 and end standard pitch control (s)')
wip(fid,fast.TurbCtrl.TPitManS(2), 'TPitManS(2) - Time to start override pitch maneuver for blade 2 and end standard pitch control (s)')
wip(fid,fast.TurbCtrl.TPitManS(3), 'TPitManS(3) - Time to start override pitch maneuver for blade 3 and end standard pitch control (s) [unused for 2 blades]')
wip(fid,fast.TurbCtrl.TPitManE(1), 'TPitManE(1) - Time at which override pitch maneuver for blade 1 reaches final pitch (s)')
wip(fid,fast.TurbCtrl.TPitManE(2), 'TPitManE(2) - Time at which override pitch maneuver for blade 2 reaches final pitch (s)')
wip(fid,fast.TurbCtrl.TPitManE(3), 'TPitManE(3) - Time at which override pitch maneuver for blade 3 reaches final pitch (s) [unused for 2 blades]')
wip(fid,fast.TurbCtrl.BlPitch(1),  'BlPitch(1)  - Blade 1 initial pitch (degrees)')
wip(fid,fast.TurbCtrl.BlPitch(2),  'BlPitch(2)  - Blade 2 initial pitch (degrees)')
wip(fid,fast.TurbCtrl.BlPitch(3),  'BlPitch(3)  - Blade 3 initial pitch (degrees) [unused for 2 blades]')
wip(fid,fast.TurbCtrl.BlPitchF(1), 'BlPitchF(1) - Blade 1 final pitch for pitch maneuvers (degrees)')
wip(fid,fast.TurbCtrl.BlPitchF(2), 'BlPitchF(2) - Blade 2 final pitch for pitch maneuvers (degrees)')
wip(fid,fast.TurbCtrl.BlPitchF(3), 'BlPitchF(3) - Blade 3 final pitch for pitch maneuvers (degrees) [unused for 2 blades]')
wip(fid,[],'---------------------- ENVIRONMENTAL CONDITIONS --------------------------------')
wip(fid,fast.Env.Gravity,          'Gravity     - Gravitational acceleration (m/s^2)')
wip(fid,[],'---------------------- FEATURE FLAGS -------------------------------------------')
wip(fid,fast.Flags.FlapDOF1,       'FlapDOF1    - First flapwise blade mode DOF (flag)')
wip(fid,fast.Flags.FlapDOF2,       'FlapDOF2    - Second flapwise blade mode DOF (flag)')
wip(fid,fast.Flags.EdgeDOF,        'EdgeDOF     - First edgewise blade mode DOF (flag)')
wip(fid,fast.Flags.TeetDOF,        'TeetDOF     - Rotor-teeter DOF (flag) [unused for 3 blades]')
wip(fid,fast.Flags.DrTrDOF,        'DrTrDOF     - Drivetrain rotational-flexibility DOF (flag)')
wip(fid,fast.Flags.GenDOF,         'GenDOF      - Generator DOF (flag)')
wip(fid,fast.Flags.YawDOF,         'YawDOF      - Yaw DOF (flag)')
wip(fid,fast.Flags.TwFADOF1,       'TwFADOF1    - First fore-aft tower bending-mode DOF (flag)')
wip(fid,fast.Flags.TwFADOF2,       'TwFADOF2    - Second fore-aft tower bending-mode DOF (flag)')
wip(fid,fast.Flags.TwSSDOF1,       'TwSSDOF1    - First side-to-side tower bending-mode DOF (flag)')
wip(fid,fast.Flags.TwSSDOF2,       'TwSSDOF2    - Second side-to-side tower bending-mode DOF (flag)')
wip(fid,fast.Flags.CompAero,       'CompAero    - Compute aerodynamic forces (flag)')
wip(fid,fast.Flags.CompNoise,      'CompNoise   - Compute aerodynamic noise (flag)')
wip(fid,[],'---------------------- INITIAL CONDITIONS --------------------------------------')
wip(fid,fast.Init.OoPDefl,         'OoPDefl     - Initial out-of-plane blade-tip displacement, (meters)')
wip(fid,fast.Init.IPDefl,          'IPDefl      - Initial in-plane blade-tip deflection, (meters)')
wip(fid,fast.Init.TeetDefl,        'TeetDefl    - Initial or fixed teeter angle (degrees) [unused for 3 blades]')
wip(fid,fast.Init.Azimuth,         'Azimuth     - Initial azimuth angle for blade 1 (degrees)')
wip(fid,fast.Init.RotSpeed,        'RotSpeed    - Initial or fixed rotor speed (rpm)')
wip(fid,fast.Init.NacYaw,          'NacYaw      - Initial or fixed nacelle-yaw angle (degrees)')
wip(fid,fast.Init.TTDspFA,         'TTDspFA     - Initial fore-aft tower-top displacement (meters)')
wip(fid,fast.Init.TTDspSS,         'TTDspSS     - Initial side-to-side tower-top displacement (meters)')
wip(fid,[],'---------------------- TURBINE CONFIGURATION -----------------------------------')
wip(fid,fast.TurbConf.TipRad,      'TipRad      - The distance from the rotor apex to the blade tip (meters)')
wip(fid,fast.TurbConf.HubRad,      'HubRad      - The distance from the rotor apex to the blade root (meters)')
wip(fid,fast.TurbConf.PSpnElN,     'PSpnElN     - Number of the innermost blade element which is still part of the pitchable portion of the blade for partial-span pitch control [1 to BldNodes] [CURRENTLY IGNORED] (-)')
wip(fid,fast.TurbConf.UndSling,    'UndSling    - Undersling length [distance from teeter pin to the rotor apex] (meters) [unused for 3 blades]')
wip(fid,fast.TurbConf.HubCM,       'HubCM       - Distance from rotor apex to hub mass [positive downwind] (meters)')
wip(fid,fast.TurbConf.OverHang,    'OverHang    - Distance from yaw axis to rotor apex [3 blades] or teeter pin [2 blades] (meters)')
wip(fid,fast.TurbConf.NacCMxn,     'NacCMxn     - Downwind distance from the tower-top to the nacelle CM (meters)')
wip(fid,fast.TurbConf.NacCMyn,     'NacCMyn     - Lateral  distance from the tower-top to the nacelle CM (meters)')
wip(fid,fast.TurbConf.NacCMzn,     'NacCMzn     - Vertical distance from the tower-top to the nacelle CM (meters)')
wip(fid,fast.TurbConf.TowerHt,     'TowerHt     - Height of tower above ground level [onshore] or MSL [offshore] (meters)')
wip(fid,fast.TurbConf.Twr2Shft,    'Twr2Shft    - Vertical distance from the tower-top to the rotor shaft (meters)')
wip(fid,fast.TurbConf.TwrRBHt,     'TwrRBHt     - Tower rigid base height (meters)')
wip(fid,fast.TurbConf.ShftTilt,    'ShftTilt    - Rotor shaft tilt angle (degrees)')
wip(fid,fast.TurbConf.Delta3,      'Delta3      - Delta-3 angle for teetering rotors (degrees) [unused for 3 blades]')
wip(fid,fast.TurbConf.PreCone(1),  'PreCone(1)  - Blade 1 cone angle (degrees)')
wip(fid,fast.TurbConf.PreCone(2),  'PreCone(2)  - Blade 2 cone angle (degrees)')
wip(fid,fast.TurbConf.PreCone(3),  'PreCone(3)  - Blade 3 cone angle (degrees) [unused for 2 blades]')
wip(fid,fast.TurbConf.AzimB1Up,    'AzimB1Up    - Azimuth value to use for I/O when blade 1 points up (degrees)')
wip(fid,[],'---------------------- MASS AND INERTIA ----------------------------------------')
wip(fid,fast.MassProp.YawBrMass,   'YawBrMass   - Yaw bearing mass (kg)')
wip(fid,fast.MassProp.NacMass,     'NacMass     - Nacelle mass (kg)')
wip(fid,fast.MassProp.HubMass,     'HubMass     - Hub mass (kg)')
wip(fid,fast.MassProp.TipMass(1),  'TipMass(1)  - Tip-brake mass, blade 1 (kg)')
wip(fid,fast.MassProp.TipMass(2),  'TipMass(2)  - Tip-brake mass, blade 2 (kg)')
wip(fid,fast.MassProp.TipMass(3),  'TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]')
wip(fid,fast.MassProp.NacYIner,    'NacYIner    - Nacelle inertia about yaw axis (kg m^2)')
wip(fid,fast.MassProp.GenIner,     'GenIner     - Generator inertia about HSS (kg m^2)')
wip(fid,fast.MassProp.HubIner,     'HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)')
wip(fid,[],'---------------------- DRIVETRAIN ----------------------------------------------')
wip(fid,fast.DrvTrn.GBoxEff,       'GBoxEff     - Gearbox efficiency (%%)')
wip(fid,fast.DrvTrn.GenEff,        'GenEff      - Generator efficiency [ignored by the Thevenin and user-defined generator models] (%%)')
wip(fid,fast.DrvTrn.GBRatio,       'GBRatio     - Gearbox ratio (-)')
wip(fid,fast.DrvTrn.GBRevers,      'GBRevers    - Gearbox reversal {T: if rotor and generator rotate in opposite directions} (flag)')
wip(fid,fast.DrvTrn.HSSBrTqF,      'HSSBrTqF    - Fully deployed HSS-brake torque (N-m)')
wip(fid,fast.DrvTrn.HSSBrDT,       'HSSBrDT     - Time for HSS-brake to reach full deployment once initiated (sec) [used only when HSSBrMode=1]')
wip(fid,fast.DrvTrn.DynBrkFi,      'DynBrkFi    - File containing a mech-gen-torque vs HSS-speed curve for a dynamic brake [CURRENTLY IGNORED] (quoted string)')
wip(fid,fast.DrvTrn.DTTorSpr,      'DTTorSpr    - Drivetrain torsional spring (N-m/rad)')
wip(fid,fast.DrvTrn.DTTorDmp,      'DTTorDmp    - Drivetrain torsional damper (N-m/s)')
wip(fid,[],'---------------------- SIMPLE INDUCTION GENERATOR ------------------------------')
wip(fid,fast.SIG.SIG_SlPc,         'SIG_SlPc    - Rated generator slip percentage (%%) [used only when VSContrl=0 and GenModel=1]')
wip(fid,fast.SIG.SIG_SySp,         'SIG_SySp    - Synchronous (zero-torque) generator speed (rpm) [used only when VSContrl=0 and GenModel=1]')
wip(fid,fast.SIG.SIG_RtTq,         'SIG_RtTq    - Rated torque (N-m) [used only when VSContrl=0 and GenModel=1]')
wip(fid,fast.SIG.SIG_PORt,         'SIG_PORt    - Pull-out ratio (Tpullout/Trated) (-) [used only when VSContrl=0 and GenModel=1]')
wip(fid,[],'---------------------- THEVENIN-EQUIVALENT INDUCTION GENERATOR -----------------')
wip(fid,fast.TEC.TEC_Freq,         'TEC_Freq    - Line frequency [50 or 60] (Hz) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_NPol,         'TEC_NPol    - Number of poles [even integer > 0] (-) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_SRes,         'TEC_SRes    - Stator resistance (ohms) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_RRes,         'TEC_RRes    - Rotor resistance (ohms) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_VLL,          'TEC_VLL     - Line-to-line RMS voltage (volts) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_SLR,          'TEC_SLR     - Stator leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_RLR,          'TEC_RLR     - Rotor leakage reactance (ohms) [used only when VSContrl=0 and GenModel=2]')
wip(fid,fast.TEC.TEC_MR,           'TEC_MR      - Magnetizing reactance (ohms) [used only when VSContrl=0 and GenModel=2]')
wip(fid,[],'---------------------- PLATFORM MODEL ------------------------------------------')
wip(fid,fast.Ptfm.PtfmModel,       'PtfmModel   - Platform model {0: none, 1: onshore, 2: fixed bottom offshore, 3: floating offshore} (switch)')
wip(fid,fast.Ptfm.PtfmFile,        'PtfmFile    - Name of file containing platform properties (quoted string) [unused when PtfmModel=0]')
wip(fid,[],'---------------------- TOWER ---------------------------------------------------')
wip(fid,fast.Twr.TwrNodes,         'TwrNodes    - Number of tower nodes used for analysis (-)')
wip(fid,fast.Twr.TwrFile,          'TwrFile - Name of file containing tower properties (quoted string)')
wip(fid,[],'---------------------- NACELLE-YAW ---------------------------------------------')
wip(fid,fast.Yaw.YawSpr,           'YawSpr      - Nacelle-yaw spring constant (N-m/rad)')
wip(fid,fast.Yaw.YawDamp,          'YawDamp     - Nacelle-yaw damping constant (N-m/rad/s)')
wip(fid,fast.Yaw.YawNeut,          'YawNeut     - Neutral yaw position--yaw spring force is zero at this yaw (degrees)')
wip(fid,[],'---------------------- FURLING -------------------------------------------------')
wip(fid,fast.Furl.Furling,         'Furling     - Read in additional model properties for furling turbine (flag)')
wip(fid,fast.Furl.FurlFile,        'FurlFile    - Name of file containing furling properties (quoted string) [unused when Furling=False]')
wip(fid,[],'---------------------- ROTOR-TEETER --------------------------------------------')
wip(fid,fast.Teet.TeetMod,         'TeetMod     - Rotor-teeter spring/damper model {0: none, 1: standard, 2: user-defined from routine UserTeet} (switch) [unused for 3 blades]')
wip(fid,fast.Teet.TeetDmpP,        'TeetDmpP    - Rotor-teeter damper position (degrees) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetDmp,         'TeetDmp     - Rotor-teeter damping constant (N-m/rad/s) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetCDmp,        'TeetCDmp    - Rotor-teeter rate-independent Coulomb-damping moment (N-m) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetSStP,        'TeetSStP    - Rotor-teeter soft-stop position (degrees) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetHStP,        'TeetHStP    - Rotor-teeter hard-stop position (degrees) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetSSSp,        'TeetSSSp    - Rotor-teeter soft-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]')
wip(fid,fast.Teet.TeetHSSp,        'TeetHSSp    - Rotor-teeter hard-stop linear-spring constant (N-m/rad) [used only for 2 blades and when TeetMod=1]')
wip(fid,[],'---------------------- TIP-BRAKE -----------------------------------------------')
wip(fid,fast.TpBr.TBDrConN,        'TBDrConN    - Tip-brake drag constant during normal operation, Cd*Area (m^2)')
wip(fid,fast.TpBr.TBDrConD,        'TBDrConD    - Tip-brake drag constant during fully-deployed operation, Cd*Area (m^2)')
wip(fid,fast.TpBr.TpBrDT,          'TpBrDT      - Time for tip-brake to reach full deployment once released (sec)')
wip(fid,[],'---------------------- BLADE ---------------------------------------------------')
wip(fid,fast.BldFile{1},           'BldFile(1) - Name of file containing properties for blade 1 (quoted string)')
wip(fid,fast.BldFile{2},           'BldFile(2) - Name of file containing properties for blade 2 (quoted string)')
wip(fid,fast.BldFile{3},           'BldFile(3) - Name of file containing properties for blade 3 (quoted string) [unused for 2 blades]')
wip(fid,[],'---------------------- AERODYN -------------------------------------------------')
wip(fid,fast.ADFile,               'ADFile     - Name of file containing AeroDyn input parameters (quoted string)')
wip(fid,[],'---------------------- NOISE ---------------------------------------------------')
wip(fid,fast.NoiseFile,            'NoiseFile   - Name of file containing aerodynamic noise input parameters (quoted string) [used only when CompNoise=True]')
wip(fid,[],'---------------------- ADAMS ---------------------------------------------------')
wip(fid,fast.ADAMSFile,            'ADAMSFile  - Name of file containing ADAMS-specific input parameters (quoted string) [unused when ADAMSPrep=1]')
wip(fid,[],'---------------------- LINEARIZATION CONTROL -----------------------------------')
wip(fid,fast.LinFile,              'LinFile    - Name of file containing FAST linearazation parameters (quoted string) [unused when AnalMode=1]')
wip(fid,[],'---------------------- OUTPUT --------------------------------------------------')
wip(fid,fast.Out.SumPrint,         'SumPrint    - Print summary data to "<RootName>.fsm" (flag)')
if isfield(fast.Out,'OutFileFmt')
  wip(fid,fast.Out.OutFileFmt,       'OutFileFmt  - Format for tabular (time-marching) output file(s) (1: text file [<RootName>.out], 2: binary file [<RootName>.outb], 3: both) (switch)')
end
wip(fid,fast.Out.TabDelim,         'TabDelim    - Generate a tab-delimited tabular output file. (flag)')
wip(fid,fast.Out.OutFmt,           'OutFmt      - Format used for tabular output except time.  Resulting field should be 10 characters. (quoted string)  [not checked for validity!]')
wip(fid,fast.Out.TStart,           'TStart      - Time to begin tabular output (s)')
wip(fid,fast.Out.DecFact,          'DecFact     - Decimation factor for tabular output {1: output every time step} (-)')
wip(fid,fast.Out.SttsTime,         'SttsTime    - Amount of time between screen status messages (sec)')
wip(fid,fast.Out.NcIMUxn,          'NcIMUxn     - Downwind distance from the tower-top to the nacelle IMU (meters)')
wip(fid,fast.Out.NcIMUyn,          'NcIMUyn     - Lateral  distance from the tower-top to the nacelle IMU (meters)')
wip(fid,fast.Out.NcIMUzn,          'NcIMUzn     - Vertical distance from the tower-top to the nacelle IMU (meters)')
wip(fid,fast.Out.ShftGagL,         'ShftGagL    - Distance from rotor apex [3 blades] or teeter pin [2 blades] to shaft strain gages [positive for upwind rotors] (meters)')
wip(fid,fast.Out.NTwGages,         'NTwGages    - Number of tower nodes that have strain gages for output [0 to 5] (-)')
wipcsv(fid,fast.Out.TwrGagNd,      'TwrGagNd    - List of tower nodes that have strain gages [1 to TwrNodes] (-) [unused if NTwGages=0]')
wip(fid,fast.Out.NBlGages,         'NBlGages    - Number of blade nodes that have strain gages for output [0 to 5] (-)')
wipcsv(fid,fast.Out.BldGagNd,      'BldGagNd    - List of blade nodes that have strain gages [1 to BldNodes] (-) [unused if NBlGages=0]')
wip(fid,'',                        'OutList     - The next line(s) contains a list of output parameters.  See OutList.txt for a listing of available output channels, (-)')
wipOutList(fid,fast.OutList);
wip(fid,[],'END of FAST input file (the word "END" must appear in the first 3 columns of this last line).')
wip(fid,[],'--------------------------------------------------------------------------------')


if fid~=1, fclose(fid); end
end


%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function wip(fid,param,descrip)
% write input file parameter
if ~any(size(param)) && ~ischar(param) 
    % do nothing if param = []
    % note: used ~any(size(param)) rather than isempty(param)
    %       so that unset parameters (size = [1 0]) will still 
    %       get through to the following elseif statements
elseif ischar(param)
    fprintf(fid,'%-16s ',param);  %output string
elseif isfloat(param)
    if numel(param)==1
        fprintf(fid,'%-16g ',param);  %output single number
    else
        str = sprintf(' %g',param);        %create list of numbers
        str = sprintf('"%s"',str(2:end));  %quote the list of numbers
        fprintf(fid,'%-16s ',str);          %output the quoted list
    end
end
fprintf(fid,'%s\n',descrip);
end

function wipcsv(fid,param,descrip)
% write input file parameter list as comma separated values
str = '';
for k = 1:length(param)
    str = strcat(str,sprintf('%d,',param(k)));
end
if length(str) > 1
    str(end) = [];
end
fprintf(fid,'%-16s %s\n',str,descrip);
end

function wiplst(fid,param,descrip)
% write input file parameter list (one per line)
for k = 1:length(param)
    fprintf(fid,'%-16s ',param{k});  %output string
    if k==1
        fprintf(fid,'%s',descrip);
    end
    fprintf(fid,'\n');
end
end

function wiptbl(fid,frmt,table,nrows)
for r=1:nrows
    for c=1:length(table)
        col = table{c};
        if iscell(col)
            param = col{r};
        else
            param = col(r);
        end
        fprintf(fid,frmt{c},param);
    end
    fprintf(fid,'\n');
end
end

function wipOutList(fid,outs)
for k=1:length(outs)
    fprintf(fid,'%s\n',outs{k});
end

%fprintf(fid,'%s',outs);
end