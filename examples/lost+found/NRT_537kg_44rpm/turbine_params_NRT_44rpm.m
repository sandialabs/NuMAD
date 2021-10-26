% Parameter script for NRT Production Controller
%
% Default Naming Conventions:
%   Un-decorated names, e.g. PitchScaleMax are in Engineering Units
%   Decorated names, e.g. PitchRateMaxPCNT are scaled [0-1];
% The unit: (%) is used to mean in the range [0-1]; i.e. 0-100%


% Information for calculating torque/speed-dependency:
% Rated Generator Speed: 1210 rpm
% Rated Power: 202.72 kW
%% From SNLV27_SFunc.fsm
% Rotor mass properties:
% 
%     Rotor Mass            (kg)         2970.678
%     Rotor Inertia         (kg-m^2)    80717.810 (KELHA: Includes 3xblades + hub)
%                                         Blade 1      Blade 2      Blade 3
%                                         -------      -------      -------
%     Mass                  (kg)          601.893      601.893      601.893
%     Second Mass Moment    (kg-m^2)    23440.352    23440.352    23440.352
%     First Mass Moment     (kg-m)       3109.373     3109.373     3109.373
%     Center of Mass        (m)             5.166        5.166        5.166
%% From  SNLV27.fst
% 27.6             GBRatio     - Gearbox ratio (-)
% 1165             HubMass     - Hub mass (kg)
% 0                TipMass(1)  - Tip-brake mass, blade 1 (kg)
% 0                TipMass(2)  - Tip-brake mass, blade 2 (kg)
% 0                TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]
% 8475.2           NacYIner    - Nacelle inertia about yaw axis (kg m^2)
% 50               GenIner     - Generator inertia about HSS (kg m^2)
% 617.216          HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)
% Inertia referred to generator: 
% Let's assume rotor inertia in SNLV27_SFunc.fsm file includes blades and hub; but 
% not the generator. THen inertia referred to the high-speed shaft would be:
% 80717.81/(27.6^2) + 50  = 156(kgm²)
% The FAST model has no information about gearbox and drive-shaft inertias
% - so presumably these are being considered negligible or are included in
% the rotor and generator inertias.


Ts_ProdCtrl=0.05; % (s)
CycleTime=0.1; % (s)

%% Generator Speed-Control Parameters
GenSpeedScaleMax = 1650; % (rpm)
GenSpeedScaleMin = 0; % (rpm)
GenSpeedSPHigh = 1210;
GenSpeedSPLow = 400;
GenSpeedSPVOG = 1500; % Set point for Vestas Overspeed Guard (VOG) Test(rpm)
GenSpeedSPSinT = 240; % Cycle period for Sin Test (s)
GenSpeedSPStepT = 60; % Step period for Step Test (s)
GenSpeedSPAmp = 1.0;  % Amplitude for Service Tests (0-1)

SCHL1_Kp = 3;   % Proportional Gain (dimensionless)
SCHL1_Ti = 7; % Integral-Action Time(s)
SCHL2_Kp = 9;   % Proportional Gain (dimensionless)
SCHL2_Ti = 7;   % Integral-Action Time(s)
SCHL1_FBFilterFreq1 = 1.0; % (Hz) tower mode
SCHL1_FBFilterFreq2 = 2.0; % (Hz) rotor flapwise asym


%% Pitch Manipulation Parameters:
PitchScaleMax = 88; % (°)
PitchScaleMin = -5; % (°)
PitchRateMax = 8; % (°/s)
PitchGainSchedMin = 5.48; % (°) From WindPitchRPM data.xlsx @ 13m/s
PitchGainSchedMax = 27.4; % (°) From WindPitchRPM data.xlsx @ 25m/s
PitchGainSchedDivisor = 5.23; % From WindPitchRPM data.xlsx


%% Torque Manipulation Parameters:
TorqueScaleMax = 1500; % (Nm) allow at least 5% for vibration control
TorqueScaleMin = 0; % (Nm)


%% Production Control Parameters
% ====== OEM V27 Rotor ======
% ProductionPitch = 1; % (°)  % ProductionPitchPCNT = (ProductionPitch-PitchScaleMin)/PitchSpan; % [0-1]
% ====== NRT Rotor ======
ProductionPitch = 0; % (°)  % ProductionPitchPCNT = (ProductionPitch-PitchScaleMin)/PitchSpan; % [0-1]
ProductionTorqueMax = 1520; % (Nm)
ProductionTorqueMin = 12;   % (Nm)
% ===== OEM V27 Rotor =====
% G = 27.5647; % Gearbox Ratio
% R = 13.5; % Rotor radius: (m)
% Cp = 0.4834; % Optimum Cp
% rho =  1.078; % Nominal Air Mass Density (kg/m³)
% lambda = 8.062; % Tip/Speed Ratio
% k=pi*rho*R^5*Cp/(2*lambda^3*G^3); % factor for load-torque Nm/(rad/s)²
% K = k* ((2*pi)/60)^2; % factor for load-torque Nm/(rpm)²
% ProductionTorqueCpMax = K* GenSpeedSPHigh^2; % Torque at GenSpeedSPHigh for Optimal power extraction.
% ProductionTorqueCpMax = 537; % Torque at GenSpeedSPHigh for Optimal power extraction (see above).
% ===== NRT Rotor =====
% K = 2.507e-4; % factor for load-torque Nm/(rpm)²
% ProductionTorqueCpMax = K* GenSpeedSPHigh^2; % Torque at GenSpeedSPHigh for Optimal power extraction.

ProductionTorqueCpMax = 367.05; % Torque at GenSpeedSPHigh for Optimal power extraction (see above).

%% FreeWheel Control Parameters
FreeWheelPitch = 45; % (°)
FreeWheelMinSpeed = 100; % (rpm)

%% Vibration-Control Parameters
VibrationControlOn=0; % True/False [0/1]
VC4_Kp = 0;
VC4_FBFilterFreq = 2.65; % (Hz)

%% Tower Resonance Avoidance Parameters
TowerResonanceCenter = 549.78; % (rpm) high-speed side, 19.945 rpm rotor
TowerResonanceHyst   =  64.78; % (rpm) high-speed side, 18.77 to 21.12 rpm (50% reduction expected)
TowerResCtrlRange    = 130.00; % (rpm) high-speed side
TowerResTorqueMax    =  50.0; % (Nm)

%% Production State-Control Parameters
%FreeWheelMinRPM=100; % rpm
%FreeWheelMaxRPM=400; % (rpm)
FreeWheelTime=15; % (s)
RunningUpTime=10; % (was 70) seconds
PowerUpTime=8; %seconds
RunningDownTime=10; % seconds
LowWindRunningDownTime=30; % seconds
% WTG States
EMERGENCY=-1;
STOP=0;
PAUSE=1;
RUN=2;



% % 
% % 
% % 
% % % Parameter script for RDV-V27 Production Controller
% % %
% % % Default Naming Conventions:
% % %   Un-decorated names, e.g. PitchScaleMax are in Engineering Units
% % %   Decorated names, e.g. PitchRateMaxPCNT are scaled [0-1];
% % % The unit: (%) is used to mean in the range [0-1]; i.e. 0-100%
% % 
% % 
% % % Information for calculating torque/speed-dependency:
% % % Rated Generator Speed: 1210 rpm
% % % Rated Power: 202.72 kW
% % %% From SNLV27_SFunc.fsm
% % % Rotor mass properties:
% % % 
% % %     Rotor Mass            (kg)         2970.678
% % %     Rotor Inertia         (kg-m^2)    80717.810 (KELHA: Includes 3xblades + hub)
% % %                                         Blade 1      Blade 2      Blade 3
% % %                                         -------      -------      -------
% % %     Mass                  (kg)          601.893      601.893      601.893
% % %     Second Mass Moment    (kg-m^2)    23440.352    23440.352    23440.352
% % %     First Mass Moment     (kg-m)       3109.373     3109.373     3109.373
% % %     Center of Mass        (m)             5.166        5.166        5.166
% % %% From  SNLV27.fst
% % % 27.6             GBRatio     - Gearbox ratio (-)
% % % 1165             HubMass     - Hub mass (kg)
% % % 0                TipMass(1)  - Tip-brake mass, blade 1 (kg)
% % % 0                TipMass(2)  - Tip-brake mass, blade 2 (kg)
% % % 0                TipMass(3)  - Tip-brake mass, blade 3 (kg) [unused for 2 blades]
% % % 8475.2           NacYIner    - Nacelle inertia about yaw axis (kg m^2)
% % % 50               GenIner     - Generator inertia about HSS (kg m^2)
% % % 617.216          HubIner     - Hub inertia about rotor axis [3 blades] or teeter axis [2 blades] (kg m^2)
% % % Inertia referred to generator: 
% % % Let's assume rotor inertia in SNLV27_SFunc.fsm file includes blades and hub; but 
% % % not the generator. THen inertia referred to the high-speed shaft would be:
% % % 80717.81/(27.6^2) + 50  = 156(kgm²)
% % % The FAST model has no information about gearbox and drive-shaft inertias
% % % - so presumably these are being considered negligible or are included in
% % % the rotor and generator inertias.
% % 
% % 
% % Ts_ProdCtrl=0.05; % (s)
% % CycleTime=0.1; % (s)
% % 
% % %% Generator Speed-Control Parameters
% % GenSpeedScaleMax = 1650; % (rpm)
% % GenSpeedScaleMin = 0; % (rpm)
% % GenSpeedSPHigh = 1210;
% % GenSpeedSPLow = 400;
% % GenSpeedSPVOG = 1500; % Set point for Vestas Overspeed Guard (VOG) Test(rpm)
% % GenSpeedSPSinT = 240; % Cycle period for Sin Test (s)
% % GenSpeedSPStepT = 60; % Step period for Step Test (s)
% % GenSpeedSPAmp = 1.0;  % Amplitude for Service Tests (0-1)
% % 
% % SCHL1_Kp = 3;   % Proportional Gain (dimensionless)
% % SCHL1_Ti = 7; % Integral-Action Time(s)
% % SCHL2_Kp = 9;   % Proportional Gain (dimensionless)
% % SCHL2_Ti = 7;   % Integral-Action Time(s)
% % SCHL1_FBFilterFreq1 = 1.0; % (Hz) tower mode
% % SCHL1_FBFilterFreq2 = 2.0; % (Hz) rotor flapwise asym
% % 
% % 
% % %% Pitch Manipulation Parameters:
% % PitchScaleMax = 88; % (°)
% % PitchScaleMin = -5; % (°)
% % PitchRateMax = 8; % (°/s)
% % PitchGainSchedMin = 5.48; % (°) From WindPitchRPM data.xlsx @ 13m/s
% % PitchGainSchedMax = 27.4; % (°) From WindPitchRPM data.xlsx @ 25m/s
% % PitchGainSchedDivisor = 5.23; % From WindPitchRPM data.xlsx
% % 
% % 
% % %% Torque Manipulation Parameters:
% % TorqueScaleMax = 1680; % (Nm) allow at least 5% for vibration control
% % TorqueScaleMin = 0; % (Nm)
% % 
% % 
% % %% Production Control Parameters
% % % ====== OEM V27 Rotor ======
% % % ProductionPitch = 1; % (°)  % ProductionPitchPCNT = (ProductionPitch-PitchScaleMin)/PitchSpan; % [0-1]
% % % ====== NRT Rotor ======
% % ProductionPitch = 0; % (°)  % ProductionPitchPCNT = (ProductionPitch-PitchScaleMin)/PitchSpan; % [0-1]
% % ProductionTorqueMax = 1520; % (Nm)
% % ProductionTorqueMin = 12;   % (Nm)
% % % ===== OEM V27 Rotor =====
% % % G = 27.5647; % Gearbox Ratio
% % % R = 13.5; % Rotor radius: (m)
% % % Cp = 0.4834; % Optimum Cp
% % % rho =  1.078; % Nominal Air Mass Density (kg/m³)
% % % lambda = 8.062; % Tip/Speed Ratio
% % % k=pi*rho*R^5*Cp/(2*lambda^3*G^3); % factor for load-torque Nm/(rad/s)²
% % % K = k* ((2*pi)/60)^2; % factor for load-torque Nm/(rpm)²
% % % ProductionTorqueCpMax = K* GenSpeedSPHigh^2; % Torque at GenSpeedSPHigh for Optimal power extraction.
% % % ProductionTorqueCpMax = 537; % Torque at GenSpeedSPHigh for Optimal power extraction (see above).
% % % ===== NRT Rotor =====
% % % K = 2edddddddddddd506e-4; % factor for load-torque Nm/(rpm)²
% % % ProductionTorqueCpMax = K* GenSpeedSPHigh^2; % Torque at GenSpeedSPHigh for Optimal power extraction.
% % ProductionTorqueCpMax = 366.9; % Torque at GenSpeedSPHigh for Optimal power extraction (see above).
% % 
% % %% FreeWheel Control Parameters
% % FreeWheelPitch = 45; % (°)
% % FreeWheelMinSpeed = 100; % (rpm)
% % 
% % %% Vibration-Control Parameters
% % VibrationControlOn=0; % True/False [0/1]
% % VC4_Kp = 0;
% % VC4_FBFilterFreq = 2.65; % (Hz)
% % 
% % %% Tower Resonance Avoidance Parameters
% % TowerResonanceCenter = 549.78; % (rpm) high-speed side, 19.945 rpm rotor
% % TowerResonanceHyst   =  64.78; % (rpm) high-speed side, 18.77 to 21.12 rpm (50% reduction expected)
% % TowerResCtrlRange    = 130.00; % (rpm) high-speed side
% % TowerResTorqueMax    =  50.0; % (Nm)
% % 
% % %% Production State-Control Parameters
% % %FreeWheelMinRPM=100; % rpm
% % %FreeWheelMaxRPM=400; % (rpm)
% % FreeWheelTime=15; % (s)
% % RunningUpTime=10; % (was 70) seconds
% % PowerUpTime=8; %seconds
% % RunningDownTime=10; % seconds
% % LowWindRunningDownTime=30; % seconds
% % % WTG States
% % EMERGENCY=-1;
% % STOP=0;
% % PAUSE=1;
% % RUN=2;





