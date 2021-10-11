%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
params.fstfn='NRT';   % without .fst
params.numadfn='NuMAD\numad.nmd';    % full path, including extension
params.fastsim = 'fast simulink';%'fast';%'adams';%
params.simulinkModel = 'Production_Control_FAST_NRELcenterSeeking.mdl';
params.simulinkModelFolder = '\\thor-storage\SCRATCH\blennis\DesignCodes\SWiFT_models\Controllers\controller_NRELcenterSeeking_newProdCtrl_TEST_NRT_44rpm';
params.operatingPoints=[5 11.7 25];  % cutin ratedspeed cutout
params.ws=[5:2:23];warning('ble: [5:2:23];')  % range of mean wind speeds for turbulent simulations
params.wd=[0]; % range of "wind direction" bias for look-up table simulations - programmed as yaw position
params.yaw = [0];    % intentional yaw misalignment, degrees (for DLC 1.1)
params.ratedSpeed=44; % rpm
params.lin=params.operatingPoints(1):1:params.operatingPoints(3);  % range of steady wind speeds for linearizations
% params.lin=10;  % range of steady wind speeds for linearizations
params.fast_path='\\thor-storage\scratch\blennis\DesignCodes\FAST_v7.02.00d\FAST.exe';
params.adams_path='call adams08r1 ru-user C:\DesignCodes\Compile_brr\FASTdll_AD_ADAMS\ADAMS08r1.dll';
params.turbsim_path='\\thor-storage\scratch\blennis\DesignCodes\TurbSim_v1.50\TurbSim.exe';
params.iecwind_path='\\thor-storage\scratch\blennis\DesignCodes\IECWind\IECWind.exe';
params.crunch_path='\\thor-storage\scratch\blennis\DesignCodes\Crunch_v3.00.00\Crunch.exe';
params.mbc_path='\\thor-storage\scratch\blennis\DesignCodes\MBC_v1.00.00a\Source';
params.sf_fat=1.380; % total fatigue safety factor
params.sf_uts=1.755; % total ultimate strength safety factor
params.sf_tow=1.755; % total tower clearance safety factor
params.numSeeds=6;warning('ble: 6;')  % number of seeds - number of 10-minute simulations - for turbulent simulations
params.delay=60;  % throw away this much simulated data at the beginning of each simulation (turbulent and otherwise)
params.SimTime=params.delay+600;warning('ble: params.delay+600;')  % total simulation time needed
params.NumGrid=10;warning('ble: [10]')  % number of grid points in turbsim 4-D wind field
params.Class=3; % turbine class: 1,2,3
params.TurbClass='C';  % turbulence class: A,B,C
params.BldGagNd=[1,2,3,4,5,6,7];  % vector of length 7; blade gage nodes (corresponding to aerodyn nodes) for moment (strain) gages in FAST computations 
% Material properties for fatigue analyses
matData(1).Name='E-LT-5500(UD)';
matData(1).E=29.38e9;
matData(1).b=10;  % from GL standard for uni-directional, epoxy laminate construction
matData(1).C=1000e6;
matData(2).Name='Newport 307 Carbon Prepreg (UD)';
matData(2).E=114.5e9;
matData(2).b=14;  % from GL standard for CFP, epoxy matrix
matData(2).C=1546e6;
matData(3).Name='SNL (Triax)';
matData(3).E=27.7e9;
matData(3).b=10;  % from GL standard for uni-directional, epoxy laminate construction
matData(3).C=700e6;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pause(10)