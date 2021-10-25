function runFASTusingSimulink(fastFilenameSFunc, modelName, YawPosBias)
%% ************************************************************************
% This function allows for a remote call of the Simulink production controller 
% The simulink model evaluates FAST in the Simulink workspace
% *************************************************************************

% Read FAST input file and set initial conditions
fst = readFastMain(fastFilenameSFunc);
runSimulation   % read in design values for the turbine operation

% pick the simulink controller model to use
if strcmp(modelName,'Production_Control_Yaw_FAST.mdl')
    % read in the assigned Yaw Parameters
    run('Yaw_parameters')
elseif strcmp(modelName,'Production_Control_FAST_noYaw.mdl')  
    % read in necessary files for the controller without yaw control
    run('Brake_parameters')
    run('Yaw_parameters')
elseif strcmp(modelName, 'Production_Control_FAST_NRELcenterSeeking.mdl')
    % read in the assigned Yaw Parameters
    run('Yaw_parameters')
    % read in the assigned generator brake parameters
    run('Brake_parameters')
    fst.TurbCtrl.HSSBrMode = 3; % this enables the Brake Torque to be added 
    fst.TurbCtrl.THSSBrDp = 0;
elseif strcmp(modelName, 'Production_Control_FAST_lookupTableYaw.mdl')
    % read in the assigned Yaw Parameters
    run('Yaw_parameters')
	run('Yaw_LUT')
    % read in the assigned generator brake parameters
    run('Brake_parameters')
    fst.TurbCtrl.HSSBrMode = 3; % this enables the Brake Torque to be added 
    fst.TurbCtrl.THSSBrDp = 0;
elseif strcmp(modelName, 'SwiftTurbine191217.mdl')
    run('Yaw_parameters')
    run('Brake_parameters')
    run('Pitch_parameters')
    run('Pitch_servo_parameters')
else
%     modelName = 'Production_Control_FAST.mdl';
end

% Save any changes to the FAST file and save as input_fast variable for the
% Read_FAST_Input function call
writeFastMain(fst,fastFilenameSFunc);
input_fast = fastFilenameSFunc;
run('Read_FAST_Input')  % Load FAST parameters in S-Function format
                        % NREL file: use same version as FAST


HydrBypassPitchRateMax = 20;  % pitch rate during hydraulic bypass (20-24 deg/s)
decimateToWorkspace = 10;


% % % this is the original definition of the startup procedure for Simulink
% % Turbine_State_ID = [0, STOP; 10, PAUSE; 20, RUN; TMax, RUN];

% forego the startup procedure in favor of a faster transient
Turbine_State_ID = [0, RUN; 10, RUN; 20, RUN; TMax, RUN];

ReleaseToFreeWheel = [0, 1; TMax, 1];

variableList = who;
for k=1:numel(variableList)
    assignin('base',variableList{k},eval(variableList{k}));
end

% assignin('base', 'input_fast', input_fast)
% assignin('base', 'fst', fst)
% evalin('base', 'Read_FAST_Input')
% evalin('base', 'V27params')
% assignin('base', 'input_fast', input_fast)
% assignin('base','HydrBypassPitchRateMax', HydrBypassPitchRateMax)
% assignin('base', 'decimateToWorkspace', decimateToWorkspace)
% assignin('base', 'Turbine_State_ID', Turbine_State_ID)
% assignin('base', 'ReleaseToFreeWheel', ReleaseToFreeWheel)
t = getCurrentTask;
if isempty(t)
    open(modelName);
end
SimOut = sim(modelName);
% ble: added a close because the optimization was freezing at a precomp
% analysis call.
bdclose

variableList = who;
for k=1:numel(variableList)
    assignin('base',variableList{k},eval(variableList{k}));
end


end