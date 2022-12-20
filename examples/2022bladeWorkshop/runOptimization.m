%% Add paths to NuMAD source
cd('..\..');
addNumadPaths;

%% Return to the example directory where the needed files are located
cd('examples\2022bladeWorkshop');

%% Create a blade object as an instance of BladeDef and read in data from the .yaml file

blade = BladeDef;
fileName = 'myBlade.yaml';
blade.readYAML(fileName);

%% Set the approximate element size for the shell mesh
blade.mesh = 0.2;

%% Import the loads table into the workspace
load('defLoadsTable.mat')

%% Set configuration parameters for the ANSYS analsysis
config.meshFile = 'master.db';
config.analysisFileName = 'bladeAnalysis';
config.np = 1;
config.analysisFlags.resultantVSspan = 1;
config.analysisFlags.mass = 0;
config.analysisFlags.deflection = 1;

%% Set basic optimization options settings
% Linear equality constraints (none)
Aeq = [];
beq = [];

% Linear inequality constraints (none)
A = [];
b = [];

% Nonlinear constraints (none)
nonlcon = [];

% Initial value for design variable, scaling factor for spar cap thickness
Xinit = 1.0;

% Lower and upper bounds for design variable
lb = 0.5;
ub = 2.0;

% Maximum iterations and function evaluations
maxIt = 10;
maxFunc = 30;

% Finite difference step
fdStep = 0.05;

% Optimizer options
options = optimoptions('fmincon','SpecifyObjectiveGradient',false,...
'FiniteDifferenceType','central','FiniteDifferenceStepSize',fdStep,...
'MaxIterations',maxIt,'MaxFunctionEvaluations',maxFunc);

%% Define the objective function
fun = @(XVar)objectiveFunction(XVar,blade,defLoadsTable,config);

%% Run optimization
[Xfinal,fval,exitflag,output] = fmincon(fun,Xinit,A,b,Aeq,beq,lb,ub,nonlcon,options);