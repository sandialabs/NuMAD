%% Add paths to NuMAD source
cd('..\..');
addNumadPaths;

%% Return to the example directory where the needed files are located
cd('examples\2022bladeWorkshop');

%% Create a blade object as an instance of BladeDef and read in data from the .yaml file

blade = BladeDef;
fileName = 'myBlade.yaml';
blade.readYAML(fileName);

%% Generate shell mesh
blade.mesh = 0.2;
includeAdhesive=0;
meshData=blade.generateShellModel('ansys',includeAdhesive);

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
Aeq = [];
beq = [];
A = [];
b = [];
nonlcon = [];
Xinit = [1.0,1.0];
lb = [0.25,0.25];
ub = [2.0,2.0];
maxIt = 15;
maxFunc = 75;
fdStep = [0.2,0.2];
options = optimoptions('fmincon','SpecifyObjectiveGradient',false,...
'FiniteDifferenceType','central','FiniteDifferenceStepSize',fdStep,...
'MaxIterations',maxIt,'MaxFunctionEvaluations',maxFunc);

%% Define the objective function
fun = @(XVar)objectiveFunction(XVar,blade,meshData,defLoadsTable,config);

%% Run optimization
[Xfinal,fval,exitflag,output] = fmincon(fun,Xinit,A,b,Aeq,beq,lb,ub,nonlcon,options);