%% Include NuMAD code/tools in MATLAB working directories
addNumadPaths;
%% Input settings
dataFolder = 'exampleBlade';  %% Folder containing the data/files for the blade model
yamlFile = 'exampleBlade.yaml';  %% Name of .yaml file giving the initial design of the blade
defLoadsPath = 'defLoadsTable.mat';  %% Name of the loads table object file for use in maximum tip-deflection analysis
failLoadsPath = 'failLoadsTable.mat'; %% Name of the loads table object file for use in the failure/buckling/fatigue analysis
rccDataPath = 'rccdata.mat';  %% Name of object storing the rain cycle counting data resulting from runIEC
IECInput = 'IECInput.inp';  %% Name of file containing input for runIEC function.  Used here just for a few reference values.
elSize = 0.4;  %% The nominal size of an element in the ANSYS shell model.
DVarInit = zeros(1,4);  %% Inital values of design variables (in this case change in thickness of 17 main blade components.
lb = -10*ones(1,4); %% Lower bound values of design variables 
ub = 10*ones(1,4); %% Upper bound values of design variables
optAl = 'objective';  %% Tag specifying optimization algorithm to use: objective=evaluate objective, gradient=gradient-based optimization, particleswarm=particle swarm optimization
maxIt = 1;  %% Maximum iterations for the optimization algorithm
numProc = 1; %% Number of processors to use for parallel optimization
swarmSize = 4; %% Swarm size, if particle swarm optimization algorithm is used.

%% Set all the fields for the analysis configuration as described in Section 5 of documentation

%% Fields related to the maximum tip deflection analysis
allConfig.defConfig.ansys.meshFile = 'master.db';
allConfig.defConfig.ansys.analysisFileName = 'bladeAnalysis';
allConfig.defConfig.ansys.np = 1;
allConfig.defConfig.ansys.analysisFlags.resultantVSspan = 0;
allConfig.defConfig.ansys.analysisFlags.mass = 1;
allConfig.defConfig.ansys.analysisFlags.deflection = 1;

%% Fields related to material rupture, buckling and fatigue analysis
allConfig.failConfig.ansys.meshFile = 'master.db';
allConfig.failConfig.ansys.analysisFileName = 'bladeAnalysis';
allConfig.failConfig.ansys.np = 1;
allConfig.failConfig.ansys.rpm = 7.9; 
allConfig.failConfig.ansys.analysisFlags.resultantVSspan = 0;
allConfig.failConfig.ansys.analysisFlags.mass = 1;
allConfig.failConfig.ansys.analysisFlags.globalBuckling = 10;
allConfig.failConfig.ansys.nBucklingModes = 10;
allConfig.failConfig.ansys.analysisFlags.failure='TWSI';
allConfig.failConfig.ansys.analysisFlags.fatigue = ["HP_TE_FLAT","HP_TE_REINF","HP_TE_PANEL","HP_SPAR","HP_LE_PANEL","HP_LE","LP_LE","LP_LE_PANEL","LP_SPAR","LP_TE_PANEL","LP_TE_REINF","LP_TE_FLAT"];


%% Fields related to natural frequency analysis
allConfig.freqConfig.ansys.meshFile = 'master.db';
allConfig.freqConfig.ansys.analysisFileName = 'bladeAnalysis';
allConfig.freqConfig.ansys.np = 1;
allConfig.freqConfig.ansys.rpm = 7.9;
allConfig.freqConfig.ansys.nFrequencyModes = 10;

%% End input 

%% Include NuMAD code/tools in MATLAB working directories
homeDir = pwd;  % Save current directory
cd ..           % Return to NuMAD root directory
cd ..
addNumadPaths;  % Establish paths
cd(homeDir);    % Return to original directory

%% Initialize parallel working directories, based on number of processors
for i = 1:numProc
    folderName = ['parworker' num2str(i)];
    mkdir(folderName);
    mkdir([folderName '\NuMAD']);
    delete([folderName '\NuMAD\objectiveHistory.txt']);
    copyfile('airfoils',[folderName '\NuMAD\airfoils']);
    copyfile([dataFolder '\' rccDataPath],[folderName '\rccdata.mat']);
end

%%  Read the .yaml file into a NuMAD blade object
mainDir = pwd;
cd(dataFolder);
blade = BladeDef;
blade.readYAML(yamlFile);
blade.mesh = elSize;

IEC = IECDef(IECInput);
IEC.setAvgWindSpeed();

%%  Load the saved loads tables into the MATLAB workspace

if(contains(defLoadsPath,'.'))
    load(defLoadsPath,'defLoadsTable');
else
    defLoadsTable = {};
end

if(contains(failLoadsPath,'.'))
    load(failLoadsPath,'loadsTable');
else
    loadsTable = {};
end

cd(mainDir);

%% Execute optimization, or just objective evaluation based on tag stored in optAl.
if(contains(optAl,'objective'))  %% Evaluate the objective
    objVal = objectiveExample(DVarInit,blade,allConfig,defLoadsTable,loadsTable,IEC);
elseif(contains(optAl,'gradient'))  %% Run gradient-based optimization
    Aeq = [];
    beq = [];
    maxFunc = 6*maxIt*length(DVarInit);
    fdStep = 0.1*ones(1,length(DVarInit));
    typicalD = 0.5*ub;
    if(numProc > 1)
        upar = true;
    else
        upar = false;
    end
    options = optimoptions('fmincon','SpecifyObjectiveGradient',false,...
        'FiniteDifferenceType','central','FiniteDifferenceStepSize',fdStep,'TypicalX',typicalD,...
        'MaxIterations',maxIt,'MaxFunctionEvaluations',maxFunc,'UseParallel',upar);
    fun = @(DVar)objectiveExample(DVar,blade,allConfig,defLoadsTable,loadsTable,IEC);
    A = [];
    b = [];
    nonlcon = [];
    [DFinal,fval,exitflag,output] = fmincon(fun,DVarInit,A,b,Aeq,beq,lb,ub,nonlcon,options);
    fid = fopen('FinalDVar.txt','wt');
    for i = 1:length(DFinal)
        fprintf(fid,'%.16e  ',DFinal(i));
    end
    fclose(fid);
    outID = fopen('completeObjectiveHistory.txt','wt');
    for i = 1:numProc
        fileName = ['parworker' int2str(i) '\NuMAD\objectiveHistory.txt'];
        inID = fopen(fileName);
        terminate = 0;
        while(terminate == 0)
            line = fgetl(inID);
            if(line == -1)
                terminate = 1;
            else
                fprintf(outID,line);
                fprintf(outID,'\n');
            end
        end
        fclose(inID);
    end
    fclose(outID);
elseif(contains(optAl,'particleswarm')) %% Run Particle swarm optimization
    if(numProc > 1)
        upar = true;
    else
        upar = false;
    end
    options = optimoptions('particleswarm','SwarmSize',swarmSize,...
            'UseParallel',upar,'MaxIterations',maxIt);

    fun = @(DVar)objectiveExample(DVar,blade,allConfig,defLoadsTable,loadsTable,IEC);
    nVar = length(DVarInit);
    [DFinal,fval,exitflag,psOutput] = particleswarm(fun,nVar,lb,ub,options);
    fid = fopen('FinalDVar.txt','wt');
    for i = 1:length(DFinal)
        fprintf(fid,'%.16e  ',DFinal(i));
    end
    fclose(fid);
    outID = fopen('completeObjectiveHistory.txt','wt');
    for i = 1:numProc
        fileName = ['parworker' int2str(i) '\NuMAD\objectiveHistory.txt'];
        inID = fopen(fileName);
        terminate = 0;
        while(terminate == 0)
            line = fgetl(inID);
            if(line == -1)
                terminate = 1;
            else
                fprintf(outID,line);
                fprintf(outID,'\n');
            end
        end
        fclose(inID);
    end
    fclose(outID);
end