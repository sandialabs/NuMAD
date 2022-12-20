%% Set initial options

DLCoptions = {'1.1','1.3','1.4','1.5','6.1','6.2','6.3'};   % Design load cases to run
inputFile = 'IECInput.inp';   % Name of the input file containing options for executing runIEC.
simFlag = 1;   %  Flag indicating whether to run fast (1), or to use existing output files (0) and generate output.
runAnsys = 0;  %  Flag indicating whether to generate max loads tables and run ansys analysis from runIEC results.

%% Establish global paths to NuMAD source code, and other modules as defined in addNumadPaths.m
%% This may need to be modified based on the path to the case/working directory, and the paths must be set for the specific machine/working environment in addNumadPaths.m

homeDir = pwd;  % Save current directory
cd ..           % Return to NuMAD root directory
cd ..
addNumadPaths;  % Establish paths
cd(homeDir);    % Return to original directory

%% Load data from IEC input file into IECDefObject

IEC = IECDef(inputFile);

%%  Load the blade object saved for the example blade

load exampleBladeObject.mat

%%  Write a NuMAD input file from the blade object

BladeDef_to_NuMADfile(blade,IEC.numadfn,'MatDBsi.txt');

%%  Execute runIEC function to perform aeroelastic analysis

IECOutput = runIEC(DLCoptions,simFlag,IEC);

if(runAnsys == 1)
    %% Make a working directory for the ANSYS analysis
    mkdir('NuMAD');
    copyfile('numad.nmd','NuMAD\numad.nmd');
    copyfile('MatDBsi.txt','NuMAD\MatDBsi.txt');
    copyfile('airfoils','NuMAD\airfoils');
    cd('NuMAD');
    
    %%  Generate the table of maximum loads seen throughout the aeroelastic simulation history

    fastGage = get_fast_gage(IEC.momentMaxRotation);
    maximumLoadsTable = FastLoads4ansys(IECOutput,fastGage,IEC);

    %%  Generate the load table corresponding to the moment of highest tip deflection throughout the aeroelastic simulation history

    maxDeflectionLoadsTable = getMaxDeflectionLoads(blade,IEC,IECOutput);

    %% Construct the configuration variables for ANSYS analysis
    
    %% Fields related to the maximum tip deflection analysis
    defConfig.meshFile = 'master.db';
    defConfig.analysisFileName = 'bladeAnalysis';
    defConfig.np = 1;
    defConfig.analysisFlags.resultantVSspan = 0;
    defConfig.analysisFlags.mass = 1;
    defConfig.analysisFlags.deflection = 1;

    %% Fields related to material rupture and buckling
    mainConfig.meshFile = 'master.db';
    mainConfig.analysisFileName = 'bladeAnalysis';
    mainConfig.np = 1;
    mainConfig.rpm = 10;
    mainConfig.analysisFlags.resultantVSspan = 0;
    mainConfig.analysisFlags.mass = 1;
    mainConfig.analysisFlags.globalBuckling = 10;
    mainConfig.nBucklingModes = 10;
    mainConfig.analysisFlags.failure='TWSI';
    
    %% Fields related to natural frequency analysis
    freqConfig.meshFile = 'master.db';
    freqConfig.analysisFileName = 'bladeAnalysis';
    freqConfig.np = 1;
    freqConfig.rpm = 10;
    freqConfig.nFrequencyModes = 10;
    
    %% Create the ANSYS model from the blade object
    blade.mesh = 0.1;
	includeAdhesive=0;
	meshData=blade.generateShellModel('ansys',includeAdhesive);
    
    %%  Run main ANSYS analysis for failure/buckling/fatigue
    
    mainAnalysisOut = mainAnsysAnalysis(blade,meshData,maximumLoadsTable,mainConfig,IEC)
    
    %%  Run ANSYS analysis for maximum tip deflection
    
    defAnalysisOut = mainAnsysAnalysis(blade,meshData,maxDeflectionLoadsTable,defConfig,IEC)
    
    %%  Run ANSYS analysis for frequency analysis
    freqAnalysisOut = getANSYSfrequency(freqConfig)
end

