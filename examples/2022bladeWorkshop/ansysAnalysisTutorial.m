%%Create a blade object from a yaml file
addNumadPaths
close all;

blade = BladeDef;
yamlFileName='myBlade';
blade.readYAML([yamlFileName '.yaml'])
%% create a shell model in ansys w/o adhesive

includeAdhesive=0;
meshData=blade.generateShellModel('ansys',includeAdhesive);

%% Load a previously build loadsTable
load('defLoadsTable.mat')

%% Set up configureation for deflection run
analysisConfig.meshFile = 'master.db';
analysisConfig.analysisFileName = 'bladeAnalysis';
analysisConfig.np = 1;
analysisConfig.analysisFlags.mass = 0;
analysisConfig.analysisFlags.deflection = 1;

ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,analysisConfig)

%% Plot the flapwise deflection
y=ansysResult.deflection{1}(:,2); %flapwise deflection

figure(1)
plot(blade.ispan,y)
hold on;
%%
%Modify Spar cap suction side
blade.components(3).name %Spar_cap_ss
figure(2)
plot(blade.components(3).cp(:,1),blade.components(3).cp(:,2))
hold on
blade.components(3).cp(:,2)=blade.components(3).cp(:,2)*1.1; %Increase component thickness by 10% throughout
plot(blade.components(3).cp(:,1),blade.components(3).cp(:,2))

blade.updateBlade %Run updateBlade with blade data was modified by user


%% Build a new ANSYS model; This time with previous mesh data. Then run deflection analysis again;
blade.generateShellModel('ansys',includeAdhesive,meshData);

ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,analysisConfig);


X=blade.ispan;
x=ansysResult.deflection{1}(:,1);
y=ansysResult.deflection{1}(:,2); %flapwise deflection
z=ansysResult.deflection{1}(:,3);

figure(1)
plot(blade.ispan,y,'--')
hold on;
%%
load('loadsTable.mat')
loadsTable
%%
failConfig.meshFile = 'master.db';
failConfig.analysisFileName = 'bladeAnalysis';
failConfig.np = 1;
failConfig.analysisFlags.mass = 1;
failConfig.analysisFlags.globalBuckling = 10;
failConfig.analysisFlags.failure='TWSI';
failConfig.analysisFlags.localFields=1;

ansysResult = mainAnsysAnalysis(blade,meshData,loadsTable,failConfig)


%% Local Fields Example

elNo=[3183,3194]; %HP spar cap (station 12)
coordSys='local'
myresult=extractFieldsThruThickness('plateStrains-all-2.txt',meshData,blade.materials,blade.stacks,blade.swstacks,elNo,coordSys)

%Plot the stress and strains though-the-thickness
scaleFactor=1e6; %Used to plot stress in MPa and strain in microstrain

figureNumbers=[4,5];
myYlabel='$y_3$ [mm]';
plotLocalFields(figureNumbers,myYlabel,myresult.(['element' num2str(elNo(1))]),scaleFactor,'k')
hold on;
plotLocalFields(figureNumbers,myYlabel,myresult.(['element' num2str(elNo(2))]),scaleFactor,'r')
