addNumadPaths
cd('C:\data\BARphase2\NuMAD Training Blade Workshop\mainANSYSanalysis')

blade = BladeDef;
yamlFileName='myBlade';
blade.readYAML([yamlFileName '.yaml'])
%%


includeAdhesive=0;
meshData=blade.generateShellModel('ansys',includeAdhesive);

%%
load('defLoadsTable.mat')

%%
config.defConfig.meshFile = 'master.db';
config.defConfig.analysisFileName = 'bladeAnalysis';
config.defConfig.np = 1;
config.defConfig.analysisFlags.resultantVSspan = 1;
config.defConfig.analysisFlags.mass = 0;
config.defConfig.analysisFlags.deflection = 1;

ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,config.defConfig);

%%
X=blade.ispan;
x=ansysResult.deflection{1}(:,1);
y=ansysResult.deflection{1}(:,2); %flapwise deflection
z=ansysResult.deflection{1}(:,3);

plot(blade.ispan,y)

%%
load('loadsTable.mat')
%%
config.failConfig.meshFile = 'master.db';
config.failConfig.analysisFileName = 'bladeAnalysis';
config.failConfig.np = 1;
%config.failConfig.rpm = 7.9; 
config.failConfig.analysisFlags.resultantVSspan = 0;
config.failConfig.analysisFlags.mass = 1;
config.failConfig.analysisFlags.globalBuckling = 10;
config.failConfig.analysisFlags.failure='TWSI';
config.failConfig.analysisFlags.localFields=1;

loadsTable=loadsTable(1:1);
ansysResult = mainAnsysAnalysis(blade,meshData,loadsTable,config.failConfig)


%%

elNo=[3227];
coordSys='local'
myresult=extractFieldsThruThickness('plateStrains-all-1.txt',meshData,blade.materials,blade.stacks,blade.swstacks,elNo,coordSys)





scaleFactor=1e6; %Used to plot stress in MPa and strain in microstrain
close all;
figureNumbers=[4,5];
myYlabel='$y_3$ [mm]';
plotLocalFields(figureNumbers,myYlabel,myresult.(['element' num2str(elNo)]),scaleFactor,'k')

