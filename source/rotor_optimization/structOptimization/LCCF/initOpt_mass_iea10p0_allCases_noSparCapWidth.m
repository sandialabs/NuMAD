function [] = initOpt_mass_iea10p0_allCases_noSparCapWidth(designNuMADfolder, bladeFile)
% close all; clear all; clc

%% ************************************************************************
% Define the starting blade file to read in design details
% *************************************************************************
% designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s1-fiberglass\NuMAD';
% designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s2-baselineCF\NuMAD';
% designNuMADfolder = 'C:\data\IEA10MW\IEA10p0-198-mk0p4-s3-heavyCF';

cd(designNuMADfolder)

% read in the starting blade file
% bladeFile='..\bladeObject_IEA10p0-198-mk0p4.mat';
load(bladeFile)


%% ************************************************************************
% Set the optimization algorithm input variables
% *************************************************************************

if contains(designNuMADfolder,'s1-fiberglass')
    % change to unidirection fiberglass (0.9mm);
    sparMaterialIndex = 1;        
    % set optimization variable limits
    ub = [100; 150; 175; 100; 50];
    lb = [20; 50; 50; 10; 1];
    % initial blade variable values
    sparSpan = [0 0.25 0.5 0.75 0.97 1]';
    sparLayers = [54 110 114 59 25 25]';  
    sparCapWidth = 900;

elseif contains(designNuMADfolder,'s2-baselineCF')
    % change to baseline CF (3mm);
    sparMaterialIndex = 6;    
    % set optimization variable limits    
    ub = [30; 75; 50; 40; 20];
    lb = [1; 10; 10; 1; 1];
    % initial blade variable values
    sparSpan = [0 0.25 0.5 0.75 0.97 1]';
    sparLayers = [1 32 23 11 3 3]'; 
    sparCapWidth = 400;
    
elseif contains(designNuMADfolder,'s3-heavyCF')
    % change to heavy textile CF (3mm);
    sparMaterialIndex = 7;    
    % set optimization variable limits  
    ub = [30; 75; 50; 40; 20];
    lb = [1; 10; 10; 1; 1];
    % initial blade variable values
    sparSpan = [0 0.25 0.5 0.75 0.97 1]';
    sparLayers = [1 49 21 20 1 1]'; 
    sparCapWidth = 400;
    
else
    error('ensure filename is correct')
end

% set optimization dependencies
A=[ 0  0  -1  1  0;
    0  0  0  -1  1];
b = [0; 0];

% set genetic algorithm optimization iteration limits
% pop=36; gen=10;
pop=36; gen=2;


%% ************************************************************************
% Flags to alter the blade design from the loaded blade object
% ************************************************************************* 
sparCmp = find(contains({blade.components.name},'spar'));
leReinfCmp = find(contains({blade.components.name},'le-reinf'));
teReinfCmp = find(contains({blade.components.name},'te-reinf'));
outerShellCmp = find(contains({blade.components.name},'outershell-bx'));
innerShellCmp = find(contains({blade.components.name},'innershell-bx'));
coreLEcmp = find(contains({blade.components.name},'core-le'));
coreTEcmp = find(contains({blade.components.name},'core-te'));
coreTipCmp = find(contains({blade.components.name},'core-tip'));
outerShellUDcmp = find(contains({blade.components.name},'outershell-ud'));
innerShellUDcmp = find(contains({blade.components.name},'innershell-ud'));

% ================= change spar cap material =================
if true
    blade.components(sparCmp).materialid = sparMaterialIndex;
end
% ================= change spar cap input span/thickness =================
if true
    blade.components(sparCmp).cp = [sparSpan sparLayers];    
    blade.sparcapwidth = sparCapWidth;
    disp(blade.components(sparCmp).name)
    disp(blade.components(sparCmp).cp)
end
% ================= change TE/LE reinf =================
if true
    % update LE reinforcement
    blade.components(leReinfCmp).cp(1,2)=50;
    blade.components(leReinfCmp).cp(2,1)=0.97;
    blade.components(leReinfCmp).cp(2,2)=50;
    blade.components(leReinfCmp).cp(3,1)=1;
    blade.components(leReinfCmp).cp(3,2)=50;
    disp(blade.components(leReinfCmp).name)
    disp(blade.components(leReinfCmp).cp)
    % update TE reinforcement
    blade.components(teReinfCmp).cp(1,2)=50;
    blade.components(teReinfCmp).cp(2,1)=0.97;
    blade.components(teReinfCmp).cp(2,2)=50;
    blade.components(teReinfCmp).cp(3,1)=1;
    blade.components(teReinfCmp).cp(3,2)=50;
    disp(blade.components(teReinfCmp).name)
    disp(blade.components(teReinfCmp).cp)
end
% ================= change panel shell thicknesses =================
if true
    shellLayers = 10;
    shellLayersTip = 5;
    % update panel outer shell thickness 
    blade.components(outerShellCmp).cp(1,1)=0;
    blade.components(outerShellCmp).cp(1,2)=shellLayers;
    blade.components(outerShellCmp).cp(2,1)=0.5;
    blade.components(outerShellCmp).cp(2,2)=shellLayers;
    blade.components(outerShellCmp).cp(3,1)=1;
    blade.components(outerShellCmp).cp(3,2)=shellLayersTip;
    disp(blade.components(outerShellCmp).name)
    disp(blade.components(outerShellCmp).cp)
    % update panel inner shell thickness 
    blade.components(innerShellCmp).cp(1,1)=0;
    blade.components(innerShellCmp).cp(1,2)=shellLayers;
    blade.components(innerShellCmp).cp(2,1)=0.5;
    blade.components(innerShellCmp).cp(2,2)=shellLayers;
    blade.components(innerShellCmp).cp(3,1)=1;
    blade.components(innerShellCmp).cp(3,2)=shellLayersTip;
    disp(blade.components(innerShellCmp).name)
    disp(blade.components(innerShellCmp).cp)
end
% ================= change panel core thicknesses =================
if true
    % update panel LE core thickness 
    blade.components(coreLEcmp).cp(1,2)=8;
    blade.components(coreLEcmp).cp(2,2)=4;
    blade.components(coreLEcmp).cp(3,2)=4;
    disp(blade.components(coreLEcmp).name)
    disp(blade.components(coreLEcmp).cp)
    % update panel TE core thickness 
    blade.components(coreTEcmp).cp(1,2)=8;
    blade.components(coreTEcmp).cp(2,2)=4;
    blade.components(coreTEcmp).cp(3,2)=4;
    disp(blade.components(coreTEcmp).name)
    disp(blade.components(coreTEcmp).cp)
    % update panel tip core thickness 
    blade.components(coreTipCmp).cp(1,2)=4;
    blade.components(coreTipCmp).cp(2,2)=4;
    disp(blade.components(coreTipCmp).name)
    disp(blade.components(coreTipCmp).cp)
end
% ================= change root transition thicknesses =================
if true
    shellUDLayers = 16;
    % update panel outer shell thickness 
    blade.components(outerShellUDcmp).cp(1,2)=shellUDLayers;
    blade.components(outerShellUDcmp).cp(2,2)=shellUDLayers;
    disp(blade.components(outerShellUDcmp).name)
    disp(blade.components(outerShellUDcmp).cp)
    % update panel inner shell thickness 
    blade.components(innerShellUDcmp).cp(1,2)=shellUDLayers;
    blade.components(innerShellUDcmp).cp(2,2)=shellUDLayers;
    disp(blade.components(innerShellUDcmp).name)
    disp(blade.components(innerShellUDcmp).cp)
end

disp(blade.materials(blade.components(sparCmp).materialid))
disp('Is this the correct material for the spar cap?')
pause(10)

%% ************************************************************************
% Call the optimization routine - first run or continuation of a run??
% *************************************************************************
useRestartFile = 0;     % flag to continue from a previous analysis or not
useParallel = true;     % flag to run the code in parallel or not

% finalize the geometry
blade.updateGeometry
blade.updateKeypoints
blade.updateBOM

if useRestartFile
    restartFile='layupDesignCandidates.txt';
else
    restartFile='';
end

structOpt_mass_IEA10p0_noSparCapWidth(blade,lb,ub,A,b,pop,gen,restartFile,useParallel)





