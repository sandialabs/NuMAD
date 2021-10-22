close all; clear all; clc
snl2

%% ************************************************************************
% Define the starting blade file to read in design details
% *************************************************************************
bladeFile='..\bladeObject_DesignA6S0.mat';

load(bladeFile)


%% ************************************************************************
% Set the optimization algorithm input variables
% *************************************************************************
lb = [100; 0.05; 10; 0.15;  1;  1;  1;  1];
ub = [700; 0.14; 40; 0.19; 13; 13; 13; 13];

A=[ 0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0;
    0  0  0  0  0  0  0  0];

b=zeros(length(A),1);

% set genetic algorithm optimization iteration limits
pop=80; gen=50;
% pop=3; gen=2;


%% ************************************************************************
% Flags to alter the blade design from the loaded blade object
% *************************************************************************

% convert to a single shear web centered on the spar?
if true 
    sw2 = 2;  % component group ID of shear web 2
    sw1 = 1;  % component group ID of shear web 2
    hpLocation = {'0.5b-c'};    % desired location of shear web 1, HP side
    lpLocation = {'0.5b-c'};    % desired loaction of shear web 1, LP side
    blade.components = blade.components([blade.components.group] ~= sw2);
    for ii = find([blade.components.group] == sw1)
        blade.components(ii).hpextents = hpLocation;
        blade.components(ii).lpextents = lpLocation;
    end
end

% convert shear and spar to carbon?
if false 
    carbonID = find(strcmp({blade.materials.name},'Carbon(UD)'));
    % convert the shear web material to carbon using carbonID
    shearwebID = find(strcmp({blade.components.name},'sw-db'));
    for ii = shearwebID
        blade.components(ii).materialid = carbonID;
    end
    % convert the spar cap material to carbon using carbonID
    sparcapID = find(strcmp({blade.components.name},'spar'));
    for ii = sparcapID
        blade.components(ii).materialid = carbonID;
    end
end

% convert to all carbon blade?    
if false     
    carbon.blade = xlsBlade('C:\data\NRT\Rotor_Design\DesignA2S2_allcarbon.xlsx');
    blade.components = carbon.blade.components;
end    


%% ************************************************************************
% Call the optimization routine - first run or continuation of a run??
% *************************************************************************
useRestartFile = 0;     % flag to continue from a previous analysis or not
useParallel = true;     % flag to run the code in parallel or not


% finalize the geometry
blade.updateGeometry
blade.updateKeypoints
blade.updateBOM

% save the preliminary blade file
% % newBladeName = [bladeFile(1:end-4) '_designIteration.mat'];
% % save(newBladeName, 'blade')

if useRestartFile
    restartFile='layupDesignCandidates.txt';
else
    restartFile='';
end

% run the optimization to maximize the stiffness
structOpt_stiff(blade,lb,ub,A,b,pop,gen,restartFile,useParallel)


% layupDesignCheck(xlsFile)

N=12; numelem=2.^linspace(1,8,N); aeSizes=numelem.^(-0.5);
% aeSizes=[0.05 0.1 0.15 0.20 0.25 0.3 0.4]';
% layupMeshCheck(xlsFile,aeSizes,0)

% layupDesignIterate(xlsFile)




