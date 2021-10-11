close all; clear all; clc

%% ************************************************************************
% Define the starting blade file to read in design details
% *************************************************************************
bladeFile='..\bladeObject_IEA10p0-198-mk0p4.mat';

load(bladeFile)


%% ************************************************************************
% Set the optimization algorithm input variables
% *************************************************************************

ub = [1300; 500; 500; 500; 500; 400];
lb = [400; 200; 100; 100; 75; 0];

A=[ 0  0  0  -1  1  0;
    0  0  0  0  -1  1];
b = [0; 0];

% set genetic algorithm optimization iteration limits
% pop=36; gen=10;
pop=36; gen=2;


%% ************************************************************************
% Flags to alter the blade design from the loaded blade object
% *************************************************************************

% convert to a single shear web centered on the spar?
if false 
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

% change spar cap input span/thickness
if true
    sparCmp = 4;
    % ensure that the component changed is the spar cap
    if ~strcmp(blade.components(sparCmp).name,'spar')
        error('blade component is not the spar cap')
    end
    
    span = [0:0.25:1]';
    thickness = [100 150 150 100 10]';
    
    blade.components(sparCmp).cp = [span thickness];
end

% SNL3.0-148: Change spar cap material
sparCmp = find(contains({blade.components.name},'spar'));
if false
    % change to baseline CF (3mm); matlInd = 6
    % change to heavy textile CF (3mm); matlInd = 7
    sparMaterialIndex = 7;    
    blade.components(sparCmp).materialid = sparMaterialIndex;
end

% ----------------- change TE/LE reinf -----------------
if true
    % update LE reinforcement
    leReinfCmp = find(contains({blade.components.name},'le-reinf'));
    disp(blade.components(leReinfCmp))
    blade.components(leReinfCmp).cp(1,2)=60;
    blade.components(leReinfCmp).cp(2,1)=1;
    blade.components(leReinfCmp).cp(2,2)=60;
    disp(blade.components(leReinfCmp).cp)
    % update TE reinforcement
    teReinfCmp = find(contains({blade.components.name},'te-reinf'));
    disp(blade.components(teReinfCmp))
    blade.components(teReinfCmp).cp(1,2)=60;
    blade.components(teReinfCmp).cp(2,1)=1;
    blade.components(teReinfCmp).cp(2,2)=60;
    disp(blade.components(teReinfCmp).cp)
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

% save the preliminary blade file
% % newBladeName = [bladeFile(1:end-4) '_designIteration.mat'];
% % save(newBladeName, 'blade')

if useRestartFile
    restartFile='layupDesignCandidates.txt';
else
    restartFile='';
end


% run the optimization to minimize the mass
% % parfor ii = 1:64
% %     x = ub;
% %     f{ii} = delete_testParFor(blade, useParallel, x)
% % end

structOpt_mass_iea10p0_s1(blade,lb,ub,A,b,pop,gen,restartFile,useParallel)


% layupDesignCheck(xlsFile)

N=12; numelem=2.^linspace(1,8,N); aeSizes=numelem.^(-0.5);
% aeSizes=[0.05 0.1 0.15 0.20 0.25 0.3 0.4]';
% layupMeshCheck(xlsFile,aeSizes,0)

% layupDesignIterate(xlsFile)




