function [damage] = layupDesign_FASTfatigue(blade,runFASTfatigue,IEC)

% disp('Creating FAST Blade file using NuMAD and PreComp...')
% numad('numad.nmd','precomp',[1 3 2])
% copyfile('FASTBlade_precomp.dat','../FASTBlade_precomp.dat')
% disp('Newly created FAST Blade file copied to FAST simulation directory')
hm=pwd;
cd ..

% copy FAST fatigue output files from the main directory
if useParallel && ~runFASTfatigue
    % move up an additional directory to copy files
    
    workerOut = dir('out/IECDLC1p2NTM*.out');
    if(isempty(workerOut))
        outFiles = dir('../out/IECDLC1p2NTM*.out');    
        for ii = 1:length(outFiles)
            [success,~,~] = copyfile(fullfile(outFiles(ii).folder,outFiles(ii).name),'out/');
        end
    end
end
tic
disp('Running FAST/AeroDyn to verify blade OoP deflection (Please verify that proper DLC''s are set up in ''runIEC.m''...')
runIEC('1.2',runFASTfatigue,IEC);
toc

[Dminer,txt]=xlsread('IECDLC_1p2_F.csv');
materialID = txt(1,2:end);
channelID = txt(2:end,1);

varFlapMoment = contains(channelID,'M') & contains(channelID,'y');
varEdgeMoment = contains(channelID,'M') & contains(channelID,'x');
varStrainHP = contains(channelID,'eps') & contains(channelID,'HP');
varStrainLP = contains(channelID,'eps') & contains(channelID,'LP');


for ii = 1:length(materialID)
    damage.flap.(strrep(materialID{ii},'-','_'))=Dminer(varFlapMoment,ii);
    damage.edge.(strrep(materialID{ii},'-','_'))=Dminer(varEdgeMoment,ii);
    damage.HP.(strrep(materialID{ii},'-','_'))=Dminer(varStrainHP,ii);
    damage.LP.(strrep(materialID{ii},'-','_'))=Dminer(varStrainLP,ii);
end

cd(hm)
