function output = layupDesign_FASTanalysis(blade,DLCoptions,runFASTsim,useParallel,IEC)
% run FAST simulations for ultimate loads and deflections
% NOTE: use layupDesign_FASTfatigue for DLC 1.2 fatigue calculations

hm=pwd;
cd ..
disp('Running FAST/AeroDyn to verify blade OoP deflection (Please verify that proper DLC''s are set up in ''runIEC.m''...')

if(useParallel && ~runFASTsim)
    workerOut = dir('out/IECDLC*.out');
    if(isempty(workerOut))
        outFiles = dir('../out/IECDLC*.out');    
        for ii = 1:length(outFiles)
            copyfile(fullfile(outFiles(ii).folder,outFiles(ii).name),'out/');
        end
    end
end

output=runIEC(DLCoptions,runFASTsim,IEC);
%% END

% for i=1:length(output.MaxOoPDefl)
%     x(i)=output.MaxOoPDefl(i).data;
% end
% designvar=max(x);

cd(hm)

end