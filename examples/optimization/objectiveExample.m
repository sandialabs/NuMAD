function [objVal] = objectiveExample(DVar,blade,config,defLoadsTable,loadsTable,IEC)
    %% Navigate to the appropriate parallel working directory
    mainDir = pwd;
    task = getCurrentTask;
    if(isempty(task))
        id = 1;
    else
        id = num2str(task.ID);
    end
    taskFolder = join([mainDir,"\parworker",id,"\NuMAD"],"");
    cd(taskFolder);
    
    %%  Assign values to the blade component thicknesses based on design variables DVar
    for i = 1:length(DVar)
        cThk = blade.components(i).cp(:,2);
        for j = 1:length(cThk)
            if(cThk(j) > 0.001)
                blade.components(i).cp(j,2) = max(0.1,(cThk(j) + DVar(i))); 
            end
        end
    end
    
    %%  Update all internal data within the blade based on new design values.
    blade.updateBlade

    %%  Generate FEA model for blade
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt');

    disp(' '); disp('Creating ANSYS model...')
    fprintf('Mesh size setting = %0.4f\n',blade.mesh)
    delete 'file.lock';
    layupDesign_ANSYSmesh(blade);
    
    %% Initialize values for objective evaluation and results.
    c1 = 60000;  %% Coefficient for constraint penalty terms, chosen based on approximate blade mass
    maxDeflection = 0;
    deflectionLimit = 20;
    maxFailIndex = 0;
    minLdFact = 100;
    maxFatDam = 0;
    
    %% Evaluate objective function as mass, plus series of penalty terms for constraints on displacement, failure, buckling, and flap frequency
    objVal = 0;
    
    %% Add displacement penalty
    if(isfield(config,'defConfig'))
        delete 'file.lock';
        ansysResult = layupDesign_ANSYSanalysis(blade,defLoadsTable,config.defConfig);
        maxDeflection = max(ansysResult.deflection{1},[],'all');
        objVal = objVal + c1*(maxDeflection/deflectionLimit)^2;
    end
    
    %% Add failure, buckling, fatigue penalties
    if(isfield(config,'failConfig'))
        delete 'file.lock';
        ansysResult = layupDesign_ANSYSanalysis(blade,loadsTable,config.failConfig,IEC);
        disp('postprocessing failure')
        if(isfield(config.failConfig.ansys.analysisFlags,'failure'))
            for i = 1:length(ansysResult.failure)
                if(ansysResult.failure{i} > maxFailIndex)
                    maxFailIndex = ansysResult.failure{i};
                end
            end
            objVal = objVal + c1*maxFailIndex^2;
        end
        disp('post.. buckling')
        if(isfield(config.failConfig.ansys.analysisFlags,'globalBuckling'))
            for i=1:length(ansysResult.globalBuckling)
                if(ansysResult.globalBuckling{i} < minLdFact)
                    minLdFact = ansysResult.globalBuckling{i};
                end
            end
            objVal = objVal + c1*(1/minLdFact)^2;
        end
        disp('...fatigue')
        if(isfield(config.failConfig.ansys.analysisFlags,'fatigue'))
            maxFatDam = 0;
            i1Max = length(ansysResult.fatigue);
            for i1 = 1:i1Max
                i1Dam = max(ansysResult.fatigue{i1}.fatigueDamage,[],'all');
                if(i1Dam > maxFatDam)
                    maxFatDam = i1Dam;
                end
            end
            objVal = objVal + c1*maxFatDam^2;
        end
    end
    
    %%  Add frequency penalty
    if(isfield(config,'freqConfig'))
        delete 'file.lock';
        rotorFreq = config.freqConfig.ansys.rpm/60;
        Freq = layupDesign_ANSYSfrequency(config.freqConfig);
        objVal = objVal + c1*(rotorFreq/(Freq(1)))^2;
    else
        Freq = [1,1];
    end
    
    %% Add total blade mass
    objVal = objVal + ansysResult.mass;
       
    %% Write results to objective history log
    fid = fopen('objectiveHistory.txt','a');
    fprintf(fid,'%s, Objective: %f, Mass: %f, Max Deflection: %f, \n',datestr(now),objVal,ansysResult.mass,maxDeflection);
    fprintf(fid,'Max Failure Index: %f \n',maxFailIndex);
    fprintf(fid,'Min Buckling Factor: %f \n',minLdFact);
    fprintf(fid,'Max Fatigue Damage: %f \n',maxFatDam);
    fprintf(fid,'Flap Frequency: %f, Edge Frequency: %f \n',Freq(1),Freq(2));
    fprintf(fid,'Design Variables: \n');
    for i = 1:length(DVar)
        fprintf(fid,'%f, ',DVar(i));
    end
    fprintf(fid,'\n');
    fclose(fid);
    cd(mainDir);
end

