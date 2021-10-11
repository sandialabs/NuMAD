function [FzMax,MzMax] = getMaxZForceMoment(runIECout)
    %% Get the maximum Z-force and Z-moment encountered at each station
    %% given the output structure from runIEC
    nGages = 10; %Assuming data is at 10 stations
    
    FzMax = -1e+10*ones(1,nGages);
    MzMax = -1e+10*ones(1,nGages);
    %Initialize
    for bb = 1:1
        FzBase = 'MaxRootFz';
        MzBase = 'MaxRootMz';
        for g =1:nGages  %Assuming data is at 10 stations
            chanName = [FzBase 'b' int2str(bb)];
            sAr = runIECout.(chanName);
            for i = 1:length(sAr)
                if(sAr(i).data > FzMax(g))
                    FzMax(g) = sAr(i).data;
                end
            end
            FzBase = ['MaxSpn' int2str(g) 'Fz'];
            
            chanName = [MzBase 'b' int2str(bb)];
            sAr = runIECout.(chanName);
            for i = 1:length(sAr)
                if(sAr(i).data > MzMax(g))
                    MzMax(g) = sAr(i).data;
                end
            end
            MzBase = ['MaxSpn' int2str(g) 'Mz'];
        end
    end
end

