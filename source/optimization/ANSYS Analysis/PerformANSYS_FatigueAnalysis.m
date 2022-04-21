function performANSYS_FatigueAnalysis(blade, config, iLoad, fid)
    %%%%%%%%%%%%%%%%%%%Outputs for fatigue analysis in MATLAB%%%%%%%%%%%%%%%%%
    fprintf(fid,'! BEGIN FATIGUE SCRIPT\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'esel,all\n');
    if iLoad==1 % Only save the following data for first load case
        %%% Material Properties %%%
        fprintf(fid,'mpwrite,Materials,txt,,\n');
        fprintf(fid,'/output,Strengths,txt,,\n');
        fprintf(fid,'TBLIST, ,ALL\n');
        fprintf(fid,'/output\n');

        %%% Section Properties %%%
        fprintf(fid,'/output, Sections,txt\n');
        fprintf(fid,'SLIST,,,,FULL\n');
        fprintf(fid,'/output\n');

        %%% Element Properties %%%
        fprintf(fid,'/output, Elements,txt\n');
        fprintf(fid,'elist,all,,,0,0 \n');
        fprintf(fid,'/output\n');
    end            
               
    if any(contains(lower(config.ansys.analysisFlags.fatigue),'all'))

        fprintf(fid,'allsel\n');
        fprintf(fid,'/prep7\n');
        fprintf(fid,'esel,all\n');
        fprintf(fid,'esel,u,type,,21  \n');
        fprintf(fid,'/POST1\n');
        fprintf(fid,'set,LAST\n');
        fprintf(fid,'RSYS,SOLU\n'); %Result in the element coordinate system

        %%% Element strains and curvatures %%%
        fprintf(fid, 'ALLSEL\n');
        fprintf(fid,'ETABLE, zcent,CENT,Z\n'); 
        fprintf(fid,'ETABLE, eps11,SMISC,9 \n');
        fprintf(fid,'ETABLE, eps22,SMISC,10 \n');
        fprintf(fid,'ETABLE, eps12,SMISC,11 \n');
        fprintf(fid,'ETABLE, kapa11,SMISC,12 \n');
        fprintf(fid,'ETABLE, kapa22,SMISC,13 \n');
        fprintf(fid,'ETABLE, kapa12,SMISC,14 \n');
        fprintf(fid,'ETABLE, gamma13,SMISC,15 \n');
        fprintf(fid,'ETABLE, gamma23,SMISC,16 \n');
        fprintf(fid,'/output,plateStrains-all-%s,txt\n',int2str(iLoad));
        fprintf(fid,'PRETAB,zcent,eps11,eps22,eps12,kapa11,kapa22,kapa12,gamma12,gamma13,gamma23\n');
        fprintf(fid, 'ETABLE,ERAS\n\n');            
    else
        TotalStations = numel(blade.ispan);
        %nSegments=numel(blade.keylabels)-1; %Number of blade regions 
        nSegments=numel(config.ansys.analysisFlags.fatigue);

        % Order of the segment names in segmentNamesReference
        % is very important. config.ansys.analysisFlags.fatigue can
        % be any order 
        segmentNamesReference=["HP_TE_FLAT","HP_TE_REINF","HP_TE_PANEL", "HP_SPAR","HP_LE_PANEL","HP_LE","LP_LE",...
                               "LP_LE_PANEL","LP_SPAR","LP_TE_PANEL","LP_TE_REINF","LP_TE_FLAT"];
        nsegmentNamesReference=numel(segmentNamesReference);          
        for i=1:nSegments
                    
            if ~ strcmpi(config.ansys.analysisFlags.fatigue(i),'webs')
                iSegment = find(strcmp(segmentNamesReference,config.ansys.analysisFlags.fatigue(i))==1);
                imax=iSegment+(TotalStations-2)*nsegmentNamesReference;
                fprintf(fid, 'ALLSEL\n');
                fprintf(fid, 'ASEL,S,AREA,,%i,%i,%i\n',iSegment,imax,nsegmentNamesReference); %Select all areas for a blade region (e.g HP_SPAR)
                fprintf(fid, 'ESLA, S\n',i,imax,nsegmentNamesReference); %Select the elements that are attatched to the selected areas
            else
                fprintf(fid, 'ALLSEL\n');
                fprintf(fid,'ESEL,S,SEC,,1000,100000   \n'); %Selects all web sections
                iSegment=nsegmentNamesReference+1;                         
            end
                    
            fprintf(fid,'/POST1\n');
            %% EMA added:
            fprintf(fid,'SET,FIRST\n');
            %% END
            fprintf(fid,'RSYS,SOLU\n'); %Result in the element coordinate system

            %%% Element strains and curvatures %%%
            fprintf(fid,'ETABLE, zcent,CENT,Z\n'); 
            fprintf(fid,'ETABLE, eps11,SMISC,9 \n');
            fprintf(fid,'ETABLE, eps22,SMISC,10 \n');
            fprintf(fid,'ETABLE, eps12,SMISC,11 \n');
            fprintf(fid,'ETABLE, kapa11,SMISC,12 \n');
            fprintf(fid,'ETABLE, kapa22,SMISC,13 \n');
            fprintf(fid,'ETABLE, kapa12,SMISC,14 \n');
            fprintf(fid,'ETABLE, gamma13,SMISC,15 \n');
            fprintf(fid,'ETABLE, gamma23,SMISC,16 \n');
            fprintf(fid,'/output,plateStrains-%s-%s,txt\n',int2str(iSegment),int2str(iLoad));
            fprintf(fid,'PRETAB,zcent,eps11,eps22,eps12,kapa11,kapa22,kapa12,gamma12,gamma13,gamma23\n');
            fprintf(fid, 'ETABLE,ERAS\n\n');
        end
    end

    fprintf(fid,'finish\n');
    fprintf(fid,'! END FATIGUE OUTPUT SCRIPT\n');
end