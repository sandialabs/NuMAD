function [designvar] = layupDesign_ANSYSanalysis(blade,loadsTable,config,IEC)
    anFlagNames = fieldnames(config.ansys.analysisFlags);
    
    global ansysPath
    ansys_product = 'ANSYS';

    if isfield(config.ansys.analysisFlags,'imperfection') && ~isempty(config.ansys.analysisFlags.imperfection) &&config.ansys.analysisFlags.globalBuckling ==0
        error('Specify number of buckling modes when performing nonlinear buckling')
    end

    % Original mesh file to analize
    if isfield(config.ansys,'meshFile') && ~isempty(config.ansys.meshFile)
        %Do nothing
    else
        config.ansys.meshFile='master';
    end
    
    % File name base name for ansys analysis files
    if isfield(config.ansys,'analysisFileName') && ~isempty(config.ansys.analysisFileName)
        ansysFilename = [config.ansys.analysisFileName];
    else
        ansysFilename=['FEmodel'];
    end
    % Number of CPUs to use
    if isfield(config.ansys,'np') && ~isempty(config.ansys.np)
        if config.ansys.np <1
            error('config.ansys.np must be greater than zero')
        else
            ncpus=config.ansys.np;
        end
    else
        ncpus=1; %Default value
    end
    %Initialize
    for i=1:length(anFlagNames)
        
        if strcmpi('globalBuckling',anFlagNames{i}) || strcmpi('resultantVSspan',anFlagNames{i}) || strcmpi('deflection',anFlagNames{i}) || strcmpi('mass',anFlagNames{i});
            if config.ansys.analysisFlags.(anFlagNames{i})~=0 ;
                designvar.(anFlagNames{i}) = cell(1,length(loadsTable));
            end
        else strcmpi('localBuckling',anFlagNames{i}) || strcmpi('failure',anFlagNames{i}) || strcmpi('fatigue',anFlagNames{i}) ||strcmpi('imperfection',anFlagNames{i}) || strcmpi('mass',anFlagNames{i});
            if ~isempty(config.ansys.analysisFlags.(anFlagNames{i}));
                designvar.(anFlagNames{i}) = cell(1,length(loadsTable));
            end
        end
    end
    
    if ~exist('designvar')
        error('no analyses are configured in config. Edit config.')
    end
    
    if isfield(config.ansys.analysisFlags,'localBuckling') && ~isempty(config.ansys.analysisFlags.localBuckling)|| isfield(config.ansys.analysisFlags,'fatigue') && ~isempty(config.ansys.analysisFlags.fatigue)
        [isoorthoInModel,compsInModel,SkinAreas,app] = getMatrialLayerInfoWithOutGUI(blade);

        bladeMatNames=cell(numel(blade.materials),1);
        for iMat=1:numel(blade.materials)
            bladeMatNames{iMat}=blade.materials(iMat).name;
        end

        matPointer=zeros(numel(isoorthoInModel),1);
        for iMat=1:numel(isoorthoInModel)
           ansysMPnumber = find(strcmp(isoorthoInModel(iMat),bladeMatNames)==1);
           matPointer(iMat)=ansysMPnumber;
        end
        
        ansysBladeMaterials=blade.materials(matPointer);
    end
    for iLoad=1:length(loadsTable)
        %% ************************************************************************
        % ================= APPLY BUCKLING LOADS TO FEA MESH =================
        forcefilename='forces';
        

        beamForceToAnsysShell('map3D_fxM0','NLIST.lis',loadsTable{iLoad},strcat(forcefilename,'.src'));

        disp('Forces mapped to ANSYS model')

        %% ************************************************************************
        % ================= PERFORM LINEAR STATIC ANALYSIS =================
        % run buckling computations in ansys
        disp(' '); disp('Running ANSYS analysis...')

        script_name='ansysAnalysis.mac';
        script_out='ansysAnalysisEcho.out';

        fid=fopen(script_name,'w+');

        fprintf(fid,'/NERR,,99999999\n');
        fprintf(fid,'resume,master,db\n');
        fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);   %From master, change the jobname
        %fprintf(fid,'resume\n');

        fprintf(fid,'! BEGIN LINEAR STATIC SCRIPT\n');
        fprintf(fid,'esel,all\n');
        fprintf(fid,'/prep7\n');
        fprintf(fid,'fdel,all\n');
        fprintf(fid,'/input,%s,src\n',forcefilename);

        %Linear Static Analysis
        fprintf(fid,'/solu\n');
        
        fprintf(fid,'antype,static\n');
%         if nonlinear
%             fprintf(fid,'nlgeom,1\n'); %%%%%%%%%%%%%%%%%%%%%%%% TEMP
%             fprintf(fid,'OUTRES,all,ALL\n');%%%%%%%%%%%%%%%%%%%%%%%% TEMP
%         else
            fprintf(fid,'pstres,on\n');
%         end
        fprintf(fid,'irlf,-1\n');
        
        fprintf(fid,'bcsoption,,incore\n');
        
        fprintf(fid,'solve\n');
        fprintf(fid,'finish\n');
        
        %Only compute mass on the first load case
        if iLoad==1 && isfield(lower(config.ansys.analysisFlags),'mass') && ~isempty(config.ansys.analysisFlags.mass)&& config.ansys.analysisFlags.mass~=0
            %Get Mass Here
            fprintf(fid,'*GET, Z_mass, ELEM, 0, MTOT, X\n');
            fprintf(fid,'/output, mass,txt\n');
            fprintf(fid,'*status,Z_mass\n');
            fprintf(fid,'/output\n');
            fprintf(fid,'finish\n');
        end
 
 %% ************************************************************************
%================= PERFORM Deflection ANALYSIS =================

        if isfield(config.ansys.analysisFlags,'deflection') && config.ansys.analysisFlags.deflection~=0
            fprintf(fid,'/POST1\n');
            fprintf(fid,'set,last\n');
            fprintf(fid,'RSYS,0\n');  %global coordinates
            deflectionFilename = 'results_deflection';

            fprintf(fid,'seltol,0.05\n');
            for i=1:numel(blade.ispan)
                fprintf(fid,'*CFOPEN, %s,out\n',[deflectionFilename '-' int2str(i)]);
                fprintf(fid,'ESEL,S,SEC,,1,999   \n');    %Selects aero shell only
                fprintf(fid,'nsle,S,   \n');    %Selects aero shell only
                fprintf(fid,'nsel,r,loc,z,%f  \n',blade.ispan(i));
    %                 fprintf(fid,'nsll,s,,\n');
                if i==numel(blade.ispan)
                    fprintf(fid,'nsel,u,node,,z_master_node_number\n');
                end
    %                 fprintf(fid,'nplot\n');
                fprintf(fid,'*GET, NsectionNodes, NODE,0,COUNT   !Get the number of nodes in the set\n');
                fprintf(fid,'*GET, node_num, NODE,0,NUM,MIN        !Get the smallest number node in the set\n');
                fprintf(fid,'*DO, i, 1, NsectionNodes                 !loop through all nodes in cross section\n');

                fprintf(fid,'*GET, xpos, NODE,node_num,loc,X\n');
                fprintf(fid,'*GET, ypos, NODE,node_num,loc,Y\n');
                fprintf(fid,'*GET, zpos, NODE,node_num,loc,Z\n');

                fprintf(fid,'*GET, u1, NODE,node_num,U,X\n');
                fprintf(fid,'*GET, u2, NODE,node_num,U,Y\n');
                fprintf(fid,'*GET, u3, NODE,node_num,U,Z\n');
                fprintf(fid,' *VWRITE,node_num,xpos,ypos,zpos,u1,u2,u3\n');
                fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
                fprintf(fid,'node_num=NDNEXT(node_num)             !Get the next higher node number in the set\n');
                fprintf(fid,'*ENDDO\n');

                fprintf(fid,'*CFCLOS\n');
                fprintf(fid,'\n \n \n');        
            end
            fprintf(fid,'finish\n');
        end     

        % calculate face stresses for wrinkling
        if isfield(config.ansys.analysisFlags,'localBuckling') && ~isempty(config.ansys.analysisFlags.localBuckling)&& ~(isfield(config.ansys.analysisFlags,'imperfection') && ~isempty(config.ansys.analysisFlags.imperfection))
             %Check for wrinkling here in a linear analysis
            [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,config.ansys.analysisFlags.localBuckling);
        end

        %%% Output resultant force and moments to file
        if isfield(config.ansys.analysisFlags, 'resultantVSspan') && config.ansys.analysisFlags.resultantVSspan~=0
            fprintf(fid,'/POST1\n');
            fprintf(fid,'set,LAST\n');
            fprintf(fid,'RSYS,0\n'); %global coordinates

%             fprintf(fid,'seltol,0.05\n');
%             fprintf(fid,'*CFOPEN, resultantVSspan,txt\n');
%             for i=1:numel(blade.ispan)
%                 fprintf(fid,'nsel,s,loc,z,0,%f  \n',blade.ispan(i));
% 
%                 if i==numel(blade.ispan)
%                     fprintf(fid,'nsel,u,node,,z_master_node_number\n');
%                 end
% 
%                 fprintf(fid,'spoint,0,%f,%f,%f\n',blade.sweep(i),blade.prebend(i),blade.ispan(i));
%                 fprintf(fid,'nplot\n');
%                 fprintf(fid,'FSUM\n');
%                 fprintf(fid,'*GET, F1, FSUM, 0, ITEM,FX\n');
%                 fprintf(fid,'*GET, F2, FSUM, 0, ITEM,FY\n');
%                 fprintf(fid,'*GET, F3, FSUM, 0, ITEM,FZ\n');
%                 fprintf(fid,'*GET, M1, FSUM, 0, ITEM,MX\n');
%                 fprintf(fid,'*GET, M2, FSUM, 0, ITEM,MY\n');
%                 fprintf(fid,'*GET, M3, FSUM, 0, ITEM,MZ\n');
% 
%                 fprintf(fid,'*VWRITE,%f,F1,F2,F3,M1,M2,M3\n',blade.ispan(i));
%                 fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
% 
%                 fprintf(fid,'\n \n \n');
%             end
%             fprintf(fid,'*CFCLOS\n');
%             fprintf(fid,'finish\n');            
            
            fprintf(fid,'/post1\n');
            fprintf(fid,'elsize=%f\n',blade.mesh);
            fprintf(fid,'nz=nint(%f/elsize) !Integer number of points to output resultant loads\n',blade.ispan(end));

            fprintf(fid,'zloc=0\n');
            fprintf(fid,'delta=0.1\n');
            fprintf(fid,'*CFOPEN, resultantVSspan,txt\n');
            fprintf(fid,'*do,I,1,nz+1\n');
            fprintf(fid,'allsel\n');
            fprintf(fid,'nsel,s,loc,z,0,zloc+delta\n');
            fprintf(fid,'spoint,0,0,0,zloc\n');
            fprintf(fid,'!nplot\n');
            fprintf(fid,'FSUM\n');
            fprintf(fid,'*GET, F1, FSUM, 0, ITEM,FX\n');
            fprintf(fid,'*GET, F2, FSUM, 0, ITEM,FY\n');
            fprintf(fid,'*GET, F3, FSUM, 0, ITEM,FZ\n');
            fprintf(fid,'*GET, M1, FSUM, 0, ITEM,MX\n');
            fprintf(fid,'*GET, M2, FSUM, 0, ITEM,MY\n');
            fprintf(fid,'*GET, M3, FSUM, 0, ITEM,MZ\n');

            fprintf(fid,'*VWRITE,zloc,F1,F2,F3,M1,M2,M3\n');
            fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
            fprintf(fid,'zloc=zloc+elsize\n');
            fprintf(fid,'*ENDDO\n');
            fprintf(fid,'*CFCLOS\n');

            fprintf(fid,'fini\n');

        end
        %% ************************************************************************
        % ================= PERFORM FATIGUE ANALYSIS =================
        if isfield(config.ansys.analysisFlags,'fatigue') && ~isempty(config.ansys.analysisFlags.fatigue)
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

                % %%% Element strains and curvatures %%%

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
    %             nSegments=numel(blade.keylabels)-1; %Number of blade regions 
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

                    % %%% Element strains and curvatures %%%


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
%% ************************************************************************
% ================= CREAT LOCAL FIELD RESULTS FOR MATLAB =================

   if isfield(config.ansys.analysisFlags,'localFields') && ~isempty(config.ansys.analysisFlags.localFields)
            %%%%%%%%%%%%%%%%%%%Outputs for fatigue analysis in MATLAB%%%%%%%%%%%%%%%%%
            fprintf(fid,'! BEGIN LOCAL FIELD SCRIPT\n');
            fprintf(fid,'allsel\n');
            fprintf(fid,'/prep7\n');
            fprintf(fid,'esel,all\n');


            
            fprintf(fid,'esel,u,type,,21  \n');
            fprintf(fid,'/POST1\n');
            fprintf(fid,'set,LAST\n');
            fprintf(fid,'RSYS,SOLU\n'); %Result in the element coordinate system
               

            % %%% Element Stress %%%

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

%             fprintf(fid,'ETABLE, zcent,CENT,Z\n'); 
%             fprintf(fid,'ETABLE, N11,SMISC,1 \n');
%             fprintf(fid,'ETABLE, N22,SMISC,2 \n');
%             fprintf(fid,'ETABLE, N12,SMISC,3 \n');
%             fprintf(fid,'ETABLE, M11,SMISC,4 \n');
%             fprintf(fid,'ETABLE, M22,SMISC,5 \n');
%             fprintf(fid,'ETABLE, M12,SMISC,6 \n');
%             fprintf(fid,'ETABLE, Q13,SMISC,7 \n');
%             fprintf(fid,'ETABLE, Q23,SMISC,8 \n');
%             %fprintf(fid,'/output,plateExamplePlateForces-all-%s,txt\n',int2str(iLoad));
%             fprintf(fid,'/output,plateForces-all-%s,txt\n',int2str(iLoad));
%             fprintf(fid,'PRETAB,zcent,N11,N22,N12,M11,M22,M12,Q12,Q13,Q23\n');
%             fprintf(fid, 'ETABLE,ERAS\n\n');
            
            
            fprintf(fid,'finish\n');
         
        end

        %% ************************************************************************
        % ================= PERFORM FAILURE ANALYSIS =================

        % Initialize GUI commands from batch operation to identify maxima
        if isfield(config.ansys.analysisFlags,'failure') && ~isempty(config.ansys.analysisFlags.failure)
            fprintf(fid,'! BEGIN FAILURE SCRIPT\n');
            fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            fprintf(fid,'!Add for PLESOL and *get,findex,PLNSOL,0,MAX to work');
            fprintf(fid,'/BATCH  \n');
            fprintf(fid,'/COM,ANSYS RELEASE Release 18.1      BUILD 18.1      UP20170403       15:49:08\n');
            fprintf(fid,'/GRA,POWER\n ');
            fprintf(fid,'/GST,ON\n ');
            fprintf(fid,'/PLO,INFO,3\n ');
            fprintf(fid,'/GRO,CURL,ON\n ');
            fprintf(fid,'/CPLANE,1   \n ');
            fprintf(fid,'/REPLOT,RESIZE  \n ');
            fprintf(fid,'WPSTYLE,,,,,,,,0\n ');
            fprintf(fid,'/SHOW\n ');
            fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
            fprintf(fid,'/POST1\n');
            fprintf(fid,'set,last\n');
            fprintf(fid,'allsel\n');
            fprintf(fid,'esel,all\n');

            fprintf(fid,'RSYS,LSYS \n'); %Result in the layer coordinate system
            fprintf(fid,'layer,fcmax\n');

            if ~contains(upper(config.ansys.analysisFlags.failure),{'PUCK' 'LARC03' 'LARC04'})
                %Do this for the failure criteria that do not distinguish between fiber
                %and matrix failure
                fc=upper(config.ansys.analysisFlags.failure);
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,findex,PLNSOL,0,MAX\n');
            else
                %Do this for the failure criteria that do distinguish between fiber
                %and matrix failure
                switch upper(config.ansys.analysisFlags.failure)
                    case 'PUCK'
                        %Fiber Failure
                        fc='PFIB';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                        %Matrix Failure
                        fc='PMAT';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
                    case 'LARC03'
                        %Fiber Failure
                        fc='L3FB';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                        %Matrix Failure
                        fc='L3MT';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
                    case 'LARC04'
                        %Fiber Failure
                        fc='L4FB';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                        %Matrix Failure
                        fc='L4MT';
                        fprintf(fid,'FCTYP,add,%s\n',fc);
                        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                        fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
                end
                %Report the higher of the fiber failure index or the matrix
                fprintf(fid,'*IF, Ffindex, GT,Mfindex, THEN\n');
                fprintf(fid,'findex=Ffindex\n');
                fprintf(fid,'*ELSE\n');
                fprintf(fid,'findex=Mfindex\n');
                fprintf(fid,'*ENDIF\n');
            end
            failureFilename = 'results_failure';
            fprintf(fid,['/output,' failureFilename ',out\n']);
            fprintf(fid,'*status,findex\n');
            fprintf(fid,'/output\n');
            fprintf(fid,'finish\n');
            fprintf(fid,'! END FAILURE SCRIPT\n');
        end


        %% ************************************************************************
        % ================= PERFORM BUCKLING ANALYSIS =================

        %Linear Buckling Analysis
        if isfield(config.ansys.analysisFlags,'globalBuckling') && config.ansys.analysisFlags.globalBuckling >0
            bucklingFilename = 'results_buckling';
            fprintf(fid,'! BEGIN BUCKLE MACRO\n');
            fprintf(fid,'allsel\n');
            fprintf(fid,'/solu\n');
            fprintf(fid,'irlf,-1\n');
            fprintf(fid,'pstres,on\n');
            fprintf(fid,'antype,buckle\n');
            fprintf(fid,strcat('bucopt,lanb,',int2str(config.ansys.analysisFlags.globalBuckling),',,,RANGE\n'));
            %fprintf(fid,strcat('MXPAND,',int2str(nmodes),',0,0,1\n'), nmodes); % Requiered for element stress/strain, etc..
            fprintf(fid,'solve\n');
            fprintf(fid,'finish\n');
            fprintf(fid,'/post1\n');
            fprintf(fid,['/output,' bucklingFilename ',out\n']);
            fprintf(fid,'set,list\n');
            fprintf(fid,'/output\n');
            fprintf(fid,'finish\n');
            fprintf(fid,'! END BUCKLE MACRO\n');   
        elseif isfield(config.ansys.analysisFlags,'globalBuckling') && config.ansys.analysisFlags.globalBuckling <0
            error('config.ansys.analysisFlags.globalBuckling must be greater than or equal to zero')
        end
        % ANSYSoutputByBladeRegion(blade,fid)


        %% ************************************************************************
        % ================= SEND COMMANDS TO ANSYS =================
        fclose(fid);

        ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',...
            ansysPath,ansys_product,script_name,script_out,int2str(ncpus));
%         KMP_STACKSIZE is 512k by default. This is not enough therefore SET
%         KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
%         simulations. EC
        % 
        while 1    
            [status,~] = system(ansys_call)  % the windows system call to run the above ansys command
            if status==0
                disp(' ')
                disp('ANSYS analysis completed')
                %delete(script_name);
                break;
            end
            fprintf('%s: Waiting for ANSYS analysis steps...\n',datestr(now))
            pause(3);
        end
  %  MATLAB POST PROCESS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%% ************************************************************************
% ================= READ MASS RESULTS INTO MATLAB =================
        if iLoad==1 && isfield(lower(config.ansys.analysisFlags),'mass') && ~isempty(config.ansys.analysisFlags.mass)&& config.ansys.analysisFlags.mass~=0
            designvar.mass=read_1_ANSYSoutput('mass.txt');
            delete mass.txt
        end

    %% ************************************************************************
% ================= READ DEFLECTION RESULTS INTO MATLAB =================       
        if isfield(config.ansys.analysisFlags,'deflection') && config.ansys.analysisFlags.deflection~=0

            nSpan=length(blade.ispan);
            data=zeros(nSpan,6);  %u1Avg,u2Avg,u3Avg,0,theta2,theta3
            for iSpan=1:nSpan
                fileName=[deflectionFilename '-' int2str(iSpan) '.out'];
                temp=txt2mat(fileName);  %node number, Displacements of the nodes on the crossection
                delete(fileName);

                %Displacement
                for k=1:3
                   data(iSpan, k) =  mean(temp(:,k+4));
                end

                [nNode,~]=size(temp);  %Number of nodes

                
                [xmax,LE]=max(temp(:,2)); %Find the max x location (column 2) 
                [xmin,TE]=min(temp(:,2)); %Find the max x location (column 2) 


                [ymax,LP]=max(temp(:,3)); %Find the max y location (column 3) 
                [ymin,HP]=min(temp(:,3)); %Find the max y location (column 3) 
%                close all;
%                plot(temp(:,2),temp(:,3),'ok')
%                hold on;

                P = temp(LE,2:4); %x,y,z coordinates for point P at leading edge
                Q = temp(TE,2:4); %x,y,z coordinates for point Q at trailing edge
                PQ=P-Q;

%                quiver(Q(1),Q(2),PQ(1),PQ(2));


%                plot(temp(:,2)+temp(LE,2),temp(:,3)+temp(LE,3),'xb')
%                axis equal;
                R=P+temp(LE,2:4);
                S=Q+temp(TE,2:4);
                RS = R-S;
%                quiver(Q(1),Q(2),RS(1),RS(2));
%                data(iSpan, 5) =  180/pi* acos(dot(RS(1:2:3),PQ(1:2:3))/(vecnorm(RS(1:2:3))*vecnorm(PQ(1:2:3))));
%                data(iSpan, 6) =  180/pi* acos(dot(RS(1:2),PQ(1:2))/(vecnorm(RS(1:2))*vecnorm(PQ(1:2))));

                index = 1:2:3;
                a = RS(index(1))*PQ(index(1));
                b = RS(index(2))*PQ(index(2));
                c = sqrt(PQ(index(1))^2+PQ(index(2))^2);
                d = sqrt(RS(index(1))^2+RS(index(2))^2);
                data(iSpan, 5) =  180/pi* acos((a+b)/(c*d));

                index = 1:2;
                a = RS(index(1))*PQ(index(1));
                b = RS(index(2))*PQ(index(2));
                c = sqrt(PQ(index(1))^2+PQ(index(2))^2);
                d = sqrt(RS(index(1))^2+RS(index(2))^2);
                
                arg=(a+b)/(c*d);
                if  arg>1
                    if round(arg,6) ==1
                      data(iSpan, 6) =  180/pi* acos(round(arg,6));
                    else
                      data(iSpan, 6) =  180/pi* acos(arg); %Imaginary results
                    end
                end

                T = temp(LP,2:4); %x,y,z coordinates for point T on suction side
                U = temp(HP,2:4); %x,y,z coordinates for point U on pressure side
                TU=T-U;

                V=T+temp(LP,2:4);
                W=U+temp(HP,2:4);
                VW = V-W;

                index = [2,3];
                a = VW(index(1))*TU(index(1));
                b = VW(index(2))*TU(index(2));
                c = sqrt(TU(index(1))^2+TU(index(2))^2);
                d = sqrt(VW(index(1))^2+VW(index(2))^2);
                data(iSpan, 4) =  180/pi* acos((a+b)/(c*d));
%               title(['ispan:' int2str(iSpan) ' theta:' num2str(data(iSpan, 6))])
            end
          
            for jj=1:6
                designvar.deflection{iLoad}=[designvar.deflection{iLoad} data(:,jj)];
            end 
        end
    %% ************************************************************************
% ================= READ STRESS RESULTANTS INTO MATLAB =================   
        if isfield(config.ansys.analysisFlags, 'resultantVSspan') && config.ansys.analysisFlags.resultantVSspan~=0
            fileName='resultantVSspan.txt';
            designvar.resultantVSspan{iLoad}=txt2mat(fileName);
            delete(fileName);
%             fileName='resultantVSspan2.txt';
%             designvar.resultantVSspan2{iLoad}=txt2mat(fileName);
%             delete(fileName);
          
        end


        
        %% ************************************************************************
        % ================= READ LINEAR BUCKLING RESULTS =================
        % read buckling results
        if isfield(config.ansys.analysisFlags,'globalBuckling') && config.ansys.analysisFlags.globalBuckling >0
            fid=fopen([bucklingFilename '.out']);
            for jj=1:5
                tline = fgetl(fid);
            end
            data=cell(1,5);
            while 1
                tline = fgetl(fid);
                if ~ischar(tline), break, end
                data=[data; textscan(tline,'%f %f %f %f %f')];
            end
            fclose(fid);
            disp(' ')
            data=cell2mat(data);
            linearLoadFactors=data(1:config.ansys.analysisFlags.globalBuckling,2); %Extract the load factors (LF) from the linear buckling analysis
            delete([bucklingFilename '.out'])
        end

        %% ************************************************************************
        % ================= PERFORM NON-LINEAR BUCKLING/WRINKLING ANALYSIS =================
        % Perform nonlinear buckling here if required (and writeANSYSgetFaceStresses 
        % at the end of the nonlinear analysis for wrikling check 
        if isfield(config.ansys.analysisFlags,'imperfection') && ~isempty(config.ansys.analysisFlags.imperfection)
            warning('output designvar. Currently does not work for nonlinear cases')
            imperfection=config.ansys.analysisFlags.imperfection./1000; %convert mm to m. 

            nonlinearLoadFactors=zeros(length(linearLoadFactors),length(imperfection)); 
            critDesignvar=zeros(length(imperfection),1);
            wrinklingLimitingElementData=zeros(length(linearLoadFactors),4,length(imperfection));
            marker={'-ok','-sk','-dk','-*k','-^k','-<k','->k','-pk','-hk'}; %Plot markers
            %SF=max(LLF); %Use one loads file for all buckling modes

            for jj=1:length(imperfection)
                for ii=1:length(linearLoadFactors) 
                   % For each load factor, create a new jobname and database and run a nonlinear static analysis 
                   nonlinearLoadFactors(ii,jj)=nonlinearBuckling(ansysFilename,ansysPath,ansys_product,config,ii,jj);
                   [wrinklingLimitingElementData(ii,:,jj)]=wrinklingForNonlinearBuckling(blade,config.ansys.analysisFlags.localBuckling,settings,ncpus,ansysFilename,ii,jj);
                end
                [minnLLF,minnLLFMode]=min(nonlinearLoadFactors(:,jj))
                [minWLF,minWLFMode]=min(wrinklingLimitingElementData(:,3,jj))
                critDesignvar(jj)=min(minnLLF,minWLF)
            end

            figure(5)
            for k=1:length(linearLoadFactors)
                %disp(strcat('-',marker(j),'k'))
                plot(imperfection*1000,nonlinearLoadFactors(k,:),marker{k})
                hold on;
            end
            legend(strcat('Mode-',num2str([1:length(linearLoadFactors)]')))
            title('Imperfection Study (Linear Elements) SNL3p0-148-mk0p2-s1-fiberglass')
            xlabel('Max Imperfection Size [mm]')
            ylabel('Buckling Load Factors [ ]')


            %wrinklingLimitingElementData - [ansysSecNumber elno lf phicr]
            designvar.globalBuckling{iLoad}= min(min(critDesignvar));
        elseif isfield(config.ansys.analysisFlags,'globalBuckling') && config.ansys.analysisFlags.globalBuckling >0
            designvar.globalBuckling{iLoad}=linearLoadFactors(1);
        end


        %% ************************************************************************
        % ================= POST-PROCESS PANEL WRINKLING FACTORS =================
        if isfield(config.ansys.analysisFlags,'localBuckling') && ~isempty(config.ansys.analysisFlags.localBuckling)

                    
            if isfield(config.ansys.analysisFlags,'imperfection') && ~isempty(config.ansys.analysisFlags.imperfection)
                
                %UNSUPPORTED AT THIS TIME 
%                 filename=strcat(ansysFilename,'-',int2str(ii),'-',int2str(jj)); %The name of the next job name
%                 %%%%%%%        Generate Wrinkling Files        %%%%%%%%%%%%%
%                 script_name=strcat('commands4-',int2str(ii),'.mac');
%                 script_out=strcat('output4-',int2str(ii),'-',int2str(jj),'.txt');
% 
%                 fid=fopen(script_name,'w+');
%                 fprintf(fid,strcat('!************   MODE-%i   ************\n'),ii);
%                 fprintf(fid,'/FILNAME,''%s'',1\n',filename);   %From master, change the jobname
%                 fprintf(fid,'resume\n');
%                 fprintf(fid,'/POST1\n');
%                 fprintf(fid,'SET,LAST\n'); 
% 
%                 [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,config.ansys.analysisFlags.localBuckling);
% 
%                 fprintf(fid,'/EXIT,NOSAVE\n');
%                 fclose(fid);
% 
%                 ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',settings.ansys_path,settings.ansys_product,script_name,script_out,int2str(np))    % KMP_STACKSIZE is 512k by default. This is not enough therefore SET
%                 % KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
%                 % simulations. EC
%                 % 
% 
% 
%                 system(ansys_call)  % the windows system call to run the above ansys command
%                 fprintf('%s: Nonlinear Mode-%s Analysis Finished\n',datestr(now),int2str(ii))
                
            end
            % perform wrinkling check
            [wrinklingLimitingElementData]=Fagerber2005wricklingCheck(app,SkinAreas,compsInModel,config.ansys.analysisFlags.localBuckling);
            designvar.localBuckling{iLoad}=wrinklingLimitingElementData(3);
            delete *faceAvgStresses.txt
        end




    %% ************************************************************************
% ================= READ FAILURE RESULTS INTO MATLAB =================   

        if isfield(config.ansys.analysisFlags,'failure') && ~isempty(config.ansys.analysisFlags.failure)
            fileName=[failureFilename '.out'];
            designvar.failure{iLoad}=read_1_ANSYSoutput(fileName);
            delete(fileName)
        end
    end
    %% ************************************************************************
    
% ================= RUN FATIGUE POST PROCESSOR ================= 
    %After all load directions are solved compute fatige damage if needed
    if isfield(config.ansys.analysisFlags,'fatigue') && ~isempty(config.ansys.analysisFlags.fatigue)
        cd ..
        [wt,rccdata]=getWindSpeedDistribution(IEC.avgws);
        cd 'NuMAD'
        designvar.fatigue=layupDesign_ANSYSfatigue(ansysBladeMaterials,wt,rccdata,IEC,loadsTable,config);
    end
end


%% Non-linear Buckling Analysis Script
function nonlinearLoadFactors=nonlinearBuckling(ansysFilename,ansysPath,ansys_product,config,ii,jj)
    warning('output designvar. Currently does not work for nonlinear cases')


    script_name=strcat('commands3-',int2str(ii),'.mac');
    script_out=strcat('output3-',int2str(ii),'-',int2str(jj),'.txt');


    fid=fopen(script_name,'w+');
    fprintf(fid,strcat('!************   MODE-%i   ************\n'),ii);
    fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);   %From master, change the jobname
    fprintf(fid,'resume\n');

    %Get Max displacement, UY
    fprintf(fid,'/POST1\n');
    fprintf(fid,'SET,1,%i\n',ii); %Read in results
    fprintf(fid,'nsel, all, node\n'); %Select all nodes
    fprintf(fid,'*get, Zncount, node,0,count\n'); %Find the number (quantity) of nodes selected
    fprintf(fid,'*dim,zNodeDisp,array,Zncount,1 \n'); %Allocate memory for an arry to hold nodal disp.

    %For each node populate array with Y-displacements
    fprintf(fid,'*DO, i,1,Zncount\n');    
    fprintf(fid,'*VGET, zNodeDisp(i,1), NODE, i, U, Y\n');
    fprintf(fid,'*ENDDO\n');

    %Find the min/max disp. value
    fprintf(fid,'*VSCFUN,zMaxUY,max,zNodeDisp\n'); 
    fprintf(fid,'*VSCFUN,zMinUY,min,zNodeDisp\n');

    fprintf(fid,'zImperfectionSF=%f/max(abs(zMaxUY),abs(zMinUY))\n',config.ansys.analysisFlags.imperfection(jj));

    %fprintf(fid,'zImperfectionSF=%f\n',config.ansys.analysisFlags.imperfection(jj));

    fprintf(fid,'/prep7\n');   
    fprintf(fid,'UPGEOM,zImperfectionSF,1,%i,''%s'',''rst''\n',ii,ansysFilename); %Introduce geometric imperfection from buckled mode shape-i                     

    fprintf(fid,'FINISH\n');  %Finish command requiered so that UPGEOM works after being in /PREP7
    filename=strcat(ansysFilename,'-',int2str(ii),'-',int2str(jj)); %The name of the next job name
    fprintf(fid,'/FILNAME,''%s'',1\n',filename);   %From master, change the jobname to master-1, master-2, etc...
    %fprintf(fid,strcat('SAVE,''',filename,''',''db'',''',strrep(pwd,'\','\\'),'''\n'));

    fprintf(fid,'\n');
    %Nonlinear Static Analysis
    fprintf(fid,'/solu\n');
    fprintf(fid,'antype,0\n');
    %fprintf(fid,'irlf,-1\n');
    fprintf(fid,'pstres,0\n');
    fprintf(fid,'NLGEOM,1\n');
    fprintf(fid,'TIME,%f\n',1);
    loadScaleFactor=5; %Make sure to scale the load arbitrarily high enought such that the solution never converges
    fprintf(fid,'allsel\n');
    fprintf(fid,'esel,s,type,,33\n'); %Select all follower elements
    fprintf(fid,'SFSCALE,PRES,%f\n',loadScaleFactor);
    fprintf(fid,'OUTRES,NSOL,ALL\n');
    fprintf(fid,'NROPT,UNSYM\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'solve\n');
    fprintf(fid,'FINISH\n');

    fprintf(fid,strcat('SAVE,''',filename,''',''db'',''',strrep(pwd,'\','\\'),'''\n'));
    fprintf(fid,'/EXIT,NOSAVE\n');

    fclose(fid);

    %%%%%% RUN ANSYS%%%%%%%
    ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',ansysPath,ansys_product,script_name,script_out,int2str(ncpus))    % KMP_STACKSIZE is 512k by default. This is not enough therefore SET
    % KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
    % simulations. EC
    % 

    
    system(ansys_call)  % the windows system call to run the above ansys command
    fprintf('%s: Nonlinear Mode-%s Analysis Finished\n',datestr(now),int2str(ii))
    

    data = readANSYSoutputs(strcat(filename,'.mntr'),10);
    a=size(data);
    nonlinearLoadFactors=data(a(1),7)*loadScaleFactor; %Extract the load level on the laster iteration before nonconvergence
end

%% Faceplate Wrinkling Post-Processing Script
function [wrinklingLimitingElementData]=wrinklingForNonlinearBuckling(blade,coreMatName,settings,np,ansysFilename,i,j)
        filename=strcat(ansysFilename,'-',int2str(i),'-',int2str(j)); %The name of the next job name
        %%%%%%%        Generate Wrinkling Files        %%%%%%%%%%%%%
        script_name=strcat('commands4-',int2str(i),'.mac');
        script_out=strcat('output4-',int2str(i),'-',int2str(j),'.txt');
        
        fid=fopen(script_name,'w+');
        fprintf(fid,strcat('!************   MODE-%i   ************\n'),i);
        fprintf(fid,'/FILNAME,''%s'',1\n',filename);   %From master, change the jobname
        fprintf(fid,'resume\n');
        fprintf(fid,'/POST1\n');
        fprintf(fid,'SET,LAST\n'); 
        
        [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,coreMatName);
        
        fprintf(fid,'/EXIT,NOSAVE\n');
        fclose(fid);
        
        ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',settings.ansys_path,settings.ansys_product,script_name,script_out,int2str(np))    % KMP_STACKSIZE is 512k by default. This is not enough therefore SET
        % KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
        % simulations. EC
        % 
        
        
        system(ansys_call)  % the windows system call to run the above ansys command
        fprintf('%s: Nonlinear Mode-%s Analysis Finished\n',datestr(now),int2str(i))
        

        
        [wrinklingLimitingElementData]=Fagerber2005wricklingCheck(app,SkinAreas,compsInModel,coreMatName);
end

%% Panel Stresses Analysis Script
function [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,coreMatName)
  
    [isoorthoInModel,compsInModel,SkinAreas,app] = getMatrialLayerInfoWithOutGUI(blade);
    %fid=fopen('getFaceStresses.mac','w+');
    TotalStations=numel(blade.ispan);
    for kStation = 1:TotalStations-1 %Loop along span
        %kPanel=find(~cellfun('isempty',strfind([SkinAreas(kStation).Material],'PANEL'))); %Array that stores the kArea index that contains 'PANEL' in the name
        %for i=1:numel(kPanel)
        for kArea = 1:numel(SkinAreas(kStation).startIB)

            %See if the section contatins Balsa/core material name (i.e find
            %the sandwhich panels)
            n = strcmp(SkinAreas(kStation).Material{kArea},app.matlist);
            mat = app.matdb(n);
            if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
                ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);

                writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)

            end
        end
    end
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    fprintf(fid,'!*************** WEB ***************\n');
    fprintf(fid,'\n');
    fprintf(fid,'\n');
    %fclose(fid);
    %Web
    TotalShearwebs = numel(app.shearweb);
    for kShearweb = 1:TotalShearwebs
        n = strcmp(app.shearweb(kShearweb).Material,app.matlist);
        mat = app.matdb(n);
        if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
            ansysSecNumber = find(strcmp({app.shearweb(kShearweb).Material},compsInModel)==1);
            ansysSecNumber=ansysSecNumber+1000;
            writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)

        end
    end
    fprintf(fid,'FINISH\n');
    fprintf(fid,'allsel\n');

end

function writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)
    %%%%%Find the face sheet%%%%
    cellMat={};
    for i=1:length(mat.layer)
        cellMat = [cellMat; {mat.layer(i).layerName}];   %Create a cell array to use "find"
    end
    kbalsa=find(strcmp(coreMatName,cellMat)==1);
    iLayer=1:(kbalsa-1); %Number of distinct materials in the top face

    % Find the number of layers in the face
    qty=0; %Quantity of layers counter
    for i=1:numel(iLayer)
        qty=qty+mat.layer(iLayer(i)).quantity;
    end

    %Loop through the top facesheet layers

    fprintf(fid, '!*************** ansysSecNumber = %i ***************\n',ansysSecNumber);
    fprintf(fid, '/POST1\n');
    fprintf(fid, '*DEL,iel\n');
    fprintf(fid, '*DEL,enum\n');
    fprintf(fid, '*DEL,nelTemp\n');
    fprintf(fid, 'RSYS, SOLU\n');
    fprintf(fid, 'ALLSEL\n');
    fprintf(fid, 'ESEL, S, SEC,,%i\n',ansysSecNumber);
    fprintf(fid, '*GET, enum, ELEM, 0, NUM, MIN, !  lowest element number in the selected set\n');
    fprintf(fid, '*get, nelTemp, ELEM,0,count\n');
    fprintf(fid, '*DIM, iel,ARRAY,nelTemp\n');
    fname=strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses'); %Text file to store averages stresses to be computed by ansys

    if isfile(strcat(fname,'.txt'))
        delete(strcat(fname,'.txt')) %Check if file exisits and delete it because following ansys commands append to file.
    end
    fprintf(fid, '*CFOPEN, %s, txt,,APPEND\n',fname);

    %Create an array with the element numbers in the selected set
    fprintf(fid, '*DO, J, 1,nelTemp  !Loop through elements\n');
    fprintf(fid, 'iel(J)=enum\n');
    fprintf(fid, 'enum =ELNEXT(enum)  !Next higher element number above N in selected set\n');
    fprintf(fid, '*ENDDO\n');

    fprintf(fid, '\n');
    fprintf(fid, 'ALLSEL\n');
    fprintf(fid, '*DO, J, 1,nelTemp  !Loop through elements\n');
    fprintf(fid, '	S11a=0 !Initialize average stress variables for each element\n');
    fprintf(fid, '	S22a=0\n');
    fprintf(fid, '	S33a=0\n');
    fprintf(fid, '	S23a=0\n');
    fprintf(fid, '	S13a=0\n');
    fprintf(fid, '	S12a=0\n');
    fprintf(fid, '	*DO, I, 1,%i    !Loop through face layers\n',qty);
    fprintf(fid, '	    LAYER,I\n');
    fprintf(fid, '		SHELL,MID   !Stress result at midlayer\n');
    fprintf(fid, '		ESEL,S,ELEM,,iel(J)\n');
    fprintf(fid, '		ETABLE,ERAS !Each element gets a new element table\n');
    fprintf(fid, '	    ETABLE,S11,S,X,AVG !AVG - Store averaged element centroid value\n');
    fprintf(fid, '	    ETABLE,S22,S,Y,AVG\n');
    fprintf(fid, '		ETABLE,S33,S,Z,AVG\n');
    fprintf(fid, '		ETABLE,S23,S,YZ,AVG\n');
    fprintf(fid, '		ETABLE,S13,S,XZ,AVG\n');
    fprintf(fid, '        ETABLE,S12,S,XY,AVG\n');
    fprintf(fid, '	    *GET,tempS11, ELEM, iel(J), ETAB, S11\n');
    fprintf(fid, '		*GET,tempS22, ELEM, iel(J), ETAB, S22\n');
    fprintf(fid, '		*GET,tempS33, ELEM, iel(J), ETAB, S33\n');
    fprintf(fid, '		*GET,tempS23, ELEM, iel(J), ETAB, S23\n');
    fprintf(fid, '		*GET,tempS13, ELEM, iel(J), ETAB, S13\n');
    fprintf(fid, '		*GET,tempS12, ELEM, iel(J), ETAB, S12\n');
    fprintf(fid, '		S11a=S11a+tempS11\n');
    fprintf(fid, '		S22a=S22a+tempS22\n');
    fprintf(fid, '		S33a=S33a+tempS33\n');
    fprintf(fid, '		S23a=S23a+tempS23\n');
    fprintf(fid, '		S13a=S13a+tempS13\n');
    fprintf(fid, '		S12a=S12a+tempS12\n');
    fprintf(fid, '	*ENDDO\n');
    fprintf(fid, '	S11a=S11a/%i\n',qty);
    fprintf(fid, '	S22a=S22a/%i\n',qty);
    fprintf(fid, '	S33a=S33a/%i\n',qty);
    fprintf(fid, '	S23a=S23a/%i\n',qty);
    fprintf(fid, '	S13a=S13a/%i\n',qty);
    fprintf(fid, '	S12a=S12a/%i\n',qty);
    fprintf(fid, '	ELNO=iel(J)  !It is needed to refer to ELNO in the command below\n');
    fprintf(fid, '*VWRITE,ELNO,S11a,S22a,S33a,S23a,S13a,S12a\n');
    fprintf(fid, '(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
    fprintf(fid, '*ENDDO\n');
    fprintf(fid, '*CFCLOS\n');
    fprintf(fid, '\n');
    fprintf(fid, '\n');
end

%% Wrinkling Analysis Post-Processing Script
function [limitingElementData]=Fagerber2005wricklingCheck(app,SkinAreas,compsInModel,coreMatName)
    %limitingElementData - [ansysSecNumber elno lf phicr]
    TotalStations = numel(app.station);
    TotalShearwebs = numel(app.shearweb);

    %%%%%%%%%%%%%%%%%   Main loop #1: loop around aero shell.   %%%%%%%%%%%%%%%%%
    LF=[];
    for kStation = 1:TotalStations-1 %Loop along span
        for kArea = 1:numel(SkinAreas(kStation).startIB)     
            %See if the section contatins Balsa/core material name (i.e find
            %the sandwhich panels)
            n = strcmp(SkinAreas(kStation).Material{kArea},app.matlist);
            mat = app.matdb(n);
            if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
                ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
                file=strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses.txt');
                avgFaceStress=txt2mat(file);  
                delete(file)
                LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName);
            end
        end
    end

     
    %%%%%%%%%%%%%%%%%   Main loop #2: loop along web.   %%%%%%%%%%%%%%%%%
    for kShearweb = 1:TotalShearwebs
        n = strcmp(app.shearweb(kShearweb).Material,app.matlist);
        mat = app.matdb(n);
        if contains([mat.layer.layerName],coreMatName)  %%%%%%%%%%%%%%%%%This logic will need to be made more general
            ansysSecNumber = find(strcmp({app.shearweb(kShearweb).Material},compsInModel)==1);
            file=strcat('section-',int2str(ansysSecNumber+1000),'-faceAvgStresses.txt');
            avgFaceStress=txt2mat(file);
            delete(file)
            LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber+1000,avgFaceStress,app,mat,coreMatName);
        end
    end

    [minLF,index]=min(LF(:,3));
    
    limitingElementData=LF(index,:);
    fprintf('\n\n The minimum wrinkling LF is: %f, wrinkle angle: %.2f°' ,minLF, LF(index,4))
    fprintf('\n and occurs in section number %i, element number %i\n, ',LF(index,1),LF(index,2))

    % [maxLF,index]=max(LF(:,3));
    % fprintf('\n\n The maximum LF is: %f, wrinkle  angle: %.2f°' ,maxLF, LF(index,4))
    % fprintf('\n and occurs in section number %i, element number %i, ',LF(index,1),LF(index,2))

end

function LF=getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName)
    %This is a recursive function for LF. It appends to the list of LF for each
    %element that has a positive LF.  EC
    [m,~]=size(avgFaceStress);
    for i=1:m %For each element in the current section number (ansysSecNumber)
        elno=avgFaceStress(i,1); %Element number
        S11a=avgFaceStress(i,2); %N/m^2 
        S22a=avgFaceStress(i,3); %N/m^2
        S12a=avgFaceStress(i,7); %N/m^2
        %Ignoring other stresses for the time being.
        %if elno==305
            %disp('press pause')
            %pause(10)

            [lf,  phicr]=checkWrinkle([S11a;S22a;S12a],mat,app,coreMatName);
            if lf>=0
               LF=[LF;ansysSecNumber elno lf phicr]; %Append load factor to list of load factors
            end
        %end
    end



    function [lf, phicr]=checkWrinkle(S_alphaBeta,mat,app,coreMatName)
        % For a single finite element, given the average in-plane stresses
        % in a face-sheet of that element, compute the load factor for that
        % element. EC

        % lf    - scalar load factor for the element
        % phicr - an angle, degrees. The direction of wrinkling
        % S11a,S22a,S12a - respective average face sheet stress
        % mat - material object
        % app - blade data


        %Locate the face sheet
        cellMat={};
        for i=1:length(mat.layer) 
            cellMat = [cellMat; {mat.layer(i).layerName}];   %Create a cell array to use "find"    
        end

        kbalsa=find(strcmp(coreMatName,cellMat)==1);
        iLayer=1:(kbalsa-1); %Index of distinct materials in the top face  
        %ilayer=(kbalsa+1):numel(cellMat)); %Number of distinct materials in the bottom face
        matCore=app.matdb(find(strcmp(coreMatName,app.matlist)==1));

        if strcmp(matCore.type,'orthotropic')   
            Ec=matCore.ez; %Note that matCore.ez is invarient for in-plane rotations only 
        elseif strcmp(matCore.type,'isotropic')
            Ec=matCore.ex;
            Gc=Ec/(2*(1+matCore.nuxy));
        else
            errordlg(sprintf('Material  "%s"  not found in database.', matCore.type),'Error');
            error('Material type "%s" not found in database.', matCore.type);
        end
        %ilayer=(kbalsa+1):numel(cellMat)); %Number of distinct materials in the bottom face
        dangle=2; %Resolution for angle (degrees)
        N=180/dangle+1; %Number of angles to evaluate
        angle = 0; 

        invLF=zeros(N,1);
        %Apt=zeros(3,3);
        %Bpt=zeros(3,3);
        Dpt=zeros(3,3);

        %Find total height of facesheet
        h=0;
        for klay = 1:numel(iLayer)
            h = h+mat.layer(iLayer(klay)).thicknessA * mat.layer(iLayer(klay)).quantity;
        end

        for kang = 1:N
           if strcmp(matCore.type,'orthotropic')   
                Gc=1/(sind(angle)^2*(1/matCore.gyz)+cosd(angle)^2*(1/matCore.gxz)); %In-plate clockwise rotation of: angle 
           end
            %Asssuming all ply angles are zero
            R_sig=[cosd(angle)^2, sind(angle)^2, -2*sind(angle)*cosd(angle);  %In-plate clockwise rotation of: angle 
                   sind(angle)^2, cosd(angle)^2, 2*sind(angle)*cosd(angle)
                   sind(angle)*cosd(angle), -sind(angle)*cosd(angle), cosd(angle)^2-sind(angle)^2];

            z1=-h/2; %The coordinate location of the bottom of the face (with the origin in the midplane of the face sheet
            for klay = 1:numel(iLayer)
                z2=z1+mat.layer(iLayer(klay)).thicknessA * mat.layer(iLayer(klay)).quantity;

                matklay = app.matdb(find(strcmp(mat.layer(iLayer(klay)).layerName,app.matlist)==1));
                % Bulid Plane Stress reduced compliance matrix for each
                % layer
                %fprintf('z1 = %f z2 = %f mat = %f   %s\n',z1,z2,matListnumber,mat.layer(klay).layerName)

                %Entries common to either isotropic or orthotropic entries
                Se=zeros(3);
                Se(1,1)=1/matklay.ex;
                Se(1,3)=0; %Valid for orthotropic materials only
                Se(2,3)=0; %Valid for orthotropic materials only
                Se(3,1)=0; %Valid for orthotropic materials only
                Se(3,2)=0; %Valid for orthotropic materials only

                if strcmp(matklay.type,'orthotropic')   
                    Se(1,2)=-matklay.prxy/matklay.ex;
                    Se(2,1)=-matklay.prxy/matklay.ex;
                    Se(2,2)=1/matklay.ey;
                    Se(3,3)=1/matklay.gxy;
                elseif strcmp(matklay.type,'isotropic')
                    Se(1,2)=-matklay.nuxy/matklay.ex;
                    Se(2,1)=-matklay.nuxy/matklay.ex;
                    Se(2,2)=1/matklay.ex;
                    Se(3,3)=2*(1+matklay.nuxy)/matklay.ex;
                else
                    errordlg(sprintf('Material  "%s"  not found in database.', matklay.type),'Error');
                    error('Material type "%s" not found in database.', matklay.type);
                end

                %Apt=Apt+R_sig*inv(Se)*R_sig'*(z2-z1);

                %Bpt=Bpt+1/2*R_sig*inv(Se)*R_sig'*(z2^2-z1^2);
                Dpt=Dpt+1/3*R_sig*inv(Se)*R_sig'*(z2^3-z1^3);
                z1=z2;
            end 
            Pcr=-3/2*(2*Dpt(1,1)*Ec*Gc)^(1/3); %N/m
            Pphi=(R_sig(1,:)*S_alphaBeta)*h; %N/m

            invLF(kang)=Pphi/Pcr;
            angle=angle+dangle;
        end
        [invlf, phicr_i]=max(invLF);
        lf=1/invlf;

        phicr=(phicr_i-1)*dangle;

%         if lf>1e6
%             figure(2)
%             plot(angle, Pcr,'k')
%             xlabel('Angle, \phi [deg]')
%             ylabel('Load [N/m]')
%             hold on; 
% 
%             plot(angle, Pphi,'r')
%             plot(angle, lf*Pphi,'b')
%             legend('P_c_r','P_\phi',strcat(num2str(lf),'P_\phi'))
%             hold off;
%             figure(3)
%             plot(angle,invLF,'r')
%             hold on; 
%             plot(phicr,invlf,'*')
%             xlabel('Angle, \phi [deg]')
%             ylabel('1/\lambda [ ]')
%             hold off;
%             fprintf('\n %8.2g  %i | %8.2g %8.2g | %8.2g %8.2g %8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
%             pause(10)
%         else
%             %fprintf('\n %8.2g  %i   %8.2g %8.2g   %8.2g %8.2g %8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
%         end
    end
end

function designvar=saveData(designvar,iLoad,airfoilSegmentName,iSpan,nodes,midNodei)
    designvar.localFields{iLoad}.(airfoilSegmentName).x(iSpan) =nodes(midNodei,2);
    designvar.localFields{iLoad}.(airfoilSegmentName).y(iSpan) =nodes(midNodei,3);
    designvar.localFields{iLoad}.(airfoilSegmentName).z(iSpan) =nodes(midNodei,4);
    designvar.localFields{iLoad}.(airfoilSegmentName).data(iSpan)=nodes(midNodei,5);
end
function midNodei=findCenterNode(nodes,direction)        
    xlist=nodes(:,direction+1); %list of nodal x values
    [xmax,~]=max(xlist);
    [xmin,~]=min(xlist);

    xmid=mean([xmin,xmax]); %x coordinate of sparcap centerline

    nnode=numel(nodes(:,1));
    minDiff=100; %initialize
    for iNode=1:nnode
        x=nodes(iNode,direction+1);    %x position of node
        diff = abs(x-xmid); %distance between the node and sparcap centerline
        if diff < minDiff
          midNodei=iNode;
          minDiff =diff;
        end

    end
end