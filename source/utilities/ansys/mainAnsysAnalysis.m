function [designvar] = mainAnsysAnalysis(blade,meshData,loadsTable,analysisConfig,varargin)
    anFlagNames = fieldnames(analysisConfig.analysisFlags);
    
    global ansysPath
    ansys_product = 'ANSYS';

    if isfield(analysisConfig.analysisFlags,'imperfection') && ~isempty(analysisConfig.analysisFlags.imperfection) &&analysisConfig.analysisFlags.globalBuckling ==0
        error('Specify number of buckling modes when performing nonlinear buckling')
    end

    % Original mesh file to analize
    if isfield(analysisConfig,'meshFile') && ~isempty(analysisConfig.meshFile)
        %Do nothing
    else
        analysisConfig.meshFile='master';
    end
    
    % File name base name for ansys analysis files
    if isfield(analysisConfig,'analysisFileName') && ~isempty(analysisConfig.analysisFileName)
        ansysFilename = [analysisConfig.analysisFileName];
    else
        ansysFilename=['FEmodel'];
    end
    % Number of CPUs to use
    if isfield(analysisConfig,'np') && ~isempty(analysisConfig.np)
        if analysisConfig.np <1
            error('analysisConfig.np must be greater than zero')
        else
            ncpus=analysisConfig.np;
        end
    else
        ncpus=1; %Default value
    end
    %Initialize
    for i=1:length(anFlagNames)
        
        if strcmpi('globalBuckling',anFlagNames{i}) || strcmpi('resultantVSspan',anFlagNames{i}) || strcmpi('deflection',anFlagNames{i}) || strcmpi('mass',anFlagNames{i});
            if analysisConfig.analysisFlags.(anFlagNames{i})~=0 ;
                designvar.(anFlagNames{i}) = cell(1,length(loadsTable));
            end
        else strcmpi('localBuckling',anFlagNames{i}) || strcmpi('failure',anFlagNames{i}) || strcmpi('fatigue',anFlagNames{i}) ||strcmpi('imperfection',anFlagNames{i}) || strcmpi('mass',anFlagNames{i});
            if ~isempty(analysisConfig.analysisFlags.(anFlagNames{i}));
                designvar.(anFlagNames{i}) = cell(1,length(loadsTable));
            end
        end
    end
    
    if ~exist('designvar')
        error('no analyses are configured in configuration st.')
    end
    
    for iLoad=1:length(loadsTable)
        %% ************************************************************************
        % ================= APPLY LOADS TO FEA MESH =================
        forcefilename='forces';
        nodeData=[meshData.outerShellNodes meshData.nodes(meshData.outerShellNodes, :)];
        if isfield(analysisConfig.analysisFlags,'FollowerForces') && ~isempty(analysisConfig.analysisFlags.FollowerForces) &&analysisConfig.analysisFlags.FollowerForces~=0 && isfield(analysisConfig.analysisFlags, 'StaticNonlinear') && ~isempty(analysisConfig.analysisFlags.StaticNonlinear) &&analysisConfig.analysisFlags.StaticNonlinear~=0
            beamForceToAnsysShellFollower('map3D_fxM0',nodeData,loadsTable{iLoad},strcat(forcefilename,'.src'));
        else
            beamForceToAnsysShell('map3D_fxM0',nodeData,loadsTable{iLoad},strcat(forcefilename,'.src'));
        end

        disp('Forces mapped to ANSYS model')

        %% ************************************************************************
        % ================= PERFORM LINEAR STATIC ANALYSIS =================
        % run buckling computations in ansys
        disp(' '); disp('Running ANSYS analysis...')

        script_name='ansysAnalysis.mac';
        script_out='ansysAnalysisEcho.out';

        fid=fopen(script_name,'w+');

        fprintf(fid,'/NERR,,99999999\n');
        fprintf(fid,'/CWD, ''%s''\n',pwd);
        fprintf(fid,'resume,master,db\n');
%         fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);   %From master, change the jobname
        fprintf(fid,'/FILNAME,''%s'',1\n',strcat(ansysFilename,'-Load', int2str(iLoad)));   %From master, change the jobname
        %fprintf(fid,'resume\n');

        fprintf(fid,'! BEGIN LINEAR STATIC SCRIPT\n');
        fprintf(fid,'esel,all\n');
        fprintf(fid,'/prep7\n');
        fprintf(fid,'fdel,all\n');
        fprintf(fid,'/input,%s,src\n',forcefilename);

        %Linear Static Analysis
        fprintf(fid,'/solu\n');
        
        fprintf(fid,'antype,static\n');
        
        if isfield(analysisConfig.analysisFlags,'StaticNonlinear') && ~isempty(analysisConfig.analysisFlags.StaticNonlinear) &&analysisConfig.analysisFlags.StaticNonlinear~=0
            fprintf(fid,'nlgeom,1\n'); %%%%%%%%%%%%%%%%%%%%%%%% TEMP
            fprintf(fid,'OUTRES,all,ALL\n');%%%%%%%%%%%%%%%%%%%%%%%% TEMP
%         else
%             fprintf(fid,'pstres,on\n');
        end
        fprintf(fid,'irlf,-1\n');
        
        fprintf(fid,'bcsoption,,incore\n');
        
        fprintf(fid,'solve\n');
        fprintf(fid,'finish\n');
        
        %Only compute mass on the first load case
        if iLoad==1 && isfield(lower(analysisConfig.analysisFlags),'mass') && ~isempty(analysisConfig.analysisFlags.mass)&& analysisConfig.analysisFlags.mass~=0
            %Get Mass Here
            fprintf(fid,'*GET, Z_mass, ELEM, 0, MTOT, X\n');
            fprintf(fid,'/output, mass,txt\n');
            fprintf(fid,'*status,Z_mass\n');
            fprintf(fid,'/output\n');
            fprintf(fid,'finish\n');
        end
 
 %% ************************************************************************
%================= PERFORM Deflection ANALYSIS =================

        if isfield(analysisConfig.analysisFlags,'deflection') && analysisConfig.analysisFlags.deflection~=0
            
            deflectionFilename = 'results_deflection';
            writeAnsysDeflections(blade, analysisConfig, iLoad, fid, deflectionFilename)
            
        end
        % calculate face stresses for wrinkling
        if isfield(analysisConfig.analysisFlags,'localBuckling') && ~isempty(analysisConfig.analysisFlags.localBuckling)&& ~(isfield(analysisConfig.analysisFlags,'imperfection') && ~isempty(analysisConfig.analysisFlags.imperfection))
             %Check for wrinkling here in a linear analysis
            [app,SkinAreas,compsInModel]=writeAnsysGetFaceStresses(blade,fid,analysisConfig.analysisFlags.localBuckling);
        end

        %%% Output resultant force and moments to file
        if isfield(analysisConfig.analysisFlags, 'resultantVSspan') && analysisConfig.analysisFlags.resultantVSspan~=0
            
            writeAnsysResultantVSSpan(blade, analysisConfig, iLoad, fid)

        end
        %% ************************************************************************
        % ================= PERFORM FATIGUE ANALYSIS =================
        if isfield(analysisConfig.analysisFlags,'fatigue') && ~isempty(analysisConfig.analysisFlags.fatigue)
            
            writeAnsysFatigue(fid,iLoad)
            
        end
%% ************************************************************************
% ================= CREAT LOCAL FIELD RESULTS FOR MATLAB =================

        if isfield(analysisConfig.analysisFlags,'localFields') && ~isempty(analysisConfig.analysisFlags.localFields)
            
            writeAnsysLocalFields(blade, analysisConfig, iLoad, fid)
         
        end

        %% ************************************************************************
        % ================= PERFORM FAILURE ANALYSIS =================

        % Initialize GUI commands from batch operation to identify maxima
        if isfield(analysisConfig.analysisFlags,'failure') && ~isempty(analysisConfig.analysisFlags.failure)
            failureFilename = 'results_failure';
            
            writeAnsysRupture(analysisConfig, iLoad, fid, failureFilename)
            
        end


        %% ************************************************************************
        % ================= PERFORM BUCKLING ANALYSIS =================

        %Linear Buckling Analysis
        if isfield(analysisConfig.analysisFlags,'globalBuckling') && analysisConfig.analysisFlags.globalBuckling >0
            bucklingFilename = 'results_buckling';
            
            writeAnsysLinearBuckling(blade, analysisConfig, iLoad, fid, bucklingFilename)
            
        elseif isfield(analysisConfig.analysisFlags,'globalBuckling') && analysisConfig.analysisFlags.globalBuckling <0
            error('analysisConfig.analysisFlags.globalBuckling must be greater than or equal to zero')
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
        if iLoad==1 && isfield(lower(analysisConfig.analysisFlags),'mass') && ~isempty(analysisConfig.analysisFlags.mass)&& analysisConfig.analysisFlags.mass~=0
            designvar.mass=read_1_ANSYSoutput('mass.txt');
            delete mass.txt
        end

    %% ************************************************************************
% ================= READ DEFLECTION RESULTS INTO MATLAB =================       
        if isfield(analysisConfig.analysisFlags,'deflection') && analysisConfig.analysisFlags.deflection~=0

            designvar.deflection=readAnsysDeflections(blade, analysisConfig, iLoad, deflectionFilename);
            
        end
    %% ************************************************************************
% ================= READ STRESS RESULTANTS INTO MATLAB =================   
        if isfield(analysisConfig.analysisFlags, 'resultantVSspan') && analysisConfig.analysisFlags.resultantVSspan~=0
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
        if isfield(analysisConfig.analysisFlags,'globalBuckling') && analysisConfig.analysisFlags.globalBuckling >0
            
            linearLoadFactors = readAnsysLinearBuckling(blade, analysisConfig, iLoad, fid, bucklingFilename);
            
        end

        %% ************************************************************************
        % ================= PERFORM NON-LINEAR BUCKLING/WRINKLING ANALYSIS =================
        % Perform nonlinear buckling here if required (and writeANSYSgetFaceStresses 
        % at the end of the nonlinear analysis for wrikling check 
        if isfield(analysisConfig.analysisFlags,'imperfection') && ~isempty(analysisConfig.analysisFlags.imperfection)
            warning('output designvar. Currently does not work for nonlinear cases')
            imperfection=analysisConfig.analysisFlags.imperfection./1000; %convert mm to m. 

            nonlinearLoadFactors=zeros(length(linearLoadFactors),length(imperfection)); 
            critDesignvar=zeros(length(imperfection),1);
            wrinklingLimitingElementData=zeros(length(linearLoadFactors),4,length(imperfection));
            marker={'-ok','-sk','-dk','-*k','-^k','-<k','->k','-pk','-hk'}; %Plot markers
            %SF=max(LLF); %Use one loads file for all buckling modes

            for jj=1:length(imperfection)
                for ii=1:length(linearLoadFactors) 
                   % For each load factor, create a new jobname and database and run a nonlinear static analysis 
                   nonlinearLoadFactors(ii,jj)=writeAnsysNonLinearBuckling(ansysFilename,ansysPath,ansys_product,analysisConfig,ii,jj, ncpus, iLoad);
                   [wrinklingLimitingElementData(ii,:,jj)]=wrinklingForNonlinearBuckling(blade,analysisConfig.analysisFlags.localBuckling,settings,ncpus,ansysFilename,ii,jj);
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
        elseif isfield(analysisConfig.analysisFlags,'globalBuckling') && analysisConfig.analysisFlags.globalBuckling >0
            designvar.globalBuckling{iLoad}=linearLoadFactors(1);
        end

        %% ************************************************************************
        % ================= POST-PROCESS PANEL WRINKLING FACTORS =================
        if isfield(analysisConfig.analysisFlags,'localBuckling') && ~isempty(analysisConfig.analysisFlags.localBuckling)

                    
            if isfield(analysisConfig.analysisFlags,'imperfection') && ~isempty(analysisConfig.analysisFlags.imperfection)
                
                %UNSUPPORTED AT THIS TIME 
                writeAnsysNonLinearLocalBuckling(blade, analysisConfig, iLoad, fid, ansysFilename, ii, jj)
                
            end
            % perform wrinkling check
            [wrinklingLimitingElementData]=writeAnsysFagerberWrinkling(app,SkinAreas,compsInModel,analysisConfig.analysisFlags.localBuckling);
            designvar.localBuckling{iLoad}=wrinklingLimitingElementData(3);
            delete *faceAvgStresses.txt
        end

    %% ************************************************************************
% ================= READ FAILURE RESULTS INTO MATLAB =================   

        if isfield(analysisConfig.analysisFlags,'failure') && ~isempty(analysisConfig.analysisFlags.failure)
            fileName=[failureFilename '.out'];
            designvar.failure{iLoad}=read_1_ANSYSoutput(fileName);
            delete(fileName)
        end
    end
    %% ************************************************************************
    
% ================= RUN FATIGUE POST PROCESSOR ================= 
    %After all load directions are solved compute fatige damage if needed
    if isfield(analysisConfig.analysisFlags,'fatigue') && ~isempty(analysisConfig.analysisFlags.fatigue)  
        if ~isempty(varargin) && isequal(class(varargin{1}),'IECDef')
            cd ..
            IEC=varargin{1};
            [wt,rccdata]=getWindSpeedDistribution(IEC.avgws);
            cd 'NuMAD'
            designvar.fatigue=postprocessANSYSfatigue(blade,meshData,wt,rccdata,IEC,loadsTable,analysisConfig);
        else
            error('IECDef required to run fatigue analysis in mainAnsysAnalysis')
        end
    end
end


function designvar=saveData(designvar,iLoad,airfoilSegmentName,iSpan,nodes,midNodei)
    designvar.localFields{iLoad}.(airfoilSegmentName).x(iSpan) =nodes(midNodei,2);
    designvar.localFields{iLoad}.(airfoilSegmentName).y(iSpan) =nodes(midNodei,3);
    designvar.localFields{iLoad}.(airfoilSegmentName).z(iSpan) =nodes(midNodei,4);
    designvar.localFields{iLoad}.(airfoilSegmentName).data(iSpan)=nodes(midNodei,5);
end