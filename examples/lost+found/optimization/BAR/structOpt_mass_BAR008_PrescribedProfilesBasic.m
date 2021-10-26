function structOpt_mass_BAR008_PrescribedProfilesBasic(blade,lb,ub,A,b,pop,gen,restartFile,useParallel,compIndex,x0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************** INPUT OPTIONS ********************************
optRoutine = 'particle swarm';%'pattern search';%'genetic algorithm';%
 

if isempty(restartFile)
    useRestartFile = false;
else
    useRestartFile = true;
end
[initPop,initScores] = setInitPop(useRestartFile,restartFile);
% % if useRestartFile
% %     options = gaoptimset( 'PopulationSize',pop,...
% %         'Generations',gen,...
% %         'InitialPopulation',initPop,...
% %         'InitialScores',initScores,...
% %         'UseParallel',useParallel,...
% %         'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
% %     
% % % %     options = optimoptions(optRoutine,'PopulationSize',pop,'Generations',gen,...
% % % %         'InitialPopulation',initPop,'InitialScores',initScores,'UseParallel',useParallel,...
% % % %         'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
% % else
% %     options = gaoptimset('PopulationSize',pop,...
% %         'Generations',gen,...
% %         'UseParallel',useParallel,...
% %         'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
% %     
% % % %     options = optimoptions(optRoutine,'PopulationSize',pop,'Generations',gen,...
% % % %         'InitialPopulation',initPop,'InitialScores',initScores,'UseParallel',useParallel,...
% % % %         'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
% % end

nVar=length(lb);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***************** DELETE OUTPUT FROM PREVIOUS RUNS **********************
delete 'layupDesignProgress.txt'
delete 'layupDesignCandidates.txt'

delete('..\carray.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************* SET UP SOME OF THE INPUTS FOR THE ANALYSIS ****************
blade.establishPaths
% Define the blade mesh
blade.mesh = 0.45; %Debugging mesh size
%blade.mesh = 0.2;

% save the gage labels for the buckling and fatigue analyses (rotated coordinate system)

%% EMA: original:
% runIEC_ipt; % load in the params structure
%% changed to:
hm = pwd;
cd ..
runIEC_ipt;
cd(hm);
%% END

thetaMomentRotation = 0:params.momentMaxRotation:180; % [deg]
thetaMomentRotation(thetaMomentRotation==180)=[]; % remove repetitive 180 deg rotation
gage_labels={}; gage_maximum=[]; gage_theta=[];
for bb = 1 % blade number(s)
    for tt = 1:length(thetaMomentRotation)
        for gg = 0:9 % number of gages
            % save the gage variable names
            if gg==0
                % create the maximum moment vector
                gage_labels{bb}{2*(tt-1)+1}{gg+1} = ['MaxRes_RootMrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                gage_maximum{bb}{2*(tt-1)+1}(gg+1) = 1;
                gage_theta{bb}{2*(tt-1)+1}(gg+1) = thetaMomentRotation(tt);
                % create the minimum moment vector
                gage_labels{bb}{2*(tt-1)+2}{gg+1}  = ['MinRes_RootMrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                gage_maximum{bb}{2*(tt-1)+2}(gg+1) = -1;
                gage_theta{bb}{2*(tt-1)+2}(gg+1) = thetaMomentRotation(tt);
            else
                % create the maximum moment vector
                gage_labels{bb}{2*(tt-1)+1}{gg+1} = ['MaxRes_Spn' num2str(gg) 'MLrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                gage_maximum{bb}{2*(tt-1)+1}(gg+1) = 1;
                gage_theta{bb}{2*(tt-1)+1}(gg+1) = thetaMomentRotation(tt);
                % create the minimum moment vector
                gage_labels{bb}{2*(tt-1)+2}{gg+1} = ['MinRes_Spn' num2str(gg) 'MLrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                gage_maximum{bb}{2*(tt-1)+2}(gg+1) = -1;
                gage_theta{bb}{2*(tt-1)+2}(gg+1) = thetaMomentRotation(tt);
            end
        end
    end
end
% save the information into the fast_gage structure: fast_gage.labels{bb}{tt}
fast_gage.labels = gage_labels;
fast_gage.maximum = gage_maximum;
fast_gage.theta = gage_theta;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************** RUN OPTIMIZATION *****************************
if 1 % call the optimization algorithm
    % ================== Choose Optmization Algorithm =====================
    switch optRoutine
        case 'genetic algorithm'
            % set up options for genetic algorithm
            if useRestartFile
                options = optimoptions('ga','PopulationSize',pop,'Generations',gen,...
                    'InitialPopulation',initPop,'InitialScores',initScores,'UseParallel',useParallel,...
                    'PlotFcn',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@optCustomPlot} );
            else
                options = optimoptions('ga','PopulationSize',pop,'Generations',gen,...
                    'UseParallel',useParallel,...
                    'PlotFcn',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@optCustomPlot} );
            end
            % call the genetic algorithm
            [x,fval,exitflag,output,population,scores] = ga(@(x)designBladeFcn(x,blade,fast_gage,useParallel,false),nVar,A,b,[],[],lb,ub,[],[],options);
            save gaResults x fval exitflag output population scores
            
        case 'pattern search'
            % call the pattern search algorithm
            clear options
% %             options = optimoptions('patternsearch','MaxFunEvals',1000,'Cache','on','UseParallel',true);
            maxIter = 100*nVar;  % default
            maxFunEval = 2000*nVar; % default
            options = optimoptions('patternsearch','MaxIterations',maxIter,'MaxFunctionEvaluations',maxFunEval,...
                'Cache','on','UseParallel',true,'UseVectorized', false,'PlotFcn',@psplotbestf);
            
            xo = (lb+ub)/2;
            
%             options = optimoptions('patternsearch','UseParallel', true, 'UseCompletePoll', true, 'UseVectorized', false,'MaxFunctionEvaluations',1000);
            [x,fval,exitflag,output] = patternsearch(@(x)designBladeFcn(x,blade,fast_gage,useParallel,false),xo,A,b,[],[],lb,ub,options);
            save patternsResults x fval exitflag output
            
        case 'particle swarm'
            % call the particle swarm algorithm
            deltaf_stop = 5; % tolerance in best case to stop optimization [kg]
            
            if 0
                options = optimoptions('particleswarm','SwarmSize',64,...
                    'UseParallel',true,'PlotFcn',{@pswplotbestf,@optCustomPlot});
            else
                options = optimoptions('particleswarm','SwarmSize',pop,'FunctionTolerance',deltaf_stop,...
                    'UseParallel',true,'PlotFcn',@pswplotbestf);
            end
            disp(' ')
            [x,fval,exitflag,output] = particleswarm(@(x)designBladeFcn(x,x0,blade,fast_gage,useParallel,false,compIndex),nVar,lb,ub,options);
            save gaResults x fval exitflag output
    end
    
else % just do manual analysis and then quit
    x=[11.44 .093 30.61 12.51 279.95 .472];
    for i=1:size(x,1)
        f(i,1) = designBladeFcn(x(i,:), blade, false, false);
        save tmp x f
    end
    disp('results summary:')
    sprintf('%i %i %i %i %i %i %.1f\n',[x f]')
end
 
end
%%
function f = designBladeFcn(x,x0,blade,fast_gage,useParallel,saveBlade,compIndex)
% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Change folder saving structure so that multiple workers in a parallel
% algorithm don't overwrite each other.

% configure the ansys simulation
config.ansys.nonlinear = 0;
config.ansys.wrinkling = 0;
config.ansys.nmodes = 8;
if useParallel
    config.ansys.np = 1;
else
   config.ansys.np = 4;
end
config.ansys.imperfection_mm = 1;
config.ansys.meshFile = 'master.db';
config.ansys.bucklingFile = 'master-V3';
config.ansys.fc='SMAX'; % EMAX,SMAX,TWSI,TWSR,HFIB,HMAT,PUCK,LARC03,LARC04




hm = pwd; % reference to the 'NuMad' folder in the main directory
if useParallel
    t = getCurrentTask;
    if isempty(t)
        id = '0'; % first iteration of genetic algorithm in generation
    else
        id = num2str(t.ID); % remaining population are sent to workers in parallel
    end  
    % distinct folder for each worker (will be overwritten once that worker is
    % finished and called back - MUST save important information to main
    % directory if you want to keep it.
    caseFolder = ['parallel_worker_' id];
    mkdir(['..\' caseFolder '\NuMAD'])
    mkdir(['..\' caseFolder '\wind'])
    mkdir(['..\' caseFolder '\out'])
    % if exist(['..\' caseFolder],'dir'), rmdir(['..\' caseFolder],'s'), end
    copyfile('airfoils', ['..\' caseFolder '\NuMAD\airfoils'])
    copyfile('..\AeroData', ['..\' caseFolder '\AeroData'])    
    copyfile('..\out\IECSweep_ramp.out', ['..\' caseFolder '\out\IECSweep_ramp.out'])
    copyfile('..\wind\ramp.wnd', ['..\' caseFolder '\wind\ramp.wnd'])
    copyfile('..\runIEC_ipt.m',['..\' caseFolder '\runIEC_ipt.m'])
    copyfile('tempFail.mac',['..\' caseFolder '\NuMAD\tempFail.mac']) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Temporary
    cd(['..\' caseFolder '\NuMAD'])
    
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Blade Design Variables and Material Properties
% material limits -- unidirectional glass fiber, 57% Vf
strain_design_GFRP_uni = [-7900 12440];
% material limits -- Zoltek PX35 pultrusions, 68% Vf
strain_design_CFRP_PX35 = [-5670 8290];
% material limits -- CFTF Kaltex pultrusions, 68% Vf
strain_design_CFRP_K20 = [-4270 4890];

% assumes the spar cap material is same for LP and HP sides
sparCmp = find(contains({blade.components.name},'Spar_cap_ss'));
sparMaterial = blade.materials(blade.components(sparCmp).materialid).name;
if strcmp(sparMaterial,'glass_uni')
    sparStrainDesign_min = strain_design_GFRP_uni(1);
    sparStrainDesign_max = strain_design_GFRP_uni(2);
    fatigueMaterial = 's1_fiberglass';
elseif strcmp(sparMaterial,'CFP-baseline')
    sparStrainDesign_min = strain_design_CFRP_PX35(1);
    sparStrainDesign_max = strain_design_CFRP_PX35(2);
    fatigueMaterial = 's2_baselineCF';
elseif strcmp(sparMaterial,'CFP-heavytextile')
    sparStrainDesign_min = strain_design_CFRP_K20(1);
    sparStrainDesign_max = strain_design_CFRP_K20(2);
    fatigueMaterial = 's3_heavyTCF';
else
    error('new material used')
end

% assumes TE and LE reinforcement is fiberglass
reinfStrainDesign_min = strain_design_GFRP_uni(1);
reinfStrainDesign_max = strain_design_GFRP_uni(2);

%%%%%%%%%%%%%%%%%%%%%%%%% set targets for design %%%%%%%%%%%%%%%%%%%%%%%%%%
designcriteria.MinHPStrain=sparStrainDesign_min;
designcriteria.MaxHPStrain=sparStrainDesign_max;
designcriteria.MinLPStrain=sparStrainDesign_min;
designcriteria.MaxLPStrain=sparStrainDesign_max;
designcriteria.MinEdgeStrain=reinfStrainDesign_min;
designcriteria.MaxEdgeStrain=reinfStrainDesign_max;
 
designcriteria.Buckle=1.68;   % it's high because I'm using a course mesh
%designcriteria.Failure=4000*1e-6;   % it's high because I'm using a course mesh
designcriteria.Failure=1/2.06; 
designcriteria.MaxOoPDefl = 17;   % m
%%%%%%%%%%%%%%%%%%%%%%% set ga case blade properties %%%%%%%%%%%%%%%%%%%%%%


% 
%     ct=numel(compIndex);
%     ub=1.1*ones(ct,1);
%     lb=0.8*ones(ct,1);
% 
% 
% %Plot upper and lower bounds
% ct=1;
% ctMatch=1;
% for i=compIndex
%     xplot=blade.components(i).cp(:,1);
%     figure(1)
%     subplot(4,4,i-1)
%     plot(xplot,blade.components(i).cp(:,2),'ok--')
%     y=lb(ct)*blade.components(i).cp(:,2);
%     y(y<0)=0;
%     plot(xplot,y,'rd') 
%     y=ub(ct)*blade.components(i).cp(:,2);
%     plot(xplot,y,'bd') 
%     ct=ct+1;
% end




disp('Creating bladeDef and writing NMD files...')

ct = 0; 

imatch=[2, 11; %Shell Skin
        12 14; %Aft web skin
        15 17]; %rear web skin
ctMatch=1;
for i=compIndex
%    if i~=imatch(ctMatch,2)
        ct = ct+1;
        blade.components(i).cp(:,2)=x(ct)*x0{ct};
        ctMatch=ctMatch+1;
%    end
end


% update numad input file
blade.updateGeometry
blade.updateKeypoints
blade.updateBOM


plotBladeDesignVar(blade,compIndex,'*-')

BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt');

if saveBlade   % save the blade file for future use - FINAL RUN ONLY.
    save blade blade
end

%%%%%%%%%%%%%%%% call precomp and generate FAST blade file %%%%%%%%%%%%%%%%
delete bmodesFrequencies.mat
%disp('Creating FAST Blade file using NuMAD and PreComp...')
% run one precomp call and generate new blade frequencies from revised
% numad.nmd input.
%numad('numad.nmd','precomp',[1 3 2])
%movefile('FASTBlade_Precomp.dat','..\')
%disp('Newly created FAST Blade file copied to FAST simulation directory')

%%%%%%%%%%%%%%%% if running ANSYS, generate FEA mesh %%%%%%%%%%%%%%%%
    
disp(' '); disp('Creating ANSYS model...')
fprintf('Mesh size setting = %0.4f\n',blade.mesh)
while 1
    %numad('numad.nmd','ansys',blade.mesh) 
      numad_multiLayer('numad.nmd','ansys',blade.mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEMP
    if exist('master.db','file')
        disp('ANSYS model created')
         
        break;
    end
    fprintf('%s: Waiting for NUMAD to make initial model in %s...\n',datestr(now), caseFolder)
    pause(3);
end

%% ************************************************************************
% ================= CREATE NODAL LISTING FROM ANSYS MESH =================
% run make_nlist.mac

ansysFilename = config.ansys.bucklingFile;
ansys_path = blade.paths.ansys;
ansys_product = 'ANSYS';

script_name='commands1.mac';
script_out='output1.txt';

fid=fopen(script_name,'w+');
fprintf(fid,'resume,master,db\n');
fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);   %From master, change the jobname

fprintf(fid,'!!! BEGIN MAKE_NLIST MACRO TEXT\n');
fprintf(fid,'ESEL,S,SEC,,1,999   \n');
fprintf(fid,'esel,u,type,,21'); %Unselect mass element
fprintf(fid,'NSLE,S  \n');
fprintf(fid,'/output,NLIST,lis\n');
fprintf(fid,'/page,1e6,,1e6,,\n');
fprintf(fid,'NLIST,ALL, , ,XYZ,NODE,NODE,NODE\n');
fprintf(fid,'/output,\n');

% fprintf(fid,'*get, nnode, node,0,count\n'); %Find the number (quantity) of nodes selected
% fprintf(fid,'*get, nel, ELEM,0,count\n'); %Find the number (quantity) of nodes selected
% 
% fprintf(fid,'/output,nel,txt\n');
% fprintf(fid,'*status, nel\n');
% fprintf(fid,'/output\n');
% 
% fprintf(fid,'/output,nnode,txt\n');
% fprintf(fid,'*status, nnode\n');
% fprintf(fid,'/output\n');

fprintf(fid,'!!! END MAKE_NLIST TEXT\n');
fclose(fid);

ansys_call = sprintf('"%s" -b -p %s -I %s -o %s',ansys_path,ansys_product,script_name,script_out);
disp(ansys_call)%ble
pse=3;
while 1
     [status,~] = dos(ansys_call);  % the windows system call to run the above ansys command
    if status==0
        disp(' ')
        disp('Node listing created')
        %delete(script_name);
        break;
    end
    fprintf('%s: Waiting for ANSYS NLIST operation in %s ...\n',datestr(now),caseFolder)
    pause(pse);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **************** CHECK RELEVANT IEC DESIGN LOAD CASES *******************
%disp('Checking Mass using FAST')
%bladeFilename = 'FASTBlade_precomp.dat';
%designvarMassPerLength=layupDesign_FASTmass_BLE(bladeFilename); 
%designVar.Mass = designvarMassPerLength * bladeLength;
% if designVar.Mass>designcriteria.Mass
%     penalty.Mass = (designVar.Mass/designcriteria.Mass)^2;
% else
%     penalty.Mass = 1;
% end

if 1
    %disp('Checking OoP deflection using FAST')
    DLCoptions = {'1.3'};
    %DLCoptions = {'1.3','1.2'};
    %% EMA: original:
%     output=layupDesign_FASTanalysis(blade,DLCoptions,useParallel);
    %% EMA: changed to:
    runFASTAnal = 0;
    output=layupDesign_FASTanalysis(blade,DLCoptions,runFASTAnal,useParallel);
    %% END
    % calculate maximum/minimum values of interest
    numDLC=length(output);
    
    fastVarNames = {'MaxOoPDefl' 'MaxHPStrain' 'MinHPStrain' ...
        'MaxLPStrain' 'MinLPStrain' 'MaxEdgeStrain' 'MinEdgeStrain'};
    
    for jj = 1:length(fastVarNames)
        % =============== FIND EXTREME FAST DATA VALUES =============== %
        designVarName = fastVarNames{jj}; fastVar=zeros(1,numDLC);
        % add the partial safety factor for critical displacement
        if contains(designVarName,'OoP')            
            FSdeflection = 1.1;
        else, FSdeflection = 1.0; 
        end
        for vv=1:length({output.(designVarName)})
            %% EMA: original:
%             normalFS = {'1p1' '1p3' '1p4' '1p5' '2p1' '3p2' '3p3' '4p2' '5p1' '6p1' '6p3'};
            %% changed to:
            normalFS = {'1p1' '1p2' '1p3' '1p4' '1p5' '2p1' '3p2' '3p3' '4p2' '5p1' '6p1' '6p3'};
            %% END
            abnormalFS = {'2p2' '2p3' '6p2' '7p1'};
            % add the specified partial safety factor for loads
            % NOTE: not adding critical deflection partial safety factor
            if contains(output(vv).(designVarName).Name, normalFS)
                FSloads = 1.35; % normal safety factor (IEC 61400-1)
            elseif contains(output(vv).(designVarName).Name, abnormalFS)
                FSloads = 1.1; % abnormal safety factor (IEC 61400-1)
            else
                error('design load case not recognized')
            end
            fastVar(vv)=output(vv).(designVarName).data * FSloads * FSdeflection;
        end
        if contains(designVarName,'max','IgnoreCase',1) || contains(designVarName,'ampl','IgnoreCase',1) % find maximum
            designVar.(designVarName)=max(fastVar);
            
        elseif contains(designVarName,'min','IgnoreCase',1) % find minimum
            designVar.(designVarName)=min(fastVar);
        else
            error('FAST variable name not correct')
        end
        % =============== ADD PENALTY =============== %
        penalty.(designVarName) = max(designVar.(designVarName)/designcriteria.(designVarName), 1)^2;
    end
else
    designVar.MaxOoPDefl = 999;
    penalty.MaxOoPDefl = 1;
end


if 0 % ================= Perform Fatigue Analysis =================
    disp('Checking fatigue values using FAST')
    runFASTfatigue = 0; % flag to re-run FAST simulations    
    damage=layupDesign_FASTFatigue(blade,runFASTfatigue,useParallel);
    designVar.flapFatigue = max(damage.flap.(fatigueMaterial)); % spar cap fatigue only
    designVar.edgeFatigue = max(damage.edge.(fatigueMaterial));     
    penalty.flapFatigue = max(1,designVar.flapFatigue^2);
    penalty.edgeFatigue = max(1,designVar.edgeFatigue^2);
else
    designVar.flapFatigue = 999;
    designVar.edgeFatigue = 999;  
    penalty.flapFatigue = 1;
    penalty.edgeFatigue = 1;
end

if 1 % ================= Perform Buckling Analysis =================


    % determine the loads to apply to the FEA mesh for maximum sectional moment   
    % determine the loads to apply to the FEA mesh for maximum sectional moment   
    halfdz=2.5; %for best results halfdz should be a multiple of L and should be <= L/2
    bladeLength = blade.ispan(end);
    r=(halfdz:2*halfdz:bladeLength)';
    loads_table = FastLoads4ansys(output,fast_gage,r);
    %loads_table = FastLoads4ansys_ec(output,fast_gage,blade);
    save loads_table
%     buckling_loadFactor=zeros(1,length(loads_table));
%     failIndex=zeros(1,length(loads_table));
    
    % perform buckling check for maxima/minima at each rotated direction
    buckling_loadFactor=zeros(length(loads_table),1);
    failIndex=zeros(length(loads_table),1);
    MaxOoPDefl =zeros(length(loads_table),1);
    forcesfile={'forces-0.forces','forces-90.forces','forces-180.forces','forces-270.forces'};
    for ii = 1:length(loads_table)
        disp('Checking buckling criteria using ANSYS')
        
        %loads = forcespforcesToLoadsTable(forcesfile{ii}) 
       % [feaOutput]=layupDesign_ANSYSbuckling_ble(blade,loads,config,ii);
        [feaOutput]=layupDesign_ANSYSbuckling_ble(blade,loads_table{ii},config,ii,caseFolder);
        
        loadVSspan=getANSYSresultantLoadsVSspan(blade.ispan(end),blade.mesh)
        
        figure(2000+ii)
        subplot(2,2,1)%mx
        plot(loadVSspan.zloc,loadVSspan.Mysum/1000,'b') %Plot ANSYS My as FAST Mx
        subplot(2,2,2)%fy
        plot(loadVSspan.zloc,-loadVSspan.Fxsum/1000,'b') 
        subplot(2,2,3)%my
        plot(loadVSspan.zloc,-loadVSspan.Mxsum/1000,'b') %Plot ANSYS Mx as FAST My
        subplot(2,2,4)%fx
        plot(loadVSspan.zloc,loadVSspan.Fysum/1000,'b') 
        %legend('Applied','Original','ANSYS Resultant Force/Moment Sum')








        buckling_loadFactor(ii)=feaOutput.globalbuckling;
        if config.ansys.wrinkling
            wrinkling_loadFactor(ii)=feaOutput.wrinkling;
        end
        
        failIndex(ii)=read_1_ANSYSoutput('fail.txt');
        
        MaxOoPDefl(ii)=read_1_ANSYSoutput('uymax.txt');
        T=table(failIndex,buckling_loadFactor,MaxOoPDefl)
    end
    % calculate the global buckling penalty using the lowest buckling load factor
    designVar.Buckle = min(buckling_loadFactor);
    designVar.Failure = max(failIndex);
    designVar.MaxOoPDefl= max(MaxOoPDefl);
    
    
    if designVar.Buckle<designcriteria.Buckle
        penalty.Buckle = (designcriteria.Buckle/designVar.Buckle)^2;
    else
        penalty.Buckle = 1;
    end
    
    if designVar.Failure>designcriteria.Failure
        penalty.Failure = (designVar.Failure/designcriteria.Failure)^2;
    else
        penalty.Failure = 1;
    end
    
    if designVar.MaxOoPDefl>designcriteria.MaxOoPDefl
        penalty.MaxOoPDefl = (designVar.MaxOoPDefl/designcriteria.MaxOoPDefl)^2;
    else
        penalty.MaxOoPDefl = 1;
    end
    
    designVar.Mass=read_1_ANSYSoutput('mass.txt');
    
    % calculate the shell wrinkling penalty using the lowest wrinkling load factor
%  SKIP WRINKLE FOR NOW
%     designVar.Wrinkle = min(wrinkling_loadFactor);
%     if designVar.Wrinkle<designcriteria.Wrinkle
%         penalty.Wrinkle = (designcriteria.Wrinkle/designVar.Wrinkle)^2;
%     else
%         penalty.Wrinkle = 1;
%     end    
else
    designVar.Buckle = 999;
    designVar.Wrinkle = 999;
    penalty.Buckle = 1;
    penalty.Wrinkle = 1;
end




%%%%%%%%%%%%%%%% END check relevant IEC Design Load Cases %%%%%%%%%%%%%%%%%
delete 'master.db'
% *** {'MaxOoPDefl' 'MaxFlapStrain' 'MinFlapStrain' 'MaxEdgeStrain' 'MinEdgeStrain' } *** %
% 
% f=designVar.Mass*penalty.MaxOoPDefl*penalty.MaxHPStrain*penalty.MinHPStrain*...
%     penalty.MaxLPStrain*penalty.MinLPStrain*penalty.MaxEdgeStrain*penalty.MinEdgeStrain*...
%     penalty.flapFatigue*penalty.Buckle*penalty.Failure;


% f=designVar.Mass*penalty.MaxOoPDefl*...
%     penalty.flapFatigue*penalty.Buckle*penalty.Failure;

f=designVar.Mass*penalty.MaxOoPDefl*...
    penalty.Buckle*penalty.Failure;

% params=[designVar.Failure,penalty.Failure,designVar.Buckle,penalty.Buckle,designVar.FFreq,penalty.FFreq,designVar.EFreq,penalty.EFreq,designVar.MaxOoPDefl,penalty.MaxOoPDefl,...
%     designVar.MaxHPStrain,penalty.MaxHPStrain,designVar.MinHPStrain,penalty.MinHPStrain,...
%     designVar.MaxLPStrain,penalty.MaxLPStrain,designVar.MinLPStrain,penalty.MinLPStrain,...
%     designVar.MaxEdgeStrain,penalty.MaxEdgeStrain,designVar.MinEdgeStrain,penalty.MinEdgeStrain,...
%     designVar.flapFatigue,penalty.flapFatigue,designVar.edgeFatigue,penalty.edgeFatigue,designVar.Mass,f];

params=[designVar.Failure,penalty.Failure,designVar.Buckle,penalty.Buckle,designVar.MaxOoPDefl,penalty.MaxOoPDefl,...
    designVar.flapFatigue,penalty.flapFatigue,designVar.edgeFatigue,penalty.edgeFatigue,designVar.Mass,f];

% save the blade object for each iteration
save_index = 0;
while save_index < 10
    try
        save blade blade
        break
    catch
        % this file may be open outside of this worker, pause and retry
        pause(0.1)
        save_index = save_index+1;
    end
    if save_index == 10
        error('saving progress in parallel is problematic')
    end
end

% write results to layupDesignProgress.txt
cd(hm)
save_index = 0;
while save_index < 10
    try
        fid=fopen('layupDesignProgress.txt','a');
         %fprintf(fid,'%s,(value/penalty),FailureIndex, %2.4f,%2.4f,Buckling LF, %2.4f,%2.4f,Flap Freq,%2.4f,%2.2f,Edge Freq,%2.4f,%2.2f,OoPDefl,%2.4f,%2.2f,HPStrain,%.0f,%2.2f,HPStrain,%.0f,%2.2f,LPStrain,%.0f,%2.2f,LPStrain,%.0f,%2.2f,EdgeStrain,%.0f,%2.2f,EdgeStrain,%.0f,%2.2f,Flap Fatigue,%2.6f,%2.2f,Edge Fatigue,%2.6f,%2.2f,Actual mass kg,%.0f,fitness function,%.0f\n',datestr(now),params);
        fprintf(fid,'\n%s,(value/penalty),FailureIndex, %2.4f,%2.4f,Buckling LF, %2.4f,%2.4f,OoPDefl,%2.4f,%2.2f,Flap Fatigue,%2.6f,%2.2f,Edge Fatigue,%2.6f,%2.2f,Actual mass kg,%.0f,fitness function,%.0f',datestr(now),params);
        fclose(fid);
        break
    catch
        % this file may be open outside of this worker, pause and retry
        pause(0.05)
        save_index = save_index+1;
    end
    if save_index == 10
        error('saving progress in parallel is problematic')
    end
end

% write the candidate to file...
save_index = 0;
candidate_data = [f x];

str='';
ct = 0;
for i=2:17
  Nx=length(blade.components(i).cp)-1 ;

  str=[str repmat(', %.2f',1,Nx) '||'];
end
while save_index < 10
    try
        fid = fopen('layupDesignCandidates.txt','a');
        %fprintf(fid,'%s,fitness function,%6.0f, designVars, %.2f, %.2f, %.2f, %.2f, %.2f, %.2f\n',datestr(now),candidate_data);
        fprintf(fid,['\n%s,fitness function,%6.0f, designVars' str],datestr(now),candidate_data);
        fclose(fid);
        break
    catch
        % this file may be open outside of this worker, pause and retry
        pause(0.05)
        save_index = save_index+1;
    end
    if save_index == 10
        error('saving progress in parallel is problematic')
    end
end

    
end
%%
function [initPop,initScores] = setInitPop(useRestartFile,restartFile)

if useRestartFile
    a=load(restartFile);
    pp= a(:, 1) == a(end, 1) ;
    aa=sortrows( a(pp,:), size(a,2) );
    initPop=aa( 1:6, 2:end-1 )
    initScores=aa( 1:6, end )
else
    % manual entry
    initPop    =[];
    initScores =[];
end

end
%%
function [state, stop] = optCustomPlot(optimValues,state)
%options, state, and flag are delivered from the MATLAB genetic algorithm functions
disp('IN OPTIMIZATION CUSTOM PLOT FUNCTION...')
% state
% optimValues
% 
% % % gennum = state.Generation;
% 
stop = false; % is this needed?
% 
% switch state
%     
%     case 'init'
%         % do not save swarm
%         nplot = size(optimValues.swarm,2); % Number of dimensions
%         for i = 1:nplot % Set up axes for plot
%             subplot(nplot,1,i);
%             tag = sprintf('psoplotrange_var_%g',i); % Set a tag for the subplot
%             semilogy(optimValues.iteration,0,'-k','Tag',tag); % Log-scaled plot
%             ylabel(num2str(i))
%         end
%         xlabel('Iteration','interp','none'); % Iteration number at the bottom
%         subplot(nplot,1,1) % Title at the top
%         title('Log range of particles by component')
%         setappdata(gcf,'t0',tic); % Set up a timer to plot only when needed
%         
%     case 'iter'
%         % save the swarm
%         nplot = size(optimValues.swarm,2); % Number of dimensions
%         for i = 1:nplot
%             subplot(nplot,1,i);
%             % Calculate the range of the particles at dimension i
%             irange = max(optimValues.swarm(:,i)) - min(optimValues.swarm(:,i));
%             tag = sprintf('psoplotrange_var_%g',i);
%             plotHandle = findobj(get(gca,'Children'),'Tag',tag); % Get the subplot
%             xdata = plotHandle.XData; % Get the X data from the plot
%             newX = [xdata optimValues.iteration]; % Add the new iteration
%             plotHandle.XData = newX; % Put the X data into the plot
%             ydata = plotHandle.YData; % Get the Y data from the plot
%             newY = [ydata irange]; % Add the new value
%             plotHandle.YData = newY; % Put the Y data into the plot
%         end
%         if toc(getappdata(gcf,'t0')) > 1/30 % If 1/30 s has passed
%             drawnow % Show the plot
%             setappdata(gcf,'t0',tic); % Reset the timer
%         end
%         
        
        % write the candidate to file...
        fid = fopen('layupDesignCandidates.txt','a');
        for ii=1:size(optimValues.swarm,1)
            fprintf( fid, ' %.0f', optimValues.iteration );
            fprintf( fid, ' %f', optimValues.swarm(ii,:) );
            fprintf( fid, ' %f', optimValues.swarmfvals(ii) );
            fprintf( fid, '\r\n' );
        end
        %         fclose( fid ) ;
        %     case 'done'
%         % do nothing
% end


% optimValues = 
% 
%   struct with fields:
% 
%            bestfval: 2.2959e+04
%               bestx: [609.4776 125.7200 91.7267 89.1016 60.4185]
%           iteration: 0
%           funccount: 64
%            meanfval: 8.8377e+04
%     stalliterations: 0
%               swarm: [64×5 double]
%          swarmfvals: [64×1 double]



end
%%
function state = gaCustomPlot(options,state,flag)
%options, state, and flag are delivered from the MATLAB genetic algorithm functions

gennum = state.Generation;

% write the candidate to file...
fid = fopen('layupDesignCandidates.txt','a');
for jj=1:size(state.Population,1)
    fprintf( fid, ' %f', state.Generation );
    for ii = 1:size(state.Population,2)
        fprintf( fid, ' %f', state.Population(jj,ii) );
    end
    fprintf( fid, ' %f', state.Score(jj) );
    fprintf( fid, '\r\n' );
end
fclose( fid ) ;

end


function plotBladeDesignVar(blade,compIndex,sym)
ctMatch=1;
    imatch=[2, 11; %Shell Skin
            12 14; %Aft web skin
            15 17]; %rear web skin
    for i=compIndex
        if i~=imatch(ctMatch,2)
        figure(1)
        subplot(4,4,i-1)
        plot(blade.components(i).cp(:,1),blade.components(i).cp(:,2),sym)
        hold on
        else
            ctMatch=ctMatch+1;
        end

    end
end
