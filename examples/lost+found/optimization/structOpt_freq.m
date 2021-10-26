function structOpt_freq(blade,lb,ub,A,b,pop,gen,restartFile,useParallel)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************** INPUT OPTIONS ********************************
if isempty(restartFile)
    useRestartFile = false;
    delete 'layupDesignCandidates.txt'
else
    useRestartFile = true;
end

[initPop,initScores] = setInitPop(useRestartFile,restartFile);
if useRestartFile
    options = gaoptimset( 'PopulationSize',pop,...
        'Generations',gen,...
        'InitialPopulation',initPop,...
        'InitialScores',initScores,...
        'UseParallel',useParallel,...
        'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
else
    options = gaoptimset('PopulationSize',pop,...
        'Generations',gen,...
        'UseParallel',useParallel,...
        'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotbestindiv,@gaCustomPlot} );
end

nVar=length(lb);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ***************** DELETE OUTPUT FROM PREVIOUS RUNS **********************
delete('..\wind\*')  % delete all files in the wind directory
delete('..\carray.mat')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ************************** RUN OPTIMIZATION *****************************
blade.mesh=0.1;

if 1 % call the genetic algorithm
    % delete the former Progress file
    delete 'layupDesignProgress.txt'
    % call the genetic algorithm
    [x,fval,exitflag,output,population,scores] = ga(@(x)designBladeFcn(x,blade,useParallel, false),...
        nVar,A,b,[],[],lb,ub,[],[],options);
    save gaResults x fval exitflag output population scores
    
else % just do manual analysis and then quit
    x=[16.977187 0.056196 60.966734 23.245862 298.698868 0.169861 2.608083];
    for i=1:size(x,1)
        f(i,1) = designBladeFcn(x(i,:), blade, false, true);
        save tmp x f
    end
    disp('results summary:')
    sprintf('%i %i %i %i %i %i %.1f\n',[x f]')
end

end
%%
function f = designBladeFcn(x, blade, useParallel, saveBlade)

% <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Change folder saving structure so that multiple workers in a parallel
% algorithm don't overwrite each other.
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
    caseFolder = sprintf(['parallel_worker_' id]);
    mkdir(['..\' caseFolder '\NuMAD'])
    % if exist(['..\' caseFolder],'dir'), rmdir(['..\' caseFolder],'s'), end
    copyfile('airfoils', ['..\' caseFolder '\NuMAD\airfoils'])
    cd(['..\' caseFolder '\NuMAD'])
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%%%%%%%%%%%%%%%%%%%%%%%%% set targets for design %%%%%%%%%%%%%%%%%%%%%%%%%%
designcriteria_OoPDefl=1;
designcriteria_Mass = 650;   % kg
designcriteria_FFreq = 43/60 * 3 * 1.1;  % for 13m rotor
designcriteria_EFreq = 43/60 * 3 * 1.1;  % for 13m rotor
% designcriteria_Buckle=1.00;  % for use with final mesh
% designcriteria_Buckle=1.62;  % for use with final mesh
% designcriteria_Buckle=1.68;  % it's high because I'm using a course mesh

%%%%%%%%%%%%%%%%%%%%%%% set ga case blade properties %%%%%%%%%%%%%%%%%%%%%%
disp('Creating bladeDef and writing NMD files...')

% spar cap width
blade.sparcapwidth=x(1);
blade.sparcapwidth
% root parameters
cmp=5;
blade.components(cmp).name
blade.components(cmp).cp(3,1)=x(2); % length
blade.components(cmp).cp(3,2)=x(3); % thickness
blade.components(cmp).cp(4,1)=x(4); % length
blade.components(cmp).cp
% spar cap thickness
cmp=7;
blade.components(cmp).name
blade.components(cmp).cp(1,2)=x(5);
blade.components(cmp).cp(2,2)=x(6);
blade.components(cmp).cp(3,2)=x(7);
blade.components(cmp).cp(4,2)=x(8);
blade.components(cmp).cp(4,1)=0.95;
blade.components(cmp).cp
% % panel parameters
% cmp=6;
% blade.components(cmp).name
% blade.components(cmp).cp(3,2)=x(1);
% blade.components(cmp).cp(4,2)=x(1);
% blade.components(cmp).cp(5,2)=x(2);
% blade.components(cmp).cp(6,2)=x(2);
% blade.components(cmp).cp
% % te-reinf parameters
% cmp=7;
% blade.components(cmp).name
% blade.components(cmp).cp(2,2)=x(6);
% blade.components(cmp).cp(3,2)=x(6);
% blade.components(cmp).cp

% update numad input file
blade.updateGeometry
blade.updateKeypoints
try
blade.updateBOM
catch
    x
    keyboard
    % why did this error occur? - ble.  Error check fixed in bladeDef.m
end

BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt');

if saveBlade   % save the blade file for future use - FINAL RUN ONLY.
    keyboard
    save blade blade
end

%%%%%%%%%%%%%%%% call precomp and generate FAST blade file %%%%%%%%%%%%%%%%
delete bmodesFrequencies.mat
disp('Creating FAST Blade file using NuMAD and PreComp...')
% run one precomp call and generate new blade frequencies from revised numad.nmd input.
numad('numad.nmd','precomp',[1 3 2])
disp('Newly created FAST Blade file copied to FAST simulation directory')
movefile('FASTBlade_Precomp.dat','..')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% **************** CHECK RELEVANT IEC DESIGN LOAD CASES *******************
disp('Checking Mass using FAST')
bladeFilename = 'FASTBlade_precomp.dat';
designvarMassPerLength=layupDesignFASTmass_BLE(bladeFilename); 
designvarMass = designvarMassPerLength * 13;
if designvarMass>designcriteria_Mass
    penaltyMass = (designvarMass/designcriteria_Mass)^2;
else
    penaltyMass = 1;
end

% disp('Checking OoP deflection using FAST')
% designvarOoPDefl=layupDesignFAST(blade,useParallel);
% if designvarOoPDefl>designcriteria_OoPDefl
%     penaltyOoPDefl = (designvarOoPDefl/designcriteria_OoPDefl)^2;
% else
%     penaltyOoPDefl = 1;
% end
designvarOoPDefl = designcriteria_OoPDefl;
penaltyOoPDefl = 1;

disp('Checking frequencies using NuMAD, PreComp, and BModes...')
[~,freqs]=layupDesignFASTBlade(blade);
designvarFFreq=freqs(1);
if designvarFFreq<designcriteria_FFreq
    penaltyFFreq = (designcriteria_FFreq/designvarFFreq)^2;
else
    penaltyFFreq = 1;
end
designvarEFreq=freqs(2);
if designvarEFreq<designcriteria_EFreq
    penaltyEFreq = (designcriteria_EFreq/designvarEFreq)^2;
else
    penaltyEFreq = 1;
end

% disp('Checking fatigue values using FAST')
% tmp=layupDesignFASTFatigue(blade);
% designvarFatigue=max(tmp(11:20,2)); % spar cap fatigue only
% if designvarFatigue>1
%     penaltyFatigue = designvarFatigue^2;
% else
%     penaltyFatigue = 1;
% end

% disp('Checking buckling criteria using ANSYS')
% [designvar,mass]=layupDesignANSYS(blade)
% designvarBuckle=designvar(1);
% if designvarBuckle<designcriteria_Buckle
%     penaltyBuckle = (designcriteria_Buckle/designvarBuckle)^2;
% else
%     penaltyBuckle = 1;
% end

%%%%%%%%%%%%%%%% END check relevant IEC Design Load Cases %%%%%%%%%%%%%%%%%
delete 'master.db'

% Maximum frequency design, ignoring fatigue, ignoring buckling, 
f=(-1*designvarFFreq)*1/penaltyMass*1/penaltyOoPDefl*1/penaltyEFreq;
params=[f,designvarFFreq,designvarOoPDefl,designvarMass];

% write results to layupDesignProgress.txt
cd(hm)
fid=fopen('layupDesignProgress.txt','a');
fprintf(fid,'%s,(value/penalty) %2.2f; Flap Freq [hz] %2.3f; OoPDefl [m] %2.3f; Actual mass [kg] %6.0f\n',datestr(now),params);
fclose(fid);

end
%%
function [initPop,initScores] = setInitPop(useRestartFile,restartFile)

if useRestartFile
    a=load(restartFile);
    pp= a(:, 1) == a(end, 1) ;
    aa=sortrows( a(pp,:), size(a,2) );
    initPop=aa( 1:7, 2:end-1 )
    initScores=aa( 1:7, end )
else
    % manual entry
    initPop    =[];
    initScores =[];
end

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