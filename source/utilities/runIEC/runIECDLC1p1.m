function output=runIECDLC1p1(params,output)

global fastPath
global adamsPath  
global turbsimPath



CaseName='IECDLC1p1NTM';

% Switches:
% % params.fastsim='fast simulink';

if params.simulate
    runturbsim=1;runfast=1;saveoutdata=1;
else
    runturbsim=0;runfast=0;saveoutdata=1;
end

load('seeds')

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
'perform a reassignment of the wind speed vector and number of seeds if';
'this should be different from the other DLC';
% % params.ws = [5:1:20];
% % params.numSeeds = 24;
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% % switch params.Class
% %     case 1
% %         avgws=0.2*50; % m/s, average wind speed of IEC Class I site (Vref=50m/s); IEC Section 6.3.1.1 Eqn (9)
% %     case 2
% %         avgws=0.2*42.5; % m/s, average wind speed of IEC Class II site (Vref=42.5m/s); IEC Section 6.3.1.1 Eqn (9)
% %     case 3
% %         avgws=0.2*37.5; % m/s, average wind speed of IEC Class III site (Vref=37.5m/s); IEC Section 6.3.1.1 Eqn (9)
% % end
% %
% %
% % % ====================== MAIN PROCESS - SIMULATIONS =======================
% % % create a counter for parallel runs
% % for ii = 1:length(params.ws)
% %     ctr.w((ii-1)*params.numSeeds+1:ii*params.numSeeds) = params.ws(ii);
% %     ctr.s((ii-1)*params.numSeeds+1:ii*params.numSeeds) = 1:params.numSeeds;
% % end
% % ctr.w1 = ctr.w; ctr.s1 = ctr.s;
% % yawMat = repmat(params.yaw, params.numSeeds*length(params.ws), 1);
% % ctr.y = reshape(yawMat, 1, numel(yawMat));
% % ctr.w = repmat(ctr.w, 1, length(params.yaw));
% % ctr.s = repmat(ctr.s, 1, length(params.yaw));

% % if runturbsim
% %     parfor ii = 1:length(ctr.w1)
% %         % Run TurbSim to generate the turbulent input files
% %         fst=readFastMain(['IEC_' params.fstfn '.fst']); % ble: WARNING - IS IT OK TO READ THE SAME FILE IN PARALLEL??
% %         tsim=readTurbSim([params.parDir 'init\wind.inp']);
% %         hm=pwd;
% %         cd([params.parDir 'wind'])
% %         % generate turbsim file
% %         tsim.RandSeed1=params.seeds(ctr.s1(ii));%ble %#ok<PFBNS>
% %         tsim.URef=ctr.w1(ii);
% %         tsim.IECturbc=params.TurbClass;
% %         tsim.IEC_WindType='NTM';
% %         tsim.RefHt= fst.TurbConf.TowerHt+fst.TurbConf.Twr2Shft+sind(-1*fst.TurbConf.ShftTilt)*-1*fst.TurbConf.OverHang;
% %         tsim.AnalysisTime=params.SimTime;
% %         tsim.UsableTime=tsim.AnalysisTime;
% %         tsim.NumGrid_Z=params.NumGrid;
% %         tsim.NumGrid_Y=tsim.NumGrid_Z;
% %         tsim.HubHt=tsim.RefHt;
% %         tsim.GridHeight=fst.TurbConf.TipRad*2*1.15;
% %         tsim.GridWidth=tsim.GridHeight;
% %         writeTurbSim(tsim,['turb_DLC1p1_' num2str(ctr.w1(ii)) 'mps_seed' num2str(ctr.s1(ii)) '.inp']);
% %         dos([turbsimPath ' turb_DLC1p1_' num2str(ctr.w1(ii)) 'mps_seed' num2str(ctr.s1(ii)) '.inp'],'-echo');
% %         cd(hm)
% %     end
% % end
% %
% %
% % parfor ii=1:numel(ctr.w)
% %     % set up a format so that files can be opened by parallel runs.
% %     t = getCurrentTask;
% %     if isempty(t)
% %         id = ''; % first iteration of genetic algorithm in generation
% %     else
% %         id = ['_parallelWorker_' num2str(t.ID)]; % remaining population are sent to workers in parallel
% %         if runfast
% %             % copy the FAST and AD input files for parallel runs
% %             disp([params.fstfn '.fst'])
% %             copyfile([params.fstfn '.fst'], [params.fstfn id '.fst'])
% %             fst=readFastMain(['IEC_' params.fstfn id '.fst']);
% %             copyfile(strrep(fst.ADFile,'"',''), [fst.ADFile(2:end-8) id '_AD.ipt'])
% %         end
% %     end
% %
% % %     if runturbsim
% % %         % Run TurbSim to generate the turbulent input files
% % %         fst=readFastMain(['IEC_' params.fstfn id '.fst']);
% % %         tsim=readTurbSim([params.parDir 'init\wind.inp']);
% % %         hm=pwd;
% % %         cd([params.parDir 'wind'])
% % %         % generate turbsim file
% % %         tsim.RandSeed1=params.seeds(ctr.s(ii));%ble %#ok<PFBNS>
% % %         tsim.URef=ctr.w(ii);
% % %         tsim.IECturbc=params.TurbClass;
% % %         tsim.IEC_WindType='NTM';
% % %         tsim.RefHt= fst.TurbConf.TowerHt+fst.TurbConf.Twr2Shft+sind(-1*fst.TurbConf.ShftTilt)*-1*fst.TurbConf.OverHang;
% % %         tsim.AnalysisTime=params.SimTime;
% % %         tsim.UsableTime=tsim.AnalysisTime;
% % %         tsim.NumGrid_Z=params.NumGrid;
% % %         tsim.NumGrid_Y=tsim.NumGrid_Z;
% % %         tsim.HubHt=tsim.RefHt;
% % %         tsim.GridHeight=fst.TurbConf.TipRad*2*1.15;
% % %         tsim.GridWidth=tsim.GridHeight;
% % %         writeTurbSim(tsim,['turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.inp']);
% % %         dos([turbsimPath ' turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.inp'],'-echo');
% % %         cd(hm)
% % %     end
% %
% %     if runfast
% %         % prepare and perform aeroelastic simulation
% %         fst=readFastMain(['IEC_' params.fstfn id '.fst']);
% %         fst.ADFile=['"' CaseName '_' params.fstfn id '_AD.ipt' '"'];% ble: change the AD filename
% %         ad.WindFile=['"' params.parDir 'wind\turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.wnd"'];
% %         fst.SimCtrl.TMax=params.SimTime;
% %         fst.ADFile=['"' CaseName '_' fst.ADFile(2:end-8) id '_AD.ipt' '"'];% ble: change the AD filename
% %         fst=runIEC_setIC(fst,ctr.w(ii),params.parDir);
% %         fst.Init.NacYaw = ctr.y(ii);
% %         thisFastName=[CaseName '_' params.fstfn id];%ble; change the FAST input filename
% %         if strcmp(fastadams, 'fast simulink')
% %             fst.TurbCtrl.YCMode   = 0; % no yaw control so that the yaw misalignment will be tested
% %             fst.TurbCtrl.PCMode   = 2;
% %             fst.TurbCtrl.VSContrl = 3;
% %             % Simulink implementation requires non-DOF initializations to 0
% %             fst.Init.OoPDefl = 0;
% %             fst.Init.IPDefl  = 0;
% %             fst.Init.TTDspFA = 0;
% %             fst.Init.TTDspSS = 0;
% %         else
% %             fst.TurbCtrl.YCMode   = 0; % no yaw control so that the yaw misalignment will be tested
% %             fst.TurbCtrl.PCMode   = 1;
% %             fst.TurbCtrl.VSContrl = 1;
% %         end
% %         writeFastMain(fst,[thisFastName '.fst']);
% %         writeFastAD(ad,strrep(fst.ADFile,'"',''));
% %
% %         switch fastadams
% %             case 'fast'
% %                 disp('Starting FAST model from Matlab.....')
% %                 dos([fastPath ' ' thisFastName '.fst']);
% %                 outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
% %                 try copyfile([thisFastName '.out'],['out/' outname]);
% %                 catch err; warning(err.message); end
% %                 try copyfile([thisFastName '.elm'],['out/' outname(1:end-4) '.elm']);
% %                 catch err; warning(err.message); end
% %                 disp('AD and FAST output files saved under new names.')
% %                 disp(' ')
% %                 disp('          *********** Done **************')
% %                 disp(' ')
% %             case 'fast simulink'
% %                 disp('Starting FAST SFunc model from Matlab.....')
% %                 fastFilenameSFunc = [thisFastName '.fst'];
% %                 runFASTusingSimulink(fastFilenameSFunc, params.simulinkModel)
% %                 outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
% %                 try movefile([thisFastName '_SFunc.out'],['out/' outname]);
% %                 catch err; warning(err.message); end
% %                 try movefile([thisFastName '_SFunc.elm'],['out/' outname(1:end-4) '.elm']);
% %                 catch err; warning(err.message); end
% %                 disp('AD and FAST_SFunc output files saved under new names.')
% %                 disp(' ')
% %                 disp('          *********** Done **************')
% %                 disp(' ')
% %             case 'adams'
% %                 fst.SimCtrl.ADAMSPrep='2';
% %                 fst.Init.OoPDefl=0;
% %                 fst.Init.IPDefl=0;
% %                 fst.Init.TeetDefl=0;
% %                 fst.Init.TTDspFA=0;
% %                 fst.Init.TTDspSS=0;
% %                 writeFastMain(fst,[thisFastName '_4ADAMS.fst']);
% %                 dos([fastPath ' ' thisFastName '_4ADAMS.fst']);
% %                 dos([adamsPath ' ' thisFastName '_4ADAMS_ADAMS.acf']);
% %                 disp('Starting ADAMS model from Matlab.....')
% %                 outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
% %                 try movefile([thisFastName '_4ADAMS_ADAMS.plt'], ['out/' outname])
% %                 catch err; warning(err.message); end
% %                 disp('AD and ADAMS output files saved under new names.')
% %                 disp(' ')
% %                 disp('          *********** Done **************')
% %                 disp(' ')
% %         end
% %     end
% % end


% ====================== MAIN PROCESS - SIMULATIONS =======================
% create a counter for parallel runs
for ii = 1:length(params.ws)
    ctr.w((ii-1)*params.numSeeds+1:ii*params.numSeeds) = params.ws(ii);
    ctr.s((ii-1)*params.numSeeds+1:ii*params.numSeeds) = 1:params.numSeeds;
end

% loop through "wind direction" as well
yawMat = repmat(params.yaw, params.numSeeds*length(params.ws), 1);
dirMat = repmat(params.wd, numel(yawMat), 1);
ctr.d = reshape(dirMat, 1, numel(dirMat));

for yaw = params.yaw
    for ii=1:(length(params.ws)*params.numSeeds)
        % set up a format so that files can be opened by parallel runs.
        t = getCurrentTask;
        if isempty(t)
            id = ''; % first iteration of genetic algorithm in generation
        else
            id = ['_parallelWorker_' num2str(t.ID)]; % remaining population are sent to workers in parallel
            if runfast
                % copy the FAST and AD input files for parallel runs
                % %                 disp([params.fstfn '.fst'])
                % %                 copyfile([params.fstfn '.fst'], [params.fstfn id '.fst'])
                % %                 fst=readFastMain(['IEC_' params.fstfn '.fst']);
                % %                 copyfile(strrep(fst.ADFile,'"',''), [fst.ADFile(2:end-8) id '_AD.ipt'])
            end
        end
        
        if runturbsim && yaw == params.yaw(1)
            % Run TurbSim to generate the turbulent input files
            fst=readFastMain(['IEC_' params.fstfn id '.fst']);
            tsim=readTurbSim([params.parDir 'init\wind.inp']);
            hm=pwd;
            cd([params.parDir 'wind'])
            % generate turbsim file
            tsim.RandSeed1=params.seeds(ctr.s(ii));%ble %#ok<PFBNS>
            tsim.URef=ctr.w(ii);
            tsim.IECturbc=params.TurbClass;
            tsim.IEC_WindType='NTM';
            tsim.RefHt= fst.TurbConf.TowerHt+fst.TurbConf.Twr2Shft+sind(-1*fst.TurbConf.ShftTilt)*-1*fst.TurbConf.OverHang;
            tsim.AnalysisTime=params.SimTime;
            tsim.UsableTime=tsim.AnalysisTime;
            tsim.NumGrid_Z=params.NumGrid;
            tsim.NumGrid_Y=tsim.NumGrid_Z;
            tsim.HubHt=tsim.RefHt;
            tsim.GridHeight=fst.TurbConf.TipRad*2*1.15;
            tsim.GridWidth=tsim.GridHeight;
            writeTurbSim(tsim,['turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.inp']);
            dos([turbsimPath ' turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.inp'],'-echo');
            cd(hm)
        end
        
        if runfast
            % prepare and perform aeroelastic simulation
            fst=readFastMain(['IEC_' params.fstfn '.fst']);%fst=readFastMain([params.fstfn id '.fst']);
            ad=readFastAD([fst.ADFile(2:end-8) '_AD.ipt']);%ad=readFastAD([fst.ADFile(2:end-8) id '_AD.ipt']);
            ad.WindFile=['"' params.parDir 'wind\turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.wnd"'];
            fst.SimCtrl.TMax=params.SimTime;
            fst.ADFile=['"' CaseName '_' params.fstfn id '_AD.ipt' '"'];% ble: change the AD filename
            if ~strcmp(params.fastsim,'fast simulink')
                fst=runIEC_setIC(fst,ctr.w(ii),params.parDir);
            end
            fst.Init.NacYaw = yaw;
            thisFastName=[CaseName '_' params.fstfn id];%ble; change the FAST input filename
            if strcmp(params.fastsim, 'fast simulink')
                fst.TurbCtrl.YCMode   = 2;
                fst.TurbCtrl.TYCOn = 0;
                fst.TurbCtrl.PCMode   = 2;
                fst.TurbCtrl.VSContrl = 3;
                % Simulink implementation requires non-DOF initializations to 0
                fst.Init.OoPDefl = 0;
                fst.Init.IPDefl  = 0;
                fst.Init.TTDspFA = 0;
                fst.Init.TTDspSS = 0;
            else
                fst.TurbCtrl.YCMode   = 0;
                fst.TurbCtrl.PCMode   = 1;
                fst.TurbCtrl.VSContrl = 1;
            end
            writeFastMain(fst,[thisFastName '.fst']);
            writeFastAD(ad,strrep(fst.ADFile,'"',''));
            
            switch params.fastsim
                case 'fast'
                    disp('Starting FAST model from Matlab.....')
                    dos([fastPath ' ' thisFastName '.fst']);
                    outname=[CaseName '_yaw' num2str(yaw) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
                    try movefile([thisFastName '.out'],['out/' outname]);
                    catch err; warning(err.message); end
                    try movefile([thisFastName '.elm'],['out/' outname(1:end-4) '.elm']);
                    catch err; warning(err.message); end
                    delete([thisFastName '.fsm'], [thisFastName '.opt']);
                    delete([thisFastName '.fst'], [thisFastName '_AD.ipt']);
                    disp('AD and FAST output files saved under new names.')
                    disp(' ')
                    disp('          *********** Done **************')
                    disp(' ')
                case 'fast simulink'
                    disp('Starting FAST SFunc model from Matlab.....')
                    fastFilenameSFunc = [thisFastName '.fst'];
                    YawPosBias = ctr.d(ii);
                    runFASTusingSimulink(fastFilenameSFunc, params.simulinkModel,YawPosBias)
                    outname=[CaseName '_yaw' num2str(yaw) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
                    try movefile([thisFastName '_SFunc.out'],['out/' outname]);
                    catch err; warning(err.message); end
                    try movefile([thisFastName '_SFunc.elm'],['out/' outname(1:end-4) '.elm']);
                    catch err; warning(err.message); end
                    delete([thisFastName '_SFunc.fsm'], [thisFastName '_SFunc.opt']);
                    delete([thisFastName '.fst'], [thisFastName '_AD.ipt']);
                    disp('AD and FAST_SFunc output files saved under new names.')
                    disp(' ')
                    disp('          *********** Done **************')
                    disp(' ')
                case 'adams'
                    fst.SimCtrl.ADAMSPrep='2';
                    fst.Init.OoPDefl=0;
                    fst.Init.IPDefl=0;
                    fst.Init.TeetDefl=0;
                    fst.Init.TTDspFA=0;
                    fst.Init.TTDspSS=0;
                    writeFastMain(fst,[params.fstfn id '_4ADAMS.fst']);
                    dos([fastPath ' ' params.fstfn id '_4ADAMS.fst']);
                    dos([adamsPath ' ' params.fstfn id '_4ADAMS_ADAMS.acf']);
                    disp('Starting ADAMS model from Matlab.....')
                    outname=[CaseName '_yaw' num2str(yaw) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.out'];
                    %                     try eval(['copyfile ' params.fstfn '_4ADAMS_ADAMS.plt ' 'out/' outname]);
                    try movefile([params.fstfn '_4ADAMS_ADAMS.plt'], ['out/' outname])
                    catch err; warning(err.message); end
                    disp('AD and ADAMS output files saved under new names.')
                    disp(' ')
                    disp('          *********** Done **************')
                    disp(' ')
            end
        end
    end
end


% ==================== 50-YEAR EXTRAPOLATION CALCULATIONS =================
% Channel labels for DLC 1.1 50-year extrapolation
labels = ['OoPDefl1' 'OoPDefl2' 'OoPDefl3' 'RootFxb1' 'RootFyb1' 'RootFzb1' 'RootMxb1' 'RootMyb1' 'RootMzb1' ...
    'RootFxb2' 'RootFyb2' 'RootFzb2' 'RootMxb2' 'RootMyb2' 'RootMzb2' ...
    'RootFxb3' 'RootFyb3' 'RootFzb3' 'RootMxb3' 'RootMyb3' 'RootMzb3' ...
    params.bladeGageLabels_FLx params.bladeGageLabels_FLy params.bladeGageLabels_FLz...
    params.bladeGageLabels_MLx params.bladeGageLabels_MLy params.bladeGageLabels_MLz];
bladeGageStations = 9;

% create a counter for reading in the data at one wind speed
for ii = 1:length(params.yaw)
    ctr2.y((ii-1)*params.numSeeds+1:ii*params.numSeeds) = params.yaw(ii);
    ctr2.s((ii-1)*params.numSeeds+1:ii*params.numSeeds) = 1:params.numSeeds;
end

warning('Change this section below to input Weibull parameters instead!!!!!')
%%%%%%%%%%%%%%%%% Weibull Wind PDF %%%%%%%%%%%%%%%%%%%%
load U_hub % ble: assign the mean wind speed.
% probability of wind and direction at hub
figure('position',[100  100 350 350])
set(gcf,'color','w');
bins = 1000;
[y, x] = hist(U_hub,bins);
area = trapz(x,y);
y = y./area;

W = repmat(1,1,bins);
bin_width = mean(diff(params.ws));
U_bin_edges = params.ws-bin_width./2;
U_bin_edges = [U_bin_edges params.ws(end)+bin_width./2];
[U_min U_rated_i] = min(abs(params.ws-11.1));
i_w = find(params.ws>10 & params.ws<14);
W(i_w) = 10;
% Weibull wind speed at hub height
f=fit(x', y',' k/l*((x-0)/l)^(k-1)*exp(-((x-0)/l)^k)', 'StartPoint', [2,7],'Weights',W); % RESULT: pdf of velocity at SWiFT.

for ww = 1:length(params.ws)
    term2(ww) = (f(U_bin_edges(ww)) + f(U_bin_edges(ww+1))).*bin_width./2;
end

N = length(params.yaw)*params.numSeeds*length(params.ws);
% Find and sort extreme values
outName={}; ctrWS = 0;
if saveoutdata % read in the FAST data files or load a saved structure with the loads
    for w = params.ws
        ctrWS = ctrWS+1;
        parfor ii = 1:length(params.yaw)*params.numSeeds
            % Read in the output files of interest for load calculations
            outName{ii}=['out/' CaseName '_yaw' num2str(ctr2.y(ii)) '_' num2str(w) 'mps_seed' num2str(ctr2.s(ii)) '.out'] %#ok<PFBNS>
            
            out=loadFASTOutData(outName{ii});  % load a file from the list
            for ll = 1:length(labels)
                dataVec = out.data(:,strcmp(out.list,labels{ll}));
                % save the maximums
                y_peakWS{ii}(ll,1) = max(dataVec);
                % save the minimums
                y_negPeakWS{ii}(ll,1) = min(dataVec);
            end
        end
        % save the variable for each wind speed
        for ll = 1:length(labels)
            tempPeak = [y_peakWS{:}];
            yMax_perWS{ll}(ctrWS,:) = tempPeak(ll,:);
            tempNegPeak = [y_negPeakWS{:}];
            yMin_perWS{ll}(ctrWS,:) = tempNegPeak(ll,:);
        end
    end
    
    % reshape the peaks so that all wind speeds are saved in one row vector
    for ll = 1:length(labels)
        yMax_AllWS{ll} = reshape(yMax_perWS{ll}',1,N);
        yMin_AllWS{ll} = reshape(yMin_perWS{ll}',1,N);
    end
    
    % save the summarized FAST output data
    saveLoads = 1;
    
    if saveLoads == 1
        save yMax_perWS yMax_perWS
        save yMin_perWS yMin_perWS
        save yMax_AllWS yMax_AllWS
        save yMin_AllWS yMin_AllWS
    elseif saveLoads == 2
        save yMax_perWS_10yaw yMax_perWS
        save yMin_perWS_10yaw yMin_perWS
        save yMax_AllWS_10yaw yMax_AllWS
        save yMin_AllWS_10yaw yMin_AllWS
    elseif saveLoads == 3
        save yMax_perWS_40yaw yMax_perWS
        save yMin_perWS_40yaw yMin_perWS
        save yMax_AllWS_40yaw yMax_AllWS
        save yMin_AllWS_40yaw yMin_AllWS
    elseif saveLoads == 4
        save yMax_perWS_ALLyaw yMax_perWS
        save yMin_perWS_ALLyaw yMin_perWS
        save yMax_AllWS_ALLyaw yMax_AllWS
        save yMin_AllWS_ALLyaw yMin_AllWS
    end
else
    % load the summarized FAST output data
    readLoads = 1;
    
    if readLoads == 1
        load yMax_perWS
        load yMin_perWS
        load yMax_AllWS
        load yMin_AllWS
    elseif readLoads == 2
        load yMax_perWS_10yaw
        load yMin_perWS_10yaw
        load yMax_AllWS_10yaw
        load yMin_AllWS_10yaw
    elseif readLoads == 3
        load yMax_perWS_40yaw
        load yMin_perWS_40yaw
        load yMax_AllWS_40yaw
        load yMin_AllWS_40yaw
    elseif readLoads == 4
        load yMax_perWS_ALLyaw
        load yMin_perWS_ALLyaw
        load yMax_AllWS_ALLyaw
        load yMin_AllWS_ALLyaw
    end
end



% find the indices corresponding to the root forces/moments for all blades
labelRootChannel = {'RootFxb' 'RootFyb' 'RootFzb' 'RootMxb' 'RootMyb' 'RootMzb'};
for ii = 1:length(labelRootChannel)
    name = labelRootChannel{ii};
    FMchannel.(name) = strfind(labels,name);
end
% find the indices corresponding to the span forces/moments for all blades
labelSpanChannel = {'FLxb' 'FLyb' 'FLzb' 'MLxb' 'MLyb' 'MLzb'};
for ss = 1:bladeGageStations
    for ii = 1:length(labelSpanChannel)
        name = ['Spn' num2str(ss) labelSpanChannel{ii}];
        FMchannel.(name) = strfind(labels,name);
    end
end
namesBldGag = fieldnames(FMchannel)
% find the indices corresponding to the out of plane deflection for all blades
labelDeflChannel = {'OoPDefl1' 'OoPDefl2' 'OoPDefl3'};
for ii = 1:length(labelDeflChannel)
    name = labelDeflChannel{ii};
    FMchannel.(name) = strfind(labels,name);
end

% find the indices corresponding to the similar forces/moments for all blades
names = fieldnames(FMchannel);
for ii = 1:length(names)
    FMchannelIndex.(names{ii}) = [];
    for ll = 1:length(labels)
        if ~isempty(FMchannel.(names{ii}){ll})
            FMchannelIndex.(names{ii})(end+1) = ll;
        end
    end
end

% call the function to calculate the 50-year extrapolated value
for ii = 1:length(names)
    dataVec{ii} = [];
    for jj = FMchannelIndex.(names{ii})
        dataVec{ii} = [dataVec{ii} yMax_AllWS{jj}];
        varStr = names{ii};
    end
    
    [load_50_Mxb.(names{ii}), e_load_50_Mxb.(names{ii}), F_50yearEstimate.(names{ii})] = performCalculationUsingFshort(dataVec{ii},...
        params, term2, N*length(FMchannelIndex.(names{ii})), varStr);
end


for ii = 1:length(namesBldGag) / 6
    for jj = 1:6
        array(ii,jj) = F_50yearEstimate.(namesBldGag{(ii-1)*6+jj})
    end
    rowNames{ii} = namesBldGag{(ii-1)*6+1}(1:4)
end
colNames = {'Fx' 'Fy' 'Fz' 'Mx' 'My' 'Mz'}
tblMax50yearFMvalues = array2table(array, 'VariableNames', colNames,'RowNames',rowNames)
writetable(tblMax50yearFMvalues, 'IECDLC_1p1_50yearLoads.xlsx','Sheet','BladeLoads','writeRowNames', true)

tblMaxOoPDefl50yearValues = table();
for ii = 1:length(labelDeflChannel)
    tblMaxOoPDefl50yearValues.(labelDeflChannel{ii}) = F_50yearEstimate.(labelDeflChannel{ii});
end
writetable(tblMaxOoPDefl50yearValues, 'IECDLC_1p1_50yearLoads.xlsx','Sheet','MaxOoPDefl','writeRowNames', true)


return
calculate the output maximums and minimums in a non-extrapolated manner
fst=readFastMain(['IEC_' params.fstfn '.fst']);
bld=readFastBlade(strrep(fst.BldFile{1},'"',''));
ad=readFastAD(strrep(fst.ADFile,'"',''));
BlLength=fst.TurbConf.TipRad-fst.TurbConf.HubRad;
spans=[fst.TurbConf.HubRad ad.RNodes(fst.Out.BldGagNd)'];
spans=spans'-spans(1)';

load carray
cs(:,1)=interp1(carray(:,2)*BlLength,carray(:,5),spans);  % edge bending c's
cs(:,2)=interp1(carray(:,2)*BlLength,carray(:,6),spans);  % flap bending c's

EI_flap=interp1(bld.prop.BlFract*BlLength,bld.prop.FlpStff,spans);
EI_edge=interp1(bld.prop.BlFract*BlLength,bld.prop.EdgStff,spans);
EIs=[EI_edge EI_flap];
EA_normal=interp1(bld.prop.BlFract*BlLength,bld.prop.EAStff,spans);

outName = {}
for w = params.ws
    for ii = 1:length(params.yaw)*params.numSeeds
        outName{end+1}=['out\' CaseName '_yaw' num2str(ctr2.y(ii)) '_' num2str(w) 'mps_seed' num2str(ctr2.s(ii)) '.out']; %#ok<AGROW,PFBNS>
    end
end

parfor cc = 1:length(outName)
    output2{cc} = sortExtremes(params,output,CaseName,outName(cc),cs,EIs,EA_normal);
end
% create the final max/min/ampl values from each of the simulations
names = fieldnames(output2{1});
for ii = 1:length(names) % loop through the channel variable
    % save the data values in a temporary vector
    for jj = 1:length(output2)
        tempVec(jj) = output2{jj}.(names{ii}).data;
    end
    % determine the index of the channel max/min/ampl value
    if ~isempty(strfind(names{ii},'Min'))
        % Compare output2 channels for minimum value
        [~,index] = min(tempVec);
    else
        % Compare output2 channels for maximum value (includes 'Ampl')
        [~,index] = max(tempVec);
    end
    output.(names{ii}) = output2{index}.(names{ii});
end


end
%%

function [load_50, e_load_50, F_50yearEstimate] = performCalculationUsingFshort(allFMdata, params, term2, N, varStr)

% probability fit switch
GEVfitSwitch = 1;

% associate a probability of occurrence for the actual loads data
y_peak_sorted = sort(allFMdata);
P_eF_d = (N:-1:1)./N;

% calculate a fit from the loads data to extapolate for the 50-year load
F = linspace(0,1.75*max(allFMdata),10000);
[y,x] = hist(allFMdata,50);

area = trapz(x,y);
y = y./area;

if GEVfitSwitch
    g = gevfit(allFMdata);
    p = gevcdf(F,g(1),g(2),g(3));
else
    g = wblfit(allFMdata);
    p = wblcdf(F,g(1),g(2));
end

term1 = 1-p;
for ws2 = 1:length(params.ws)
    Integrand(ws2,:) = term1.*term2(ws2);
end
for i = 1:length(F)
    P_eF(i) = trapz(params.ws,Integrand(:,i));
end
P_eF = fliplr(P_eF);
F = fliplr(F);

if GEVfitSwitch
    pp = gevpdf(F,g(1),g(2),g(3));
else
    pp = wblpdf(F,g(1),g(2));
end

% figure
% plot(x,y), hold on
% plot(F,pp)
% set(gcf,'color','w');
% grid on
% xlabel(varStr)
% ylabel('Probability')
% grid on

figure
semilogy(F,P_eF); hold on
plot(y_peak_sorted,P_eF_d,'.'); hold on
load_50 = [0 max(F)];
e_load_50 = 3.8E-7.*ones(1,2);
plot(load_50,e_load_50,'r-')
set(gcf,'color','w');
xlabel(varStr)
ylabel('Probability of Exceedance')
grid on
ylim([10^-7 10^0])

try
    % find the 50-year estimate
    i_interp = find(P_eF>=1E-8 & P_eF<=1E-1);
    % calculate the 50 year estimated load
    F_50yearEstimate = interp1(P_eF(i_interp),F(i_interp),e_load_50(1),'pchip');
catch
    % %     keyboard;
end

end