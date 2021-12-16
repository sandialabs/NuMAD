function [output,outputAllFiles1p2,Dminer]=runIECDLC1p1_1p2(params,output,DLCoptions)
global fastPath
global adamsPath  
global turbsimPath

CaseName='IECDLC1p1NTM';


if params.simulate
    runturbsim=1;runfast=1;runcrunch=1;saverccdata=1;
else
    % re-run Crunch and save RCCDATA when EI and/or c change for calculated strains
    runturbsim=0;runfast=0;runcrunch=1;saverccdata=1;
    % do not need to re-run Crunch and save RCCDATA when EI and c do not change
% %     runturbsim=0;runfast=0;runcrunch=0;saverccdata=0;
end

switch params.Class
    case 1
        avgws=0.2*50; % m/s, average wind speed of IEC Class I site (Vref=50m/s); IEC Section 6.3.1.1 Eqn (9)
    case 2
        avgws=0.2*42.5; % m/s, average wind speed of IEC Class II site (Vref=42.5m/s); IEC Section 6.3.1.1 Eqn (9)
    case 3
        avgws=0.2*37.5; % m/s, average wind speed of IEC Class III site (Vref=37.5m/s); IEC Section 6.3.1.1 Eqn (9)
end

% Channel labels for crunch RF analysis - Moments only
labels = ['RootMxb1' params.bladeGageLabels_MLx...
        'RootMyb1' params.bladeGageLabels_MLy];
if 0 % add turbine LSS moments into the crunch analysis
    labelsLSS = {'LSShftMxa' 'LSSTipMys' 'LSSTipMzs' 'LSShftFxa' 'LSShftFys' 'LSShftFzs'};
    params.resultantMoments = 0;
    labels = labelsLSS;
end

if params.combinedStrain
    labelsEPS = ['RootMxb1' params.bladeGageLabels_MLx...
        'RootMyb1' params.bladeGageLabels_MLy ...
        'RootFzb1' params.bladeGageLabels_FLz];
    labels = labelsEPS;
end

% ====================== MAIN PROCESS - SIMULATIONS =======================
% create a counter for parallel runs
for ii = 1:length(params.ws)
    ctr.w((ii-1)*params.numSeeds+1:ii*params.numSeeds) = params.ws(ii);
    ctr.s((ii-1)*params.numSeeds+1:ii*params.numSeeds) = 1:params.numSeeds;
end
ctr.w1 = ctr.w; ctr.s1 = ctr.s;
yawMat = repmat(params.yaw, params.numSeeds*length(params.ws), 1);
ctr.y = reshape(yawMat, 1, numel(yawMat));
ctr.w = repmat(ctr.w, 1, length(params.yaw));
ctr.s = repmat(ctr.s, 1, length(params.yaw));
% loop through "wind direction" as well
dirMat = repmat(params.wd, numel(yawMat), 1);
ctr.d = reshape(dirMat, 1, numel(dirMat));
ctr.y = repmat(ctr.y, 1, length(params.wd));
ctr.w = repmat(ctr.w, 1, length(params.wd));
ctr.s = repmat(ctr.s, 1, length(params.wd));


if runturbsim
    parfor ii = 1:length(ctr.w1)
        fst=readFastMain(['IEC_' params.fstfn '.fst']);
        tsim=readTurbSim([params.parDir 'init\wind.inp']);
        hm=pwd;
        cd('wind')
        % generate turbsim file
        tsim.RandSeed1=params.seeds(ctr.s1(ii));%ble %#ok<PFBNS>
        tsim.URef=ctr.w1(ii);
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
        writeTurbSim(tsim,['turb_DLC1p1_' num2str(ctr.w1(ii)) 'mps_seed' num2str(ctr.s1(ii)) '.inp']);
        dos([turbsimPath ' turb_DLC1p1_' num2str(ctr.w1(ii)) 'mps_seed' num2str(ctr.s1(ii)) '.inp'],'-echo');
        cd(hm)
    end
end


parfor ii=1:numel(ctr.w)
    % ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % set up a format so that files can be opened by parallel runs.
    t = getCurrentTask;
    if isempty(t)
        id = ''; % first iteration of genetic algorithm in generation
    else
        id = ['_parallelWorker_' num2str(t.ID)]; % remaining population are sent to workers in parallel
        if runfast
            % copy the FAST and AD input files for parallel runs
% %             disp([params.fstfn '.fst'])
% %             copyfile([params.fstfn '.fst'], [params.fstfn id '.fst'])
% %             fst=readFastMain(['IEC_' params.fstfn id '.fst']);
% %             copyfile(strrep(fst.ADFile,'"',''), [fst.ADFile(2:end-8) id '_AD.ipt'])
        end
    end
    % ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    if runfast
        % prepare and perform aeroelastic simulation
        fst=readFastMain(['IEC_' params.fstfn '.fst']);%fst=readFastMain(['IEC_' params.fstfn id '.fst']);
        ad=readFastAD([fst.ADFile(2:end-8) '_AD.ipt']);%ad=readFastAD([fst.ADFile(2:end-8) id '_AD.ipt']);
        ad.WindFile=['"wind\turb_DLC1p1_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '.wnd"'];
        fst.SimCtrl.TMax=params.SimTime;
        fst.ADFile=['"' CaseName '_' params.fstfn id '_AD.ipt' '"'];% ble: change the AD filename
        if ~strcmp(params.fastsim,'fast simulink')
            fst=runIEC_setIC(fst,params.operatingPoints(2),'');
        end
        fst.Init.NacYaw = ctr.y(ii);
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
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '_dir' num2str(ctr.d(ii)) '.out'];
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
                YawPosBias = ctr.d(ii)
                runFASTusingSimulink(fastFilenameSFunc, params.simulinkModel,YawPosBias)
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '_dir' num2str(ctr.d(ii)) '.out'];
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
                writeFastMain(fst,[thisFastName '_4ADAMS.fst']);
                dos([fastPath ' ' thisFastName '_4ADAMS.fst']);
                dos([adamsPath ' ' thisFastName '_4ADAMS_ADAMS.acf']);
                disp('Starting ADAMS model from Matlab.....')
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '_dir' num2str(ctr.d(ii)) '.out'];
                try movefile([thisFastName '_4ADAMS_ADAMS.plt'], ['out/' outname])
                catch err; warning(err.message); end
                disp('AD and ADAMS output files saved under new names.')
                disp(' ')
                disp('          *********** Done **************')
                disp(' ')
        end
    end
end

% create the blade structural constants for strain calculations
fst=readFastMain(['IEC_' params.fstfn '.fst']);
bld=readFastBlade(strrep(fst.BldFile{1},'"',''));
ad=readFastAD(strrep(fst.ADFile,'"',''));

% store the blade gage locations, [root, spanwise gages]
spans=[fst.TurbConf.HubRad ad.RNodes(fst.Out.BldGagNd)'];
spans=spans'-spans(1)';
BlLength=fst.TurbConf.TipRad-fst.TurbConf.HubRad;

% load the neutral axis information
load carray
cs(:,1)=interp1(carray(:,2)*BlLength,carray(:,5),spans);  % edge bending c's
cs(:,2)=interp1(carray(:,2)*BlLength,carray(:,6),spans);  % flap bending c's

EI_flap=interp1(bld.prop.BlFract*BlLength,bld.prop.FlpStff,spans);
EI_edge=interp1(bld.prop.BlFract*BlLength,bld.prop.EdgStff,spans);
EIs=[EI_edge EI_flap];
EA_normal=interp1(bld.prop.BlFract*BlLength,bld.prop.EAStff,spans);


% perform the design load case 1.1 50-year extrapolation
if any(contains(DLCoptions,'1.1')) || any(contains(DLCoptions,'all'))
    % EMA: Set the Weibull distribution parameters from the wind speed
    % wSWblP = [16,2.5];
    % wSWblP = [5.5,4.0];
        
    weibullParams.shape = 2; % scale factor of two is a Rayleigh distribution
    weibullParams.scale = avgws / gamma(1+1/weibullParams.shape); % most probable wind speed formulation used for calculation of scale factor
    
    [output,outputAllFiles1p2] = IECDLC1p1_extrapolation(ctr,CaseName,params,weibullParams,cs,EIs,EA_normal,output);
end

% perform the design load case 1.2 fatigue life calculation
Dminer = 999;
if any(contains(DLCoptions,'1.2')) %|| any(contains(DLCoptions,'all'))
    [Dminer] = IECDLC1p2_fatigue(CaseName,params,ctr,labels,avgws,runcrunch,saverccdata,fst,cs,EIs,EA_normal,output);
end

end

%%
function [output,outputAllFiles] = IECDLC1p1_extrapolation(ctr,CaseName,params,weibullParams,cs,EIs,EA_normal,output)
% =========================================================================
% Perform 50-year load extrapolation according to IEC design load case 1.1
% =========================================================================

% ==================== 50-YEAR EXTRAPOLATION CALCULATIONS =================

% Find peaks in the files according to the standard; "Select the largest
% value between successive upcrossings of the mean plus 1.4 times the
% standard deviation of the load time series
output2 = cell(1,length(ctr.w));
parfor ii = 1:length(ctr.w)
    % identify the peaks from each FAST simulation output file
    outName=['out/' CaseName '_yaw' num2str(ctr.y(ii)) '_' num2str(ctr.w(ii)) 'mps_seed' num2str(ctr.s(ii)) '_dir' num2str(ctr.d(ii)) '.out'];
    
    % initialize the output2 structure which will store the peaks for each
    % simulation output file (outName)
    output2{ii}.MaxBladeFlapBendingMoment=[]; % set up output structure
    output2{ii}.MaxFlapStrain=[];
    % call the functions to set up the calculated channels and 
    output2{ii} = sortExtremes(params,output2{ii},CaseName,{outName},cs,EIs,EA_normal);
end


% determine the unique wind speed simulation peaks
[ws,wsIndex] = unique(ctr.w);
wsIndexRange = [wsIndex [wsIndex(2:end)-1; length(ctr.w)]];


% combined the peaks into bins (currently only for wind speed, could be
% added for yaw, wind direction, etc.)
names = fieldnames(output2{1});
parfor ii = 1:length(names)
    nChannels = sum(contains(fieldnames(output2{1}.(names{ii})),'peaks','IgnoreCase',1));
    for jj = 1:size(wsIndexRange,1)
        data = [output2{wsIndexRange(jj,1):wsIndexRange(jj,2)}];
        data2 = [data.(names{ii})];
        if ctr.w(wsIndexRange(jj,1)) ~= ctr.w(wsIndexRange(jj,2))
            error('wind speed not the same')
        end
        peaks{ii}(jj).variable = names{ii};
        peaks{ii}(jj).windSpd = ws(jj);
        for kk = 1:nChannels
% %             output2{wsIndexRange(jj,1)}.(names{ii}).(['peaks_' num2str(kk)])
            peaks{ii}(jj).(['peaks_' num2str(kk)]) = [data2.(['peaks_' num2str(kk)])];   
        end
    end
end

% Use the identified peaks (max and min) within the time series output
% % files to extrapolate to the 50-year probability 
% Get the 50 year extrapolated value for max/min of each output channel
simTime_sec = params.SimTime - params.delay;
totalSimTime_wspd = simTime_sec*params.numSeeds*length(params.yaw)*length(params.wd);
fiftyYears_sec = 50*365.25*24*60*60;
timeRatio = totalSimTime_wspd / fiftyYears_sec; % NOTE: assumes even distribution of yaw and wind direction variables (e.g., probably not right)

maxPeakFit=cell(1,length(names)); extremeFiftyYear=nan(1,length(names));
parfor ii = 1:length(names)
    maxFiftyYear = [];
    nChannels = sum(contains(fieldnames(peaks{ii}(1)),'peaks','IgnoreCase',1));
    isMinimum = contains(peaks{ii}(1).variable,'min','IgnoreCase',1);
    if isMinimum
        gains = -1;
    else
        gains = 1;
    end    
    for cc = 1:nChannels
        peakChan = ['peaks_' num2str(cc)];
        for jj = 1:length(peaks{ii})
            % fit the labels data channel for each wind speed            
            peaksMax = gains .* peaks{ii}(jj).(peakChan);
            peaksMax_mean = nanmean(peaksMax);
            peaksMax_stdev = nanstd(peaksMax);
% % 			% fit using the generalized extreme value distribution for each wind speed
% %             fitParams = gevfit(peaksMax);    
% %             maxPeakFit{ii}(jj,1:3) = fitParams;
% %             maxPeakFit{ii}(jj,4) = length(peaksMax); 
            % fit using a Gaussian distribution for each wind speed
            maxPeakFit{ii}(jj,1:2) = [peaksMax_mean,peaksMax_stdev];
            maxPeakFit{ii}(jj,4) = length(peaksMax);
            if(any(isnan(peaksMax)))
                warning(['peak data has NaN values ' peaks{ii}(jj).variable])
            end
        end
        % calculate the 50 year value for the data channel in labels
        tempMax = get50YearLoad(maxPeakFit{ii},weibullParams,params,timeRatio);
        maxFiftyYear = max([maxFiftyYear tempMax]);
    end
    extremeFiftyYear(ii) = gains .* maxFiftyYear;
end

% assign the 50 year extrapolation value to the final output structure
for ii = 1:length(names)    
    output_new.(names{ii}).data = extremeFiftyYear(ii);
    output_new.(names{ii}).Chan = output2{1}.(names{ii}).Chan;
    output_new.(names{ii}).File = 'Fifty year extrapolation';
    output_new.(names{ii}).time = fiftyYears_sec;
    output_new.(names{ii}).Name = output2{1}.(names{ii}).Name;
end
% add DLC1.1 50 year extrapolation extreme values to the output structure
if size(fieldnames(output),1) == size(fieldnames(output_new),1)
    output = [output; output_new];
else
    output = output_new;
end


% =========================================================================
% NOTE: for now save max/min from DLC1.1 time series in addition to fifty
% year extrapolation for comparison
% =========================================================================

% Find and sort extreme values
outName={};
for yaw = params.yaw
    for w=params.ws
        for s=1:params.numSeeds
            for d=params.wd
                outName{end+1}=['out\' CaseName '_yaw' num2str(yaw) '_' ...
                    num2str(w) 'mps_seed' num2str(s) '_dir' num2str(d) '.out'];
            end
        end
    end
end
% initialize the output2 variable
output3 = cell(1,length(outName));
for cc = 1:length(outName)
%     output2{cc} = output;
    output3{cc}.MaxBladeFlapBendingMoment=[]; % set up output structure
    output3{cc}.MaxFlapStrain=[];
end
% determine the maximum/minimum loads for each output file
parfor cc = 1:length(outName)
    output3{cc} = sortExtremes(params,output3{cc},'IECDLC1p2NTM',outName(cc),cs,EIs,EA_normal);
end

% create the final max/min/ampl values from each of the simulations
names = fieldnames(output3{1});
for ii = 1:length(names) % loop through the channel variable
    % save the data values in a temporary vector
    for jj = 1:length(output3)        
        tempVec(jj) = output3{jj}.(names{ii}).data;
    end
    % determine the index of the channel max/min/ampl value
    if ~isempty(strfind(names{ii},'Min'))
        % Compare output2 channels for minimum value
        [~,index] = min(tempVec);
    else
        % Compare output2 channels for maximum value (includes 'Ampl')
        [~,index] = max(tempVec);
    end
    output_new.(names{ii}) = output3{index}.(names{ii});
end
% add DLC1.2 extreme values to the output structure (for comparison with DLC 1.1)
if size(fieldnames(output),1) == size(fieldnames(output_new),1)
    output = [output; output_new];
else
    output = output_new;
end

% save output2 results in a 2x2 matrix that is more intuitive to search
ct=1; outputAllFiles = cell(length(params.yaw),length(params.ws));
for ii=1:length(params.yaw)
    for jj=1:length(params.ws)
        seedIndex = (ct-1)*params.numSeeds+(1:params.numSeeds);
        outputAllFiles{ii,jj} = [output3{seedIndex}];
        ct = ct+1;
    end
end

end

function F_50year = get50YearLoad(peakFit,weibullParams,params,timeRatio)
    % 
    X1 = 0.0;
    Y1 = probEx(X1,peakFit,weibullParams,params) - timeRatio;
    it = 0;
    while(Y1 <= 0 && it < 50)
        X1 = X1 - 2.0^it;
        Y1 = probEx(X1,peakFit,weibullParams,params) - timeRatio;
        it = it + 1;
    end
    if(it == 50)
        error('Couldnt find suitable X1')
    end
    X2 = 1.0;
    Y2 = probEx(X2,peakFit,weibullParams,params) - timeRatio;
    it = 0;
    while(Y2 >= 0 && it < 50)
        X2 = X2 + 2.0^it;
        Y2 = probEx(X2,peakFit,weibullParams,params) - timeRatio;
        it = it + 1;
    end
    if(it == 50)
        error('Couldnt find suitable Y2')
    end
    Xmid = 0.5*(X1 + X2);
    Ymid = probEx(Xmid,peakFit,weibullParams,params) - timeRatio;
    it = 0;
    while(abs(Ymid) > 0.0001*timeRatio && it < 50)
        if(Ymid > 0)
            X1 = Xmid;
        else
            X2 = Xmid;
        end
        Xmid = 0.5*(X1 + X2);
        Ymid = probEx(Xmid,peakFit,weibullParams,params) - timeRatio;
        it = it + 1;
    end
    if(it >= 50)
        warning('Root not found')
    end
    F_50year = Xmid;
end

function P_e = probEx(F,peakFit,weibullParams,params)
    Nws = length(params.ws);
    integrand = zeros(1,Nws);
    for i = 1:Nws
        %PexV = 1.0 - gevcdf(F,peakFit(i,1),peakFit(i,2),peakFit(i,3))^peakFit(i,4);
        PexV = 1.0 - cdf('Normal',F,peakFit(i,1),peakFit(i,2))^peakFit(i,4);
        %PexV = 1.0 - cdf('GeneralizedExtremeValue',F,peakFit(i,1),peakFit(i,2),peakFit(i,3))^peakFit(i,4);
        integrand(i) = PexV*wblpdf(params.ws(i),weibullParams.scale,weibullParams.shape);
    end
    P_e = trapz(params.ws,integrand);
end

%%
function [Dminer] = IECDLC1p2_fatigue(CaseName,params,ctr,labels,avgws,runcrunch,saverccdata,fst,cs,EIs,EA_normal,output)
% =========================================================================
% Perform fatigue calculations according to the IEC design load case 1.2
% =========================================================================
global crunchPath

% ble.  runcrunch is only operated on the wind speed counter, once all of
% the seeds have been simulated.
labelsCalcChanNames = {};
if runcrunch % Setup and run Crunch for rainflow cycle counting
    mkdir('rcc'); % ensure that the rcc/ directory exists
    parfor w=1:length(params.ws)
        % read in one example FAST output file
        fastoutname=[CaseName '_yaw' num2str(params.yaw(1)) '_' num2str(params.ws(1)) 'mps_seed1_dir' num2str(ctr.d(1))];  %temps
        fastout=loadFASTOutData(['out\' fastoutname '.out']);
        cru=readCrunch([params.parDir 'init\example.cru']);
        % Prepare crunch information
        % fill in aggregate file name
        cru.JobOpt.AggRoot=[CaseName '_' num2str(params.ws(w)) 'mps'];
        cru.IptDataLay.CTRow=0;
        cru.IptDataLay.CURow=0;
        cru.IptDataLay.FDRow=9;
        cru.IptDataLay.TStartTEnd=[params.delay,params.SimTime];
        % ble: this change is needed for simulink operation because the
        % time saved is not [30, 630] but is [30.018, 629.933]
        if strcmp(params.fastsim, 'fast simulink')
            cru.IptDataLay.TStartTEnd=[params.delay,params.SimTime-0.1];
        end
        timeChannelOffset = 1; % adding a time channel
        % fill in Channel information
        cru.ChanInfo.NumInCols=length(fastout.list);
        cru.ChanInfo.NumCols=length(labels)+timeChannelOffset;   % Time variable is the first channel
        
        % find columns for desired channels        
        cru.ChanInfo.OrigChan(1,1)=strmatch('Time',fastout.list);
        cru.ChanInfo.ChanTitle{1,1}=['"' fastout.list{strmatch('Time',fastout.list)} '"'];
        cru.ChanInfo.ChanUnits{1,1}=['"' fastout.units{strmatch('Time',fastout.list)} '"'];
        cru.ChanInfo.Scale(1,1)=1;
        cru.ChanInfo.Offset(1,1)=0;
        for i=1:length(labels)
            labels(i)
            cru.ChanInfo.OrigChan(i+1,1)=strmatch(labels{i},fastout.list);
            cru.ChanInfo.ChanTitle{i+1,1}=['"' fastout.list{cru.ChanInfo.OrigChan(i+1,1)} '"'];
            cru.ChanInfo.ChanUnits{i+1,1}=['"' fastout.units{cru.ChanInfo.OrigChan(i+1,1)} '"'];
            cru.ChanInfo.Scale(i+1,1)=1;
            cru.ChanInfo.Offset(i+1,1)=0;
        end
        cru.TimeCol=1;
        % fill in NumRFCols
        cru.NumRFCols=length(labels);
        % fill in column ranges
        cru.ColNums=(2:cru.ChanInfo.NumCols)';
        cru.HalfCycMult=ones(cru.NumRFCols,1)*0.5;
        % produce list of files
        cru.NumFiles=params.numSeeds*length(params.yaw);
        ctrFiles = 1;
        for yaw = params.yaw
            for s=1:params.numSeeds
                for d = params.wd
                    cru.InFiles{ctrFiles}=sprintf('"%s_yaw%i_%imps_seed%i_dir%i.out"',CaseName,yaw,params.ws(w),s,d);
                    ctrFiles = ctrFiles+1;
                end
            end
        end
        cru.EEGrps=[];
        
        % Add calculated channels for the rainflow cycle counting
%         labelsCalcChan = labels; 
        cru.CalcChan.ColTitle = {};
        cru.CalcChan.Units = {}; cru.CalcChan.Equation = {};
        
        % =========== calculations for station angular rotations ==========
        if params.momentMaxRotation <= 90 
            thetaMomentRotation = 0:params.momentMaxRotation:180; % [deg]
            thetaMomentRotation(thetaMomentRotation==180)=[]; % remove repetitive 180 deg rotation
            for ii = 1:(fst.Out.NBlGages+1)
                % perform the resultant calculations for x labels
                if isempty(strfind(labels{ii},'yb'))
                    xChan = labels{ii};
                    yChan = strrep(labels{ii},'xb','yb');
                    xChanCol = find(strcmp(cru.ChanInfo.ChanTitle,['"' xChan '"']));
                    cru.ChanInfo.ChanTitle{xChanCol}
                    yChanCol = find(strcmp(cru.ChanInfo.ChanTitle,['"' yChan '"']));
                    cru.ChanInfo.ChanTitle{yChanCol}                   
                    
                    if contains(cru.ChanInfo.ChanTitle(yChanCol),'Spn')
                        structTwist= params.bladeGageCoordinateRotation(ii-1);
                    else
                        structTwist=0;  %Root already correct
                    end
                    
                    for dd = thetaMomentRotation
                        cru.CalcChan.ColTitle{end+1,1} = ['"' xChan(1:strfind(xChan,'x')-1) 'r' num2str(dd) '"'];
                        cru.CalcChan.Units{end+1,1} = '"kN-m"';
                        cru.CalcChan.Equation{end+1,1} = ['"C' num2str(xChanCol) '*COSD(' num2str(dd+structTwist) ')+C' ...
                            num2str(yChanCol) '*SIND(' num2str(dd+structTwist) ')"'];
                    end
                end
            end
            cru.CalcChan.NumCChan = length(cru.CalcChan.ColTitle);            
            % add the calculated channels to the labels file            
            labelsCalcChan = [labels strrep(cru.CalcChan.ColTitle,'"','')'];
        end
        
        % ========== calculations for combined bending/normal strain ======       
        if params.combinedStrain
            gageCtr = 0;
            for ii = 1:length(labels)
                % perform the resultant calculations for Myb labels
                if contains(labels{ii},'yb')
                    disp(labels{ii})
                    gageCtr = gageCtr+1;
                    myChan = labels{ii};
                    mxChan = strrep(labels{ii},'yb','xb');
                    fzChan = strrep(strrep(labels{ii},'yb','zb'),'M','F');
                    myChanCol = find(strcmp(labels,myChan))+timeChannelOffset;
                    mxChanCol = find(strcmp(labels,mxChan))+timeChannelOffset;
                    fzChanCol = find(strcmp(labels,fzChan))+timeChannelOffset;
                    
                    % define the constants to convert load to strain
                    ea_const = EA_normal(gageCtr);
                    cx_const = cs(gageCtr,2);
                    eiflap_const = EIs(gageCtr,2);
                    % conversions
                    Nm_kNm = 1000; 
                    micro = 1e6; % calculate in micro-strain for precision
                    
                    % high pressure side strain of center of spar cap
                    cru.CalcChan.ColTitle{end+1,1} = ['"' myChan(1:strfind(myChan,'M')-1) 'epsHP"'];
                    cru.CalcChan.Units{end+1,1} = '"micro-strain"';
                    cru.CalcChan.Equation{end+1,1} = ['"(C' num2str(fzChanCol) '*' num2str(Nm_kNm) '/' num2str(ea_const) '+C' ...
                        num2str(myChanCol) '*' num2str(Nm_kNm) '*' num2str(cx_const) '/' num2str(eiflap_const) ')*' num2str(micro) '"'];
                    % low pressure side strain of center of spar cap
                    cru.CalcChan.ColTitle{end+1,1} = ['"' myChan(1:strfind(myChan,'M')-1) 'epsLP"'];
                    cru.CalcChan.Units{end+1,1} = '"micro-strain"';
                    cru.CalcChan.Equation{end+1,1} = ['"(C' num2str(fzChanCol) '*' num2str(Nm_kNm) '/' num2str(ea_const) '-C' ...
                        num2str(myChanCol) '*' num2str(Nm_kNm) '*' num2str(cx_const) '/' num2str(eiflap_const) ')*' num2str(micro) '"'];
                                        
                    % TEST -- match fatigue life from stress calculations
                    cru.CalcChan.ColTitle{end+1,1} = ['"' myChan(1:strfind(myChan,'M')-1) 'epsTST"'];
                    cru.CalcChan.Units{end+1,1} = '"micro-strain"'; % use micro-strain for precision concerns
                    cru.CalcChan.Equation{end+1,1} = ['"(C' num2str(myChanCol) '*' num2str(Nm_kNm) ...
                        '*' num2str(cx_const) '/' num2str(eiflap_const) ')*' num2str(micro) '"'];
                end
            end
            cru.CalcChan.NumCChan = length(cru.CalcChan.ColTitle);
            % add the calculated channels to the labels file
            labelsCalcChan = [labels strrep(cru.CalcChan.ColTitle,'"','')'];
        end
        
        cru.NumRFCols=length(labelsCalcChan);
        % fill in column ranges
        cru.ColNums=[cru.ColNums; cru.ColNums(end)+[1:cru.CalcChan.NumCChan]'];
        cru.HalfCycMult=ones(cru.NumRFCols,1)*0.5;
        
        labelsCalcChanNames{w}.cell = labelsCalcChan;
        
        % write new crunch file
        hm=pwd;
        cd('out')
        writeCrunch(cru,[cru.JobOpt.AggRoot '.cru']);
        disp('...Crunch file generated')
        
        % run crunch
        dos([crunchPath ' ' cru.JobOpt.AggRoot '.cru']);
        cd(hm);
        
        % Move rcc files from out/ to rcc/
        for jj=1:length(labelsCalcChan)
            filename=[CaseName '_' num2str(params.ws(w)) 'mps_' labelsCalcChan{jj} '.rcc'];
            copyfile(['out/' filename],['rcc/' filename]);
            delete(['out/' filename]);
            disp(['out/' filename ' has been moved to rcc/' filename]);
        end        
        disp('          .....Finished moving *.rcc files.')
    end
end

% save the raw cycle count (RCC) data
if saverccdata
    for w=1:length(params.ws)
        % read in the channel variable names from the label cell
        labelsCalcChan = labelsCalcChanNames{w}.cell;
        for jj=1:length(labelsCalcChan)
            filename=[CaseName '_' num2str(params.ws(w)) 'mps_' labelsCalcChan{jj} '.rcc'];
            a=txt2mat(['rcc/' filename]);
            if contains(labelsCalcChan{jj},'eps')
                % strain are calculated in micro-strain for precision
                rccdata{jj,w}.cycles=a(:,1);
                rccdata{jj,w}.ranges=a(:,2)./1e6;  % crunch saves range, not amplitude [strain]
                rccdata{jj,w}.amplitudes=a(:,2)./1e6/2;  % crunch saves range, not amplitude [strain]
                rccdata{jj,w}.means=a(:,3)./1e6;         % [Nm]
                rccdata{jj,w}.windspeed=params.ws(w);
                rccdata{jj,w}.label=labelsCalcChan{jj};
            else % convert forces and moments from FAST to SI base units (e.g., N, N-m)
                rccdata{jj,w}.cycles=a(:,1);
                rccdata{jj,w}.ranges=a(:,2).*1000;  % crunch saves range, not amplitude [N-m]
                rccdata{jj,w}.amplitudes=a(:,2).*1000/2;  % crunch saves range, not amplitude [N-m]
                rccdata{jj,w}.means=a(:,3).*1000;         % [Nm]
                rccdata{jj,w}.windspeed=params.ws(w);
                rccdata{jj,w}.label=labelsCalcChan{jj};
            end
        end
    end
    % save the raw cycle counts in rccdata
    if runcrunch
        save rccdata rccdata -v7.3;
    end
end

% set up and run MFatigue
fprintf('Total fatigue safety factor: %f',params.sf_fat)
params.simtime=params.numSeeds*(fst.SimCtrl.TMax-params.delay); % simulated and rainflow counted time, seconds

% ble: fix the MFatigue calculations to determine the correct span location




%% Matlab file 'rccdata.mat' is created by running the aeroelastic
% simulation scripts followed by Crunch analysis

[wt,rccdata]=getWindSpeedDistribution(avgws);

params.fatigueCriterion = 'Shifted Goodman';%'Goodman';%  % fatigue failure criterion
params.fatigueStress = 'Equivalent';%'Amplitude Only';%  % stress to use for fatigue failure ('Equivalent','Amplitude Only',...)

Dminer=MFatigue(wt,rccdata,cs,EIs,params);

% % Find and sort extreme values
% outName={};
% for yaw = params.yaw
%     for w=params.ws
%         for s=1:params.numSeeds
%             for d=params.wd
%                 outName{end+1}=['out\' CaseName '_yaw' num2str(yaw) '_' ...
%                     num2str(w) 'mps_seed' num2str(s) '_dir' num2str(d) '.out'];
%             end
%         end
%     end
% end
% 
% % initialize the output2 variable
% output2 = cell(1,length(outName));
% for cc = 1:length(outName)
% %     output2{cc} = output;
%     output2{cc}.MaxBladeFlapBendingMoment=[]; % set up output structure
%     output2{cc}.MaxFlapStrain=[];
% end
% % determine the maximum/minimum loads for each output file
% parfor cc = 1:length(outName)
%     output2{cc} = sortExtremes(params,output2{cc},'IECDLC1p2NTM',outName(cc),cs,EIs,EA_normal);
% end
% % create the final max/min/ampl values from each of the simulations
% names = fieldnames(output2{1});
% for ii = 1:length(names) % loop through the channel variable
%     % save the data values in a temporary vector
%     for jj = 1:length(output2)        
%         tempVec(jj) = output2{jj}.(names{ii}).data;
%     end
%     % determine the index of the channel max/min/ampl value
%     if ~isempty(strfind(names{ii},'Min'))
%         % Compare output2 channels for minimum value
%         [~,index] = min(tempVec);
%     else
%         % Compare output2 channels for maximum value (includes 'Ampl')
%         [~,index] = max(tempVec);
%     end
%     output_new.(names{ii}) = output2{index}.(names{ii});
% end
% if size(fieldnames(output),1) == size(fieldnames(output_new),1)
%     output = [output; output_new];
% else
%     output = output_new;
% end
% 
% % save output2 results in a 2x2 matrix that is more intuitive to search
% ct=1; outputAllFiles = cell(length(params.yaw),length(params.ws));
% for ii=1:length(params.yaw)
%     for jj=1:length(params.ws)
%         seedIndex = (ct-1)*params.numSeeds+(1:params.numSeeds);
%         outputAllFiles{ii,jj} = [output2{seedIndex}];
%         ct = ct+1;
%     end
% end


end