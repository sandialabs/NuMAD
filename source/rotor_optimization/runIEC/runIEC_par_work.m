function output=runIEC_par_work(DLCoptions,simFlag,parFlag)
% options.  String array specifying which DLC's to run.  For example,
%    DLCoptions = '1.2 1.4 6.1 2.1';
% or
%    DLCoptions = 'all'
% simFlag = 0/1; (1) call FAST and perform simulations or (0) just process existing data
%

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if exist('parFlag','var')
	if parFlag == 1
		run('..\runIEC_ipt.m');
		copyfile(['..\' params.fstfn '.fst'], [params.fstfn '.fst'])
		fst=readFastMain([params.fstfn '.fst']);
		copyfile(['..\' strrep(fst.ADFile,'"','')], strrep(fst.ADFile,'"',''))
		copyfile('..\pitch.ipt', 'pitch.ipt')        
		params.parDir = '..\';
	else
		runIEC_ipt
        fst=readFastMain([params.fstfn '.fst']); 
		params.parDir = '';
%         delete('wind\*')
	end
else
	runIEC_ipt
	fst=readFastMain([params.fstfn '.fst']);   
	params.parDir = '';
end

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
params.fullLoads = 1;                   % perform full loads analysis?
params.fullLoadsResultantMoments = 1;   % perform resultant moment calculation and storage?
gageSetCase = 'set1';%'set2';%
params.resultantMoments = 0;            % perform fatigue calculations on resultant moment stations?
params.simulate=simFlag;
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>



fst.Out.TStart=params.delay;
% read AD input data
ad=readFastAD(strrep(fst.ADFile,'"',''));

% set up blade gages and locations
fst.Out.BldGagNd=params.BldGagNd;  % Units: -
fst.Out.NBlGages=length(fst.Out.BldGagNd);  % Units: -
for ss = 1:fst.Out.NBlGages
    params.bladeGageLabels_MLx{ss} = ['Spn' num2str(fst.Out.BldGagNd(ss)) 'MLxb1']; %#ok<*AGROW>
    params.bladeGageLabels_MLy{ss} = ['Spn' num2str(fst.Out.BldGagNd(ss)) 'MLyb1'];
    params.bladeGageLabels_MLz{ss} = '';
    params.bladeGageLabels_FLx{ss} = '';
    params.bladeGageLabels_FLy{ss} = '';
    params.bladeGageLabels_FLz{ss} = '';
end

params.bladeGageLabels = [params.bladeGageLabels_MLx params.bladeGageLabels_MLy];%...
    %params.bladeGageLabels_MLz params.bladeGageLabels_FLx...
    %params.bladeGageLabels_FLy params.bladeGageLabels_FLz];

if params.fullLoadsResultantMoments
    % this script must be run twice for the first 9 gages, and then for the
    % additional gages.
    switch gageSetCase
    case 'set1'
        setGag = 1:9; % *set 1*
    case 'set2'
        setGag = 10:14; % *set 2*
    end
    fst.Out.NBlGages = length(setGag);
    totalBladeGageNumber = 9;   % actual number of gages to be used
    BldGagMaxSpanLoc = 0.95;
    BldGagMinSpanLoc = 0;     % minimum span location of a strain gage (no gage here currently; RootMxb1...)
    drSpan = (BldGagMaxSpanLoc - 0) / totalBladeGageNumber;
    gagSpan = BldGagMinSpanLoc + drSpan : drSpan : BldGagMaxSpanLoc;
    bladeCoordinateSpanStations = (ad.RNodes - fst.TurbConf.HubRad) ./ (fst.TurbConf.TipRad-fst.TurbConf.HubRad);
    fst.Out.BldGagNd = 0; 
    params.bladeGageLabels_MLx = {}; params.bladeGageLabels_MLy = {}; params.bladeGageLabels_MLz = {};
    params.bladeGageLabels_FLx = {}; params.bladeGageLabels_FLy = {}; params.bladeGageLabels_FLz = {};
    params.bladeGageLabels = {};
    for bb = 1%:3 %number of blades
    for ss = 1:length(setGag)%fst.Out.NBlGages
        [~, rGagSpan(ss)] = min(abs(bladeCoordinateSpanStations - gagSpan(setGag(ss))));
        fst.Out.BldGagNd(ss) = rGagSpan(ss);
        if bb == 1
            params.bladeGageLabels_MLx{ss} = ['Spn' num2str(ss) 'MLxb' num2str(bb)]; %#ok<*AGROW>
            params.bladeGageLabels_MLy{ss} = ['Spn' num2str(ss) 'MLyb' num2str(bb)];
            params.bladeGageLabels_MLz{ss} = ['Spn' num2str(ss) 'MLzb' num2str(bb)];
            params.bladeGageLabels_FLx{ss} = ['Spn' num2str(ss) 'FLxb' num2str(bb)];
            params.bladeGageLabels_FLy{ss} = ['Spn' num2str(ss) 'FLyb' num2str(bb)];
            params.bladeGageLabels_FLz{ss} = ['Spn' num2str(ss) 'FLzb' num2str(bb)];
            rowName{ss} = ['gage' num2str(ss)];
        else
            params.bladeGageLabels_MLx{end+1} = ['Spn' num2str(ss) 'MLxb' num2str(bb)]; %#ok<*AGROW>
            params.bladeGageLabels_MLy{end+1} = ['Spn' num2str(ss) 'MLyb' num2str(bb)];
            params.bladeGageLabels_MLz{end+1} = ['Spn' num2str(ss) 'MLzb' num2str(bb)];
            params.bladeGageLabels_FLx{end+1} = ['Spn' num2str(ss) 'FLxb' num2str(bb)];
            params.bladeGageLabels_FLy{end+1} = ['Spn' num2str(ss) 'FLyb' num2str(bb)];
            params.bladeGageLabels_FLz{end+1} = ['Spn' num2str(ss) 'FLzb' num2str(bb)];
        end
    end
    end
    params.bladeGageLabels = [params.bladeGageLabels_MLx params.bladeGageLabels_MLy...
        params.bladeGageLabels_MLz params.bladeGageLabels_FLx...
        params.bladeGageLabels_FLy params.bladeGageLabels_FLz];
end
tbl = array2table([gagSpan(setGag)' bladeCoordinateSpanStations(rGagSpan) rGagSpan']);
tbl.Properties.VariableNames = {'GageSpanLocation_bladeCoordinates' 'ActualSpanLocations_bladeCoordinates' 'ADlocation'};
tbl.Properties.RowNames = rowName;
% disp(tbl)
% disp('Press F5 to confirm the strain gage locations and proceed...')
% keyboard 


if false(params.fullLoads)
    % define the necessary FAST outputs
    fst.OutList=[{'WindVxi','WindVyi','WindVzi',...
        'HorWindV','NacYawErr',...  % add these channels for JCB's yaw controller
        'GenPwr','HSShftTq','HSSBrTq','GenSpeed','RotSpeed','TSR',...
        'TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS',...
        'NcIMUTAxs','NcIMUTAys',...
        'BldPitch2',...
        'OoPDefl1','OoPDefl2','OoPDefl3',...
        'TipClrnc1','TipClrnc2','TipClrnc3',...
        'TipDxb1','TipDxb2','TipDxb3',...
        'RootMxb1','RootMyb1','RootMzb1',...
        'RootMxb2','RootMyb2','RootMzb2',...
        'RootMxb3','RootMyb3','RootMzb3',...
        'RootMxc1','RootMyc1',...
        'RootMxc2','RootMyc2',...
        'RootMxc3','RootMyc3'} ...
        params.bladeGageLabels...
        {'LSShftFxa','LSShftMxa'}];   
else
    % set up analysis for complete loads determination to ensure turbine
    % operation within structure loads envelope
    fst.OutList=[{'WindVxi','WindVyi','WindVzi',...
        'HorWindV','NacYawErr',...  % add these channels for JCB's yaw controller
        'GenPwr','HSShftTq','HSSBrTq','GenSpeed','RotSpeed','TSR',...
        'TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS',...
        'NcIMUTAxs','NcIMUTAys',...
        'BldPitch2',...
        'OoPDefl1','OoPDefl2','OoPDefl3',...
        'TipClrnc1','TipClrnc2','TipClrnc3',...
        'TipDxb1','TipDxb2','TipDxb3'}...
        params.bladeGageLabels...
        ... SOME NEW CHANNELS TO ADD:
        {'RootFxc1','RootFxc2','RootFxc3',...
        'RootFyc1','RootFyc2','RootFyc3',...
        'RootFxb1','RootFxb2','RootFxb3',...
        'RootFyb1','RootFyb2','RootFyb3',...
        'RootFzb1','RootFzb2','RootFzb3',...
        'RootMxc1','RootMxc2','RootMxc3',...
        'RootMyc1','RootMyc2','RootMyc3',...
        'RootMxb1','RootMxb2','RootMxb3',...
        'RootMyb1','RootMyb2','RootMyb3',...
        'RootMzb1','RootMzb2','RootMzb3',...
        'YawBrFxn','YawBrFyn','YawBrFzn',...
        'YawBrFxp','YawBrFyp',...
        'YawBrMxn','YawBrMyn','YawBrMzn',...
        'YawBrMxp','YawBrMyp',...
        'RotThrust','LSShftTq',...
        'LSShftFxa','LSShftFya','LSShftFza',...
        'LSShftFys','LSShftFzs',...
        'LSShftMxa',...
        'LSSTipMya','LSSTipMza',...
        'LSSTipMys','LSSTipMzs',...
        'TwrBsFxt','TwrBsFyt','TwrBsFzt',...
        'TwrBsMxt','TwrBsMyt','TwrBsMzt',...
        }];
end
%CK
% disp(fst.OutList)
% pause(5)

for i=1:length(ad.PrnElm)
    ad.PrnElm{i}='PRINT';
end

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if parFlag == 1
    if ~strcmp(ad.FoilNm{1}(2:4),'..\')
    for ii = 1:ad.NumFoil, ad.FoilNm{ii} = ['"..\' ad.FoilNm{ii}(2:end)]; end
    end
    if ~strcmp(fst.Twr.TwrFile(2:4),'..\')
    fst.Twr.TwrFile = ['"..\' fst.Twr.TwrFile(2:end)];
    end
end
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% Write input files
writeFastMain(fst,[params.fstfn '.fst']);
writeFastAD(ad,strrep(fst.ADFile,'"',''));

% find distances from neutral axes for this blade
if ~exist('carray.mat','file')
    make_c_array
    plot_c_array
else
    disp('Using previously saved neutral axis information, c.')
    load carray
end

% Set random seeds
if ~exist([params.parDir 'seeds.mat'],'file')
    seeds=randi(123456,1,params.numSeeds);
    save([params.parDir 'seeds'],'seeds')
else
    load([params.parDir 'seeds'])
    % ble: added this check - seeds file can exist but not be correct.
    if length(seeds) ~= params.numSeeds
        clear seeds
        seeds=randi(123456,1,params.numSeeds);
        save([params.parDir 'seeds'],'seeds')
    end
    % ble: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
end
params.seeds = seeds; %ble: seeds needs to be saved this way for parallel operation.

output.MaxFlapBendingMoment=[]; % set up output structure
output.MaxRotorFlapBendingMoment=[];
output.MaxFlapStrain=[];
output.MaxEdgeStrain=[];
output.MaxOoPDefl=[];
output.MaxRotorTh=[];
output.MinFlapBendingMoment=[];
output.MinRotorFlapBendingMoment=[];
output.MinFlapStrain=[];
output.MinEdgeStrain=[];
output.MinOoPDefl=[];
output.MinRotorTh=[];

iecwind=makeIECWind(params); % make IECWind files

if ~isempty(strfind(DLCoptions,'sweep')) %ble|| ~isempty(strfind(DLCoptions,'all'))
    output=runIECSweep(params,output); % run wind speed sweep
end
if ~isempty(strfind(DLCoptions,'1.1')) %ble|| ~isempty(strfind(DLCoptions,'all'))
    output=runIECDLC1p1(params,output); % run DLC 1.1 50-year occurrence load
end
if ~isempty(strfind(DLCoptions,'1.2')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output1p2]=runIECDLC1p2(params,output); % run DLC 1.2 Fatigue
    save output1p2 output1p2
end
if ~isempty(strfind(DLCoptions,'1.3')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output1p3]=runIECDLC1p3(params,output); % run DLC 1.3 UTS, ETM
    save output1p3 output1p3
end
if ~isempty(strfind(DLCoptions,'1.4')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output1p4]=runIECDLC1p4_par_work(params,output); % run DLC 1.4 ECD with direction change
    save output1p4 output1p4
end
if ~isempty(strfind(DLCoptions,'1.5')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output1p5]=runIECDLC1p5(params,output); % run DLC 1.5 EWS
    save output1p5 output1p5
end
if ~isempty(strfind(DLCoptions,'6.1')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output6p1]=runIECDLC6p1(params,output); % run DLC 6.1 EWM 50-year parked
    save output6p1 output6p1
end
if ~isempty(strfind(DLCoptions,'6.2')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output6p2]=runIECDLC6p2(params,output); % run DLC 6.2 EWM 50-year parked, no yaw control
    save output6p2 output6p2
end
if ~isempty(strfind(DLCoptions,'6.3')) || ~isempty(strfind(DLCoptions,'all'))
    [output,output6p3]=runIECDLC6p3(params,output); % run DLC 6.3 EWM 1-year parked
    save output6p3 output6p3
end

if ~isempty(strfind(DLCoptions,'loads')) 
    % this script does not run an analysis, but returns loads from a
    % performed analysis of DLC's 1.2 and 1.3
    [results1p2]=runIECDLC1p2_returnLoads(params,output); % read DLC 1.2 loads
    % plot results
    figure
    subplot 211
    histogram(results1p2.Fx_Err_kN,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fx Error [kN]')
    title('DLC 1p2')
    grid on
    subplot 212
    histogram(results1p2.Fy_Err_kN,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fy Error [kN]')
    grid on
    edges = [-100:1:100]
    figure
    subplot 211
    histogram(results1p2.Fx_PercErr,edges,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fx % Error')
    title('DLC 1p2')
    grid on
    subplot 212
    histogram(results1p2.Fy_PercErr,edges,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fy % Error')
    grid on
    % return loads from DLC 1.3 analysis
    [results1p3]=runIECDLC1p3_returnLoads(params,output); % read DLC 1.2 loads
    figure
    subplot 211
    histogram(results1p3.Fx_Err_kN,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fx Error [kN]')
    title('DLC 1p3')
    grid on
    subplot 212
    histogram(results1p3.Fy_Err_kN,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fy Error [kN]')
    grid on
    edges = [-100:1:100]
    figure
    subplot 211
    histogram(results1p3.Fx_PercErr,edges,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fx % Error')
    title('DLC 1p3')
    grid on
    subplot 212
    histogram(results1p3.Fy_PercErr,edges,'DisplayStyle','stairs','Normalization','probability')
    xlabel('Fy % Error')
    grid on
end

% % runLinear(params); % run linearization

% Write results to Excel file

% xlsfn='IECDLC_Results.xlsx';
% writeResult(output,xlsfn);
% save output output

end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function writeResult(output,xlsfn)

names=fieldnames(output);
numDLC = size(output,1);
xlsData=[];

for i=1:length(names)
    xlsData=[xlsData; {'DLC Name',names{i},'Occurring at time (sec)','on Channel','in File'}];
    combinedOut.(names{i}) = [output.(names{i})];
    
    for j=1:numDLC
        xlsData=[xlsData;...
            {getfield(combinedOut,names{i},{j},'Name'),...
            getfield(combinedOut,names{i},{j},'data'),...
            getfield(combinedOut,names{i},{j},'time'),...
            getfield(combinedOut,names{i},{j},'Chan'),...
            getfield(combinedOut,names{i},{j},'File')}];
    end
    %     xlswrite(xlsfn,xlsData,names{i});
end
xlswrite(xlsfn,xlsData,1);

% names=fieldnames(output);
% xlsData=[];
% 
% for i=1:length(names)
%     data=getfield(output,names{i});
%     xlsData=[xlsData; {'DLC Name',names{i},'Occurring at time (sec)','on Channel','in File'}];
%     for j=1:length(data)
%         xlsData=[xlsData;...
%             {getfield(output,names{i},{j},'Name'),...
%             getfield(output,names{i},{j},'data'),...
%             getfield(output,names{i},{j},'time'),...
%             getfield(output,names{i},{j},'Chan'),...
%             getfield(output,names{i},{j},'File')}];
%     end
%     %     xlswrite(xlsfn,xlsData,names{i});
% end
% xlswrite(xlsfn,xlsData,1);

end