function output=runIEC(DLCoptions,simFlag,parFlag)
% This function runs IEC design load cases and creates an
% ``output`` structure
%
% Parameters
% ------------
%     DLCoptions : string
%         String array specifying which design load cases to run, 
%         Options = {'1.2' '1.4' '6.1' '2.1'} or 'all', Default = ?
%     
%     simFlag : int
%         On/Off flag, Options: 1 = call FAST and perform simulations, 
%         0 = process existing data, Default = ?
%     
%     parFlag : int
%         On/Off flag, Options: 0 = ?, 1 = ?, Default = ?
%
% Returns
% ------------
%     output : struct
%         output structure      
%

if exist('parFlag','var')
	if parFlag == 1
        iecInputFile
		copyfile(['..\' params.fstfn '.fst'], [params.fstfn '.fst'])
		fst=readFastMain([params.fstfn '.fst']);
		copyfile(['..\' strrep(fst.ADFile,'"','')], strrep(fst.ADFile,'"',''))
		copyfile('..\pitch.ipt', 'pitch.ipt')
		params.parDir = '..\';
	else
		iecInputFile
        fst=readFastMain([params.fstfn '.fst']); 
		params.parDir = '';
	end
else
	iecInputFile
	fst=readFastMain([params.fstfn '.fst']);   
	params.parDir = '';
end

% Run IECDef methods to check inputs and set average wind speed
params.checkInputs
params.setAvgWindSpeed

% Read AD input data
ad=readFastAD(strrep(fst.ADFile,'"',''));

% Read Blade input data
bld1=readFastBlade(strrep(fst.BldFile{1},'"',''));

% Add Simulink controller path -- removed at end of script so it can change
addpath(genpath(params.simulinkModelFolder));

% Run IECDef method to set simulation flag
params.setSimFlag(simFlag)

    %Moved to IECDef Class
    % params.fullLoads = 1;                   % perform full loads analysis?
    % params.combinedStrain = 1;              % perform fatigue calculations on combined bending and normal strain for spar
    % params.gageSetCase = 'set1';%'set2';%

% Change FAST input file parameters for analysis
fst.Out.TStart=params.delay;

% Change IMU location to match the safety system accelerometer switch location
fst.Out.NcIMUxn = 0.965;
fst.Out.NcIMUyn = 0.00;
fst.Out.NcIMUzn = 0.63;

% Run IECDef method to set up blade gages and locations
params.setGageLabels(fst,ad)
    
    % Moved to IEfDec, replaced with params
    % fst.Out.BldGagNd=params.BldGagNd;  % Units: -
	% fst.Out.NBlGages=length(params.BldGagNd);  % Units: -

% Calculate the coordinate transformation from local blade coordinate
% system used for spanwise gages to the root blade coordinate system
gagStrcTwst = interp1(bld1.prop.BlFract, bld1.prop.StrcTwst, gagSpan);    

% Run IECDef method to set gage coordinate rotation
params.setBladeGageCoordinateRotation(gagStrcTwst);

% Run IECDef method to run full loads or
params.runFullLoads

t = getCurrentTask;
if isempty(t)
    pause(5)
end

for i=1:length(ad.PrnElm)
    ad.PrnElm{i}='PRINT';
end

% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if parFlag == 1 %ble: changing this so the airfoil polars and tower file are in the same directory as FAST input file
%     if ~strcmp(ad.FoilNm{1}(2:4),'..\')
%     for ii = 1:ad.NumFoil, ad.FoilNm{ii} = ['"..\' ad.FoilNm{ii}(2:end)]; end
%     end
%     if ~strcmp(fst.Twr.TwrFile(2:4),'..\')
%     fst.Twr.TwrFile = ['"..\' fst.Twr.TwrFile(2:end)];
%     end
    copyfile(['..\' fst.Twr.TwrFile(2:end-1)], fst.Twr.TwrFile(2:end-1))
end
% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Rewrite the input files with IEC designation to not overwrite previous values
fst.ADFile = ['"IEC_' fst.ADFile(2:end)];
writeFastMain(fst,['IEC_' params.fstfn '.fst']);
writeFastAD(ad,strrep(fst.ADFile,'"',''));

% find distances from neutral axes for this blade
if ~exist('carray.mat','file')
    make_c_array
    plot_c_array
else
    disp('Using previously saved neutral axis information, c.')
    load carray
end

% Run IECDef method to set random seeds
params.setRandomSeeds

% Generate output structure
output.MaxBladeFlapBendingMoment=[]; % set up output structure
output.MinBladeFlapBendingMoment=[];
output.MaxFlapStrain=[];
output.MinFlapStrain=[];
output.MaxEdgeStrain=[];
output.MinEdgeStrain=[];
output.MaxOoPDefl=[];
output.MinOoPDefl=[];
output.MaxRotorTh=[];
output.MinRotorTh=[];

iecwind=makeIECWind(params); % make IECWind files

if any(contains(DLCoptions,'sweep')) %ble|| any(contains(DLCoptions,'all'))
    output=runIECSweep(params,output); % run wind speed sweep
end
if any(contains(DLCoptions,'1.1')) || any(contains(DLCoptions,'1.2')) || any(contains(DLCoptions,'all'))
    [output,output1p2]=runIECDLC1p1_1p2(params,matData,output,DLCoptions); % run DLC 1.1 50-year occurrence load
    save output1p2 output1p2
end
if any(contains(DLCoptions,'1.3')) || any(contains(DLCoptions,'all'))
    [output,output1p3]=runIECDLC1p3(params,output); % run DLC 1.3 UTS, ETM
    save output1p3 output1p3
end
if any(contains(DLCoptions,'1.4')) || any(contains(DLCoptions,'all'))
    [output,output1p4]=runIECDLC1p4(params,output); % run DLC 1.4 ECD with direction change
    save output1p4 output1p4
end
if any(contains(DLCoptions,'1.5')) || any(contains(DLCoptions,'all'))
    [output,output1p5]=runIECDLC1p5(params,output); % run DLC 1.5 EWS
    save output1p5 output1p5
end
if any(contains(DLCoptions,'6.1')) || any(contains(DLCoptions,'all'))
    [output,output6p1]=runIECDLC6p1(params,output); % run DLC 6.1 EWM 50-year parked
    save output6p1 output6p1
end
if any(contains(DLCoptions,'6.2')) || any(contains(DLCoptions,'all'))
    [output,output6p2]=runIECDLC6p2(params,output); % run DLC 6.2 EWM 50-year parked, no yaw control
    save output6p2 output6p2
end
if any(contains(DLCoptions,'6.3')) || any(contains(DLCoptions,'all'))
    [output,output6p3]=runIECDLC6p3(params,output); % run DLC 6.3 EWM 1-year parked
    save output6p3 output6p3
end

if contains(DLCoptions,'loads') 
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

% Save the structure output file
save output output

% Write results to Excel file
xlsfn='IECDLC_Results.xlsx';

%ble:<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% fix errors from optimization run by not saving data
try writeResult(output,xlsfn);
catch, error('error with saving') % do nothing
end
%ble:>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% Remove Simulink controller path -- removed at end of script so it can change
rmpath(genpath(params.simulinkModelFolder));
end

%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
function writeResult(output,xlsfn)
    names=fieldnames(output);
    numDLC = size(output,1);
    xlsData=[];

    csvfn = strrep(xlsfn,'.xlsx','.csv');
    fid  = fopen(csvfn,'w');
    for i=1:length(names)
    % %     xlsData=[xlsData; {'DLC Name',names{i},'Occurring at time (sec)','on Channel','in File'}];
        combinedOut.(names{i}) = [output.(names{i})];
        hdrTxt = {'DLC Name',names{i},'Occurring at time (sec)','on Channel','in File'};
        fprintf(fid, '%s,%s,%s,%s,%s\n', hdrTxt{1},hdrTxt{2},hdrTxt{3},hdrTxt{4},hdrTxt{5});
        for j=1:numDLC
    % %         xlsData=[xlsData;...
    % %             {getfield(combinedOut,names{i},{j},'Name'),...
    % %             getfield(combinedOut,names{i},{j},'data'),...
    % %             getfield(combinedOut,names{i},{j},'time'),...
    % %             getfield(combinedOut,names{i},{j},'Chan'),...
    % %             getfield(combinedOut,names{i},{j},'File')}];
            dlcTxt = {getfield(combinedOut,names{i},{j},'Name'),...
                getfield(combinedOut,names{i},{j},'data'),...
                getfield(combinedOut,names{i},{j},'time'),...
                getfield(combinedOut,names{i},{j},'Chan'),...
                getfield(combinedOut,names{i},{j},'File')};
            fprintf(fid, '%s,%f,%f,%s,%s\n', dlcTxt{1},dlcTxt{2},dlcTxt{3},dlcTxt{4},dlcTxt{5});
        end
    end
    fclose(fid);

    % xlswrite(xlsfn,xlsData,1);

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
