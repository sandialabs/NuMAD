function [output,outputAllFiles]=runIECDLC1p5(params,output)
global fastPath
global adamsPath  %not used?

CaseName='IECDLC1p5EWS';

% Switches:
% params.fastsim='fast simulink';

% create the wind file labels
ws=(params.operatingPoints(1):1:params.operatingPoints(3));
windLabels={};
for i=1:length(ws)
    windLabels{end+1}=sprintf('EWSV+%3.1f.wnd',ws(i));
    windLabels{end+1}=sprintf('EWSV-%3.1f.wnd',ws(i));
end

% ====================== MAIN PROCESS - SIMULATIONS =======================
% create a counter for parallel runs
for ii = 1:length(params.yaw)
    ctr.windLabels((ii-1)*length(windLabels)+1:ii*length(windLabels)) = windLabels;
    ctr.y((ii-1)*length(windLabels)+1:ii*length(windLabels)) = params.yaw(ii);
end
% loop through the "wind direction" as well
ctr.y1 = ctr.y; 
dirMat = repmat(params.wd, numel(ctr.windLabels), 1);
ctr.d = reshape(dirMat, 1, numel(dirMat));
ctr.windLabels = repmat(ctr.windLabels, 1, length(params.wd));
ctr.y = repmat(ctr.y, 1, length(params.wd));


parfor ii=1:length(ctr.windLabels)
    % ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % set up a format so that files can be opened by parallel runs.
    t = getCurrentTask;
    if isempty(t)
        id = ''; % first iteration of genetic algorithm in generation
    else
        id = ['_parallelWorker_' num2str(t.ID)]; % remaining population are sent to workers in parallel
        if params.simulate
            % copy the FAST and AD input files for parallel runs
% %             disp([params.fstfn '.fst'])
% %             copyfile([params.fstfn '.fst'], [params.fstfn id '.fst'])
% %             fst=readFastMain([params.fstfn id '.fst']);
% %             copyfile(strrep(fst.ADFile,'"',''), [fst.ADFile(2:end-8) id '_AD.ipt'])
        end
    end
    % ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
    if params.simulate
        fst=readFastMain(['IEC_' params.fstfn '.fst']);%fst=readFastMain([params.fstfn id '.fst']);
        ad=readFastAD([fst.ADFile(2:end-8) '_AD.ipt']);%ad=readFastAD([fst.ADFile(2:end-8) id '_AD.ipt']);
        % prepare and perform aeroelastic simulation
        ad.WindFile=['"' params.parDir 'wind\' ctr.windLabels{ii} '"'];
        fst.SimCtrl.TMax=params.delay+100;
        fst.ADFile=['"' CaseName '_' params.fstfn id '_AD.ipt' '"'];% ble: change the AD filename
        if ~strcmp(params.fastsim,'fast simulink')
            fst=runIEC_setIC(fst,params.operatingPoints(2),params.parDir);
        end
        fst.Init.NacYaw = ctr.y(ii);
        thisFastName=[CaseName '_' params.fstfn id];
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
            fst.TurbCtrl.YCMode   = 0; % no yaw control so that the yaw misalignment will be tested
            fst.TurbCtrl.PCMode   = 1;
            fst.TurbCtrl.VSContrl = 1;            
        end
        writeFastMain(fst,[thisFastName '.fst']);
        writeFastAD(ad,strrep(fst.ADFile,'"',''));
        
        switch params.fastsim
            case 'fast'
                disp('Starting FAST model from Matlab.....')
                dos([fastPath ' ' thisFastName '.fst']);
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' ctr.windLabels{ii}(1:end-4) '_dir' num2str(ctr.d(ii)) '.out'];
                try movefile([thisFastName '.out'],[params.parDir 'out/' outname]);
                catch err; warning(err.message); end
                try movefile([thisFastName '.elm'],[params.parDir 'out/' outname(1:end-4) '.elm']);
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
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' ctr.windLabels{ii}(1:end-4) '_dir' num2str(ctr.d(ii)) '.out'];
                try movefile([thisFastName '_SFunc.out'],[params.parDir 'out/' outname]);
                catch err; warning(err.message); end
                try movefile([thisFastName '_SFunc.elm'],[params.parDir 'out/' outname(1:end-4) '.elm']);
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
                outname=[CaseName '_yaw' num2str(ctr.y(ii)) '_' ctr.windLabels{ii}(1:end-4) '.out'];
                try movefile([thisFastName '_4ADAMS_ADAMS.plt'], [params.parDir 'out/' outname])
                catch err; warning(err.message); end
                disp('AD and ADAMS output files saved under new names.')
                disp(' ')
                disp('          *********** Done **************')
                disp(' ')

        end
    end    
end

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

% Find extreme values for tower clearance and root bending moments
outName={};
for yaw = params.yaw
    for n=1:length(windLabels)
        for d=params.wd
            outName{end+1}=[params.parDir 'out\' CaseName '_yaw' num2str(yaw) '_' ...
                windLabels{n}(1:end-4) '_dir' num2str(ctr.d(ii)) '.out'];
        end
    end
end

% initialize the output2 variable
output2 = cell(1,length(outName));
for cc = 1:length(outName)
%     output2{cc} = output;
    output2{cc}.MaxBladeFlapBendingMoment=[]; % set up output structure
    output2{cc}.MaxFlapStrain=[];
end
% determine the maximum/minimum loads for each output file
parfor cc = 1:length(outName)
    output2{cc} = sortExtremes(params,output2{cc},CaseName,outName(cc),cs,EIs,EA_normal);
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
    output_new.(names{ii}) = output2{index}.(names{ii});
end
if size(fieldnames(output),1) == size(fieldnames(output_new),1)
    output = [output; output_new];
else
    output = output_new;
end

% save output2 results in a 2x2 matrix that is more intuitive to search
ct=1; outputAllFiles = cell(length(params.yaw),length(windLabels));
for ii = 1:length(params.yaw)
    for jj=1:length(windLabels)
        outputAllFiles{ii,jj} = output2{ct};
        ct = ct+1;
    end
end

end