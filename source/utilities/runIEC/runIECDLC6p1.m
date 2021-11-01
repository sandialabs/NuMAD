function [output,outputAllFiles]=runIECDLC6p1(params,output)
global fastPath
global adamsPath
CaseName='IECDLC6p1EWM50';

% Switches:
% params.fastsim='fast simulink';

% MAIN PROCESS
if false % run full set
    windLabels={'EWM50+00.wnd',...
        'EWM50+05.wnd',...
        'EWM50+10.wnd',...
        'EWM50+15.wnd',...
        'EWM50-05.wnd',...
        'EWM50-10.wnd',...
        'EWM50-15.wnd'};
else % subset, for efficiency
    windLabels={'EWM50+15.wnd' 'EWM50-15.wnd'};
    warning('runIECDLC6p1 is running only a single case in this configuration.  Please verify that it is the correct case.')
end

parfor n=1:length(windLabels)
    
    if params.simulate
        % set up a format so that files can be opened by parallel runs.
        t = getCurrentTask;
        if isempty(t)
            id = ''; % first iteration of genetic algorithm in generation
        else
            id = ['_parallelWorker_' num2str(t.ID)]; % remaining population are sent to workers in parallel
            % copy the FAST and AD input files for parallel runs
% %             disp([params.fstfn '.fst'])
% %             copyfile([params.fstfn '.fst'], [params.fstfn id '.fst'])
% %             fst=readFastMain(['IEC_' params.fstfn '.fst']);
% %             copyfile(strrep(fst.ADFile,'"',''), [fst.ADFile(2:end-8) id '_AD.ipt'])
        end
        
        fst=readFastMain(['IEC_' params.fstfn '.fst']);%fst=readFastMain(['IEC_' params.fstfn id '.fst']);
        ad=readFastAD([fst.ADFile(2:end-8) '_AD.ipt']);%ad=readFastAD([fst.ADFile(2:end-8) id '_AD.ipt']);
        % prepare and perform aeroelastic simulation
        ad.WindFile=['"wind\' windLabels{n} '"'];
        fst.SimCtrl.TMax=300;
        fst.Out.TStart=250;
        
        % DLC 6.1-specific settings:
        ptch=90;
        fst.TurbCtrl.BlPitch(1)=ptch;
        fst.TurbCtrl.BlPitch(2)=ptch;
        fst.TurbCtrl.BlPitch(3)=ptch;
        fst.TurbCtrl.BlPitchF(1)=ptch;
        fst.TurbCtrl.BlPitchF(2)=ptch;
        fst.TurbCtrl.BlPitchF(3)=ptch;
        fst.TurbCtrl.TPitManS(1)=0;
        fst.TurbCtrl.TPitManS(2)=0;
        fst.TurbCtrl.TPitManS(3)=0;
        fst.TurbCtrl.TPitManE(1)=0;
        fst.TurbCtrl.TPitManE(2)=0;
        fst.TurbCtrl.TPitManE(3)=0;
        fst.Flags.GenDOF='False';
        fst.Init.RotSpeed=0;
        ad.IndModel='NONE';
        %         ad.StallMod='STEADY';
        % end DLC 6.1-specific settings
        
        fst.ADFile=['"' CaseName '_' params.fstfn id '_AD.ipt' '"'];% ble: change the AD filename
        thisFastName=[CaseName '_' params.fstfn id];
        if strcmp(params.fastsim, 'fast simulink')
            fst.TurbCtrl.YCMode   = 0; % do not yaw for this Design Load Case
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
                outname=[CaseName '_' windLabels{n}(1:end-4) '.out'];
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
                runFASTusingSimulink(fastFilenameSFunc, params.simulinkModel)
                outname=[CaseName '_' windLabels{n}(1:end-4) '.out'];
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
                writeFastMain(fst,[params.fstfn '_4ADAMS.fst']);
                dos([fastPath ' ' params.fstfn '_4ADAMS.fst']);
                dos([adamsPath ' ' params.fstfn '_4ADAMS_ADAMS.acf']);
                disp('Starting ADAMS model from Matlab.....')
                outname=[CaseName '_' windLabels{n}(1:end-4) '.out'];
                try movefile([thisFastName '_4ADAMS_ADAMS.plt'], ['out/' outname])
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
for  n=1:length(windLabels)
    outName{end+1}=['out\' CaseName '_' windLabels{n}(1:end-4) '.out'];
end

% initialize the output2 variable
output2 = cell(1,length(outName));
for cc = 1:length(outName)
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
ct=1; outputAllFiles = cell(1,length(windLabels));
for ii=1
    for jj=1:length(windLabels)
        outputAllFiles{ii,jj} = output2{ct};
        ct = ct+1;
    end
end

end