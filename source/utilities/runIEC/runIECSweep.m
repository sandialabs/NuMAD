function output=runIECSweep(params,output)
global fastPath
global adamsPath
CaseName='IECSweep';

% Switches:
% params.fastsim='fast simulink';

% MAIN PROCESS
windLabels={'ramp.wnd'};

for n=1:length(windLabels)
    
    if params.simulate
        fst=readFastMain(['IEC_' params.fstfn '.fst']);
        ad=readFastAD(strrep(fst.ADFile,'"',''));
        % prepare and perform aeroelastic simulation
        ad.WindFile=['"' params.parDir 'wind\' windLabels{n} '"'];
        fst.SimCtrl.TMax=1000;
        fst.ADFile=['"' CaseName '_' params.fstfn '_AD.ipt' '"'];% ble: change the AD filename
        thisFastName=[CaseName '_' params.fstfn];
        if strcmp(params.fastsim, 'fast simulink')
            fst.TurbCtrl.YCMode = 2;
            fst.TurbCtrl.TYCOn = 0;
            fst.TurbCtrl.PCMode = 2;
            fst.TurbCtrl.VSContrl = 3;
        else
            fst.TurbCtrl.YCMode = 0;
            fst.TurbCtrl.PCMode = 1;
            fst.TurbCtrl.VSContrl = 1;            
        end
        writeFastMain(fst,[thisFastName '.fst']);
        writeFastAD(ad,strrep(fst.ADFile,'"',''));

        switch params.fastsim
            case 'fast'
                disp('Starting FAST model from Matlab.....')
                dos([fastPath ' ' thisFastName '.fst']);
                outname=[CaseName '_' windLabels{n}(1:end-4) '.out'];
                movefile([thisFastName '.out'],[params.parDir 'out/' outname]);
                delete([thisFastName '.elm'], [thisFastName '.fsm'], [thisFastName '.opt']);
                delete([thisFastName '.fst'], [thisFastName '_AD.ipt']);
                disp('AD and FAST output files saved under new names.')
                disp(' ')
                disp('          *********** Done **************')
                disp(' ')
            case 'fast simulink'
                disp('Starting FAST SFunc model from Matlab.....')
                fastFilenameSFunc = [thisFastName '.fst'];
                YawPosBias = params.wd(1);
                runFASTusingSimulink(fastFilenameSFunc, params.simulinkModel, YawPosBias)
                outname=[CaseName '_' windLabels{n}(1:end-4) '.out'];
                movefile([thisFastName '_SFunc.out'],[params.parDir 'out/' outname]);
                delete([thisFastName '_SFunc.elm'], [thisFastName '_SFunc.fsm'], [thisFastName '_SFunc.opt']);
                delete([thisFastName '.fst'], [thisFastName '_AD.ipt']);
                disp('AD and FAST_SFunc output files saved under new names.')
                disp(' ')
                disp('          *********** Done **************')
                disp(' ')
            case 'adams'
                fast.SimCtrl.ADAMSPrep='2';
                fast.Init.OoPDefl=0;
                fast.Init.IPDefl=0;
                fast.Init.TeetDefl=0;
                fast.Init.TTDspFA=0;
                fast.Init.TTDspSS=0;
                writeFastMain(fast,[params.fstfn '_4ADAMS.fst']);
                dos([fastPath ' ' params.fstfn '_4ADAMS.fst']);
                dos([adamsPath ' ' params.fstfn '_4ADAMS_ADAMS.acf']);
                disp('Starting ADAMS model from Matlab.....')
                outname=[CaseName '_' windLabels '.out'];
                eval(['movefile ' params.fstfn '_4ADAMS_ADAMS.plt ' params.parDir 'out/' outname]);
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
save spans spans

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
    outName{end+1}=[params.parDir 'out/' CaseName '_' windLabels{n}(1:end-4) '.out'];
end
output = sortExtremes(params,output,CaseName,outName,cs,EIs,EA_normal);

end