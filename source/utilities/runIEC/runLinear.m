function output=runLinear(params)
global fastPath

CaseName='Linear';

if params.simulate
    fprintf('IECSweep_%s.out \n',params.fstfn)
    SFunc = '_SFunc';
    %ble: change so that it can read SFunc outputs.
    out=loadFASTOutData(sprintf(['IECSweep_%s' SFunc '.out'],params.fstfn));
    for ws=[0 params.lin]
        fst=readFastMain([params.fstfn '.fst']);
        ad=readFastAD(strrep(fst.ADFile,'"',''));
        
        ad.WindFile=sprintf('"wind\\steady_%d.wnd"',ws);
        fst.SimCtrl.TMax=1000;
        fst.SimCtrl.AnalMode=2;
        fst.TurbCtrl.YCMode=0;
        fst.TurbCtrl.PCMode=0;
        fst.TurbCtrl.GenTiStr='True';
        fst.TurbCtrl.GenTiStp='True';
        fst.TurbCtrl.TimGenOn=0;
        fst.TurbCtrl.TimGenOf=9999.9;
        fst.TurbCtrl.THSSBrDp=9999.9;
        fst.TurbCtrl.TiDynBrk=9999.9;
        fst.TurbCtrl.TTpBrDp(1)=9999.9;
        fst.TurbCtrl.TTpBrDp(2)=9999.9;
        fst.TurbCtrl.TTpBrDp(3)=9999.9;
        fst.TurbCtrl.TYawManS=9999.9;
        fst.TurbCtrl.TPitManS(1)=9999.9;
        fst.TurbCtrl.TPitManS(2)=9999.9;
        fst.TurbCtrl.TPitManS(1)=9999.9;
        fst.TurbCtrl.TBDepISp(1)=9999.9;
        fst.TurbCtrl.TBDepISp(2)=9999.9;
        fst.TurbCtrl.TBDepISp(3)=9999.9;
        fst.Flags.CompNoise='False';
        ad.StallMod='STEADY';
        ad.InfModel='EQUIL';
        
        if ws==0 % if doing the static turbine linearization
            fst.TurbConfig.GenDOF='False';
            spd=0;
            fst.Init.RotSpeed=spd;
            fst.Flags.CompAero='False';
            ad.IndModel='NONE';
        else
            fst.TurbConfig.GenDOF='True';
            id1=strmatch('WindVxi',out.list);
            id2=strmatch('RotSpeed',out.list);
            pointer=min(find(out.data(:,id1)>=ws));
            spd=out.data(pointer,id2);
            fst.Init.RotSpeed=spd;
            fst.Flags.CompAero='True';
            ad.IndModel='SWIRL';
        end
        
        disp(sprintf('Linearization: target rotor speed at %f4.1 m/s is %5.2f rpm',ws,spd));
        fst.ADFile=['"' CaseName '_' strrep(fst.ADFile,'"','') '"'];
        thisFastName=[CaseName '_' params.fstfn]
        writeFastMain(fst,[thisFastName '.fst']);
        writeFastAD(ad,strrep(fst.ADFile,'"',''));
        
        disp('Starting FAST model from Matlab.....')
        dos([fastPath ' ' thisFastName '.fst']);
        outname=sprintf('%s_%dmps.lin',CaseName,ws);
        copyfile([thisFastName '.lin'],['out/' outname]);
        disp('LINearization output files saved under new names.')
        disp(' ')
        disp('          *********** Done **************')
        disp(' ')
    end
end

end