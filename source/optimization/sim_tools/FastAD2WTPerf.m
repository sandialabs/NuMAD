function wt = FastAD2WTPerf(fast,ad)
%FASTAD2WTPERF  Create WTPerf data structurefrom information in the Fast 
%       and AeroDyn inputs  
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   wt = FastAD2WTPerf(fast,ad)
%   
%   fast = FAST data structure as created by readFastMain()
%   ad = AeroDyn data structure as created by readFastAD()
% 

wt.title= {'WT_Perf input file title',  'Line two of title'};
wt.Echo= 'False';
wt.DimenInp= 'True';
wt.Metric= 'True';
wt.NumSect= 1;
wt.MaxIter= 30;
wt.ATol= 1.0000e-006;
wt.SWTol= 1.0000e-006;
wt.NSplit=30;
wt.TipLoss= 'True';
wt.HubLoss= 'True';
wt.Swirl= 'True';
wt.SkewWake= 'True';
wt.AdvBrake= 'True';
wt.IndType= 'True';
wt.AIDrag= 'True';
wt.TIDrag= 'True';
wt.TISingularity= 'True';
wt.DAWT= 'False';
wt.Cavitation= 'False';
wt.PressAtm=101325;
wt.PressVapor=2500;
wt.CavSF=1.0;
wt.WatDepth=33.0; 
wt.UnfPower='False';
wt.ConvFlag=1;
wt.Beep='False';

wt.NumBlade= fast.SimCtrl.NumBl;
wt.RotorRad= fast.TurbConf.TipRad;
wt.HubRad= fast.TurbConf.HubRad;
wt.PreCone= fast.TurbConf.PreCone(1);
wt.Tilt= fast.TurbConf.ShftTilt;
wt.Yaw= 0;
wt.HubHt= fast.TurbConf.TowerHt + fast.TurbConf.Twr2Shft+fast.TurbConf.OverHang*sind(fast.TurbConf.ShftTilt);
wt.NumSeg= ad.BldNodes;
wt.seg.RElm=ad.RNodes;
wt.seg.Twist=ad.AeroTwst;
wt.seg.Chord=ad.Chord;
wt.seg.AFfile=ad.NFoil;
for i=1:wt.NumSeg
    wt.seg.PrntElem{i}='True';
end
wt.Rho= ad.Rho;
wt.KinVisc= ad.KinVisc;
wt.ShearExp= 0;

if strcmpi(ad.UseCm,'no_cm')
    wt.UseCm= 'False';
else
    wt.UseCm= 'True';
end
wt.UseCpmin='False';

wt.NumAF= ad.NumFoil;
for i=1:wt.NumAF
%     wt.AF_File{i}=strrep(ad.FoilNm{i},'AD_','WTP_');
    wt.AF_File{i}=ad.FoilNm{i};
end
wt.TabDel= 'True';
wt.KFact= 'True';
wt.WriteBED= 'True';
wt.InputTSR= 'True';
wt.SpdUnits= '"mps"';
wt.NumCases= 0;
wt.OutMaxCp='False';

% ParRow:                    Row parameter    (1-rpm, 2-pitch, 3-tsr/speed).
% ParCol:                    Column parameter (1-rpm, 2-pitch, 3-tsr/speed).
% ParTab:                    Table parameter  (1-rpm, 2-pitch, 3-tsr/speed).
% OutPwr:                    Request output of rotor power?
% OutCp:                     Request output of Cp?
% OutTrq:                    Request output of shaft torque?
% OutFlp:                    Request output of flap bending moment?
% OutThr:                    Request output of rotor thrust?
% PitSt, PitEnd, PitDel:     First, last, delta blade pitch (deg).
% OmgSt, OmgEnd, OmgDel:     First, last, delta rotor speed (rpm).
% SpdSt, SpdEnd, SpdDel:     First, last, delta speeds.

wt.par.ParRow=3;
wt.par.ParCol=2;
wt.par.ParTab=1;
wt.par.OutPwr= 'True';
wt.par.OutCp= 'True';
wt.par.OutTrq= 'True';
wt.par.OutFlp= 'True';
wt.par.OutThr= 'True';
wt.par.Pit= [2 8 1];
wt.par.Omg= [10 14 .5];
wt.par.Spd= [3 12.5 0.50];

end
