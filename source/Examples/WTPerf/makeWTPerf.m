delete SNLV27_WTPerf.oup

fast=readFastMain('SNLV27.fst');
ad=readFastAD('SNLV27_AD.ipt');

wt = FastAD2WTPerf(fast,ad);

wt.NumSect=10;    
wt.MaxIter=50;
wt.PreCone= 0;
wt.Tilt= 0;
wt.Yaw= 0;
wt.par.OutFlp='False';
wt.par.OutPwr='False';
wt.par.OutThr='False';
wt.par.OutTrq='False';
wt.par.Pit=[-2 8 0.5];
rotspd = 30;  %
wt.par.Omg=[rotspd rotspd 0];
wt.par.Spd=[1 11 0.5];

fn='SNLV27_WTPerf.wtp';
writeWTPerf(wt,fn);
wtp=['C:\DesignCodes\WT_Perf_v3.10\WT_Perf.exe ' fn];
dos(wtp);
Cp_surfs


