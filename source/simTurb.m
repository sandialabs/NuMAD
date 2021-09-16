function output = simTurb(in)
% SIMTURB(IN) used in conjuction with other scripts to automate the process
% of running turbulent 3d wind field for fatigue damage calculations.
%**********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
% Examine this script for definition of 'in' data structure.
%  

% determine random seeds for turbulent wind files
try
    load ranset
catch
    disp('No "ranset" variable found...generating random seeds for use');
    ranset=rand(in.numSeeds,1)*1e6;
    save ranset ranset
end

for w=in.ws  % loop on all wind speeds
    for s=1:in.numSeeds  % loop on all wind seeds
        tsim=readTurbSim(sprintf('%s.inp',in.wind_name));
        fst=readFastMain(sprintf('%s.fst',in.fst_name));
        ad=readFastAD(sprintf('%s.ipt',in.ad_name));
        
        % TurbSim stuff
        tsim.HubHt=fst.TurbConf.TowerHt+fst.TurbConf.Twr2Shft+fst.TurbConf.OverHang*-1*sind(-1*fst.TurbConf.ShftTilt);
        tsim.GridWidth=1.1*2*fst.TurbConf.TipRad;
        tsim.GridHeight=1.1*2*fst.TurbConf.TipRad;  % (maybe round it up)
        tsim.RefHeight=tsim.HubHt;
        tsim.title='tsimtmp';
        tsim.URef=w;
        tsim.RandSeed1=ranset(s);
        writeTurbSim(tsim,sprintf('%s.inp',in.wind_name));
        dos(sprintf('%s %s.inp',in.tsim_path,in.wind_name));
        
        fst.SimCtrl.TMax=in.simTime;
        fst.Out.SttsTime=round(in.simTime/10);
        ad.WindFile=sprintf('"%s.wnd"',in.wind_name);
        writeFastAD(ad,sprintf('%s.ipt',in.ad_name));
        writeFastMain(fst,sprintf('%s.fst',in.fst_name));
        dos(sprintf('%s %s.fst',in.fst_path,in.fst_name));
        
        fn{s}=sprintf('%i_%s_%i.out',w,in.fst_name,s);
        copyfile(sprintf('%s.out',in.fst_name),fn{s});
    end
    prepCrunch('crunch.cru',sprintf('%s.out',in.fst_name),in.RFids,fn,sprintf('%i_%s',w,in.fst_name),in.EEGrps);
    dos([in.crunch_path ' cru.cru']);
end

end
