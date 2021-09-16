function output = simBasic(in)
% SIMBASIC(IN) used in conjuction with other scripts to automate the process
% of running simple simulations involving FAST and AeroDyn input files, as 
% well as Crunch analysis.
%**********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
% Examine this script for definition of 'in' data structure.
%  

fst=readFastMain(sprintf('%s.fst',in.fst_name));
ad=readFastAD(sprintf('%s.ipt',in.ad_name));

fst.SimCtrl.TMax=in.simTime;
fst.Out.SttsTime=round(in.simTime/10);
ad.WindFile=sprintf('"%s.wnd"',in.wind_name);
writeFastAD(ad,sprintf('%s.ipt',in.ad_name));
writeFastMain(fst,sprintf('%s.fst',in.fst_name));
dos(sprintf('%s %s.fst',in.fst_path,in.fst_name));

fn{1}=sprintf('%s.out',in.fst_name);
prepCrunch('crunch.cru',fn{1},in.RFids,fn,in.fst_name,in.EEGrps);
dos([in.crunch_path ' cru.cru']);

end


