function output = prepCrunch(crunchfile,fastoutfile,RFcols,fns,AggRoot,EEGrps)
%PREPCRUNCH  Set up the Crunch data structure and write the Crunch input file
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% output = prepCrunch(crunchfile,fastoutfile,RFcols,fns,AggRoot,EEGrps)
%   
%   See simTurb() for usage example.
% 
%   crunchfile = name of Crunch input file to be written
%   fastoutfile = file name of an example Fast output file; used to gather
%           information about the data content that will be processed
%   RFcols = array of column id's corresponding to data in the Fast .out
%           file
%   

data=loadFastOutData(fastoutfile);
crunch=readCrunch(crunchfile);

crunch.ChanInfo='';
crunch.ChanInfo.NumInCols=length(data.list);
if RFcols(1)~=1  % make sure that the first column from the out file (Time) is used in this analysis
    RFcols=[1 RFcols];
end
crunch.ChanInfo.NumCols=length(RFcols);
for i=1:crunch.ChanInfo.NumCols
    crunch.ChanInfo.ChanTitle{i}=['"' data.list{RFcols(i)} '"'];
    crunch.ChanInfo.ChanUnits{i}=['"' data.units{RFcols(i)} '"'];
    crunch.ChanInfo.OrigChan(i)=RFcols(i);
    crunch.ChanInfo.Scale(i)=1;
    crunch.ChanInfo.Offset(i)=0;
end
crunch.JobOpt.AggRoot=AggRoot;
crunch.ColNums=(2:length(RFcols));
crunch.NumRFCols=crunch.ChanInfo.NumCols-1;
crunch.HalfCycMult=ones(crunch.NumRFCols,1)*0.5;
crunch.EEGrps=EEGrps;
crunch.NumFiles=length(fns);  
crunch.InFiles=fns;
writeCrunch(crunch,'cru.cru');
end