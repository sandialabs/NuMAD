function varargout = readWTPerf(input_file)
%READWTPERFINPUT  Read a WT_Perf input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   wtpinput = ReadWTPerfInput(input_file_name)
% 
%   Works with v3.04.00c of WT Perf
%   Returns a structure containing the data in a WT_Perf input file.

% open file
fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"\n',input_file);
end

line=fgetl(fid);
wtpinput.title{1}=fgetl(fid);
wtpinput.title{2}=fgetl(fid);
line=fgetl(fid);
wtpinput.Echo=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.DimenInp=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.Metric=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
line=fgetl(fid);
wtpinput.NumSect=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
wtpinput.MaxIter=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
wtpinput.ATol=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.SWTol=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
line=fgetl(fid);
wtpinput.TipLoss=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.HubLoss=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.Swirl=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.SkewWake=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.AdvBrake=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.IndProp=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.AIDrag=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.TIDrag=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
line=fgetl(fid);
wtpinput.NumBlade=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
wtpinput.RotorRad=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.HubRad=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.PreCone=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.Tilt=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.Yaw=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.HubHt=fscanf(fid,'%g',[1 1]); line=fgetl(fid);

wtpinput.NumSeg=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
line=fgetl(fid);
for n=1:wtpinput.NumSeg
    wtpinput.seg.RElm(n)=fscanf(fid,'%g',[1 1]);
    wtpinput.seg.Twist(n)=fscanf(fid,'%g',[1 1]);
    wtpinput.seg.Chord(n)=fscanf(fid,'%g',[1 1]);
    wtpinput.seg.AFfile(n)=fscanf(fid,'%d',[1 1]);
    wtpinput.seg.PrntElem{n}=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
end

line=fgetl(fid);
wtpinput.Rho=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.KinVisc=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.ShearExp=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
wtpinput.UseCm=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.UseCpmin=fscanf(fid,'%s',[1 1]); line=fgetl(fid);

wtpinput.NumAF=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
for n=1:wtpinput.NumAF
    wtpinput.AF_File{n}=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
end

line=fgetl(fid);

wtpinput.UnfPower=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.TabDel=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.OutNines=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.Beep=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.KFact=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.WriteBED=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.InputTSR=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.OutMaxCp=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
wtpinput.SpdUnits=fscanf(fid,'%s',[1 1]); line=fgetl(fid);

line=fgetl(fid);
wtpinput.NumCases=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
line=fgetl(fid);
for n=1:wtpinput.NumCases
    wtpinput.case{n}.WSorTSR=fscanf(fid,'%g',[1 1]);
    wtpinput.case{n}.RotSpd=fscanf(fid,'%g',[1 1]);
    wtpinput.case{n}.Pitch=fscanf(fid,'%g',[1 1]); line=fgetl(fid);
end

if (wtpinput.NumCases == 0)
    line=fgetl(fid);
    wtpinput.par.ParRow=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
    wtpinput.par.ParCol=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
    wtpinput.par.ParTab=fscanf(fid,'%d',[1 1]); line=fgetl(fid);
    wtpinput.par.OutPwr=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
    wtpinput.par.OutCp=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
    wtpinput.par.OutTrq=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
    wtpinput.par.OutFlp=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
    wtpinput.par.OutThr=fscanf(fid,'%s',[1 1]); line=fgetl(fid);
    wtpinput.par.Pit=fscanf(fid,'%g,',[1 3]); line=fgetl(fid);
    wtpinput.par.Omg=fscanf(fid,'%g,',[1 3]); line=fgetl(fid);
    wtpinput.par.Spd=fscanf(fid,'%g,',[1 3]); line=fgetl(fid);
end

% close file
status=fclose(fid);

if (nargout==0)
    subplot(3,1,1);
    plot(wtpinput.seg.RElm,wtpinput.seg.Twist);
    title('Twist');
    subplot(3,1,2);
    plot(wtpinput.seg.RElm,wtpinput.seg.Chord);
    title('Chord');
    subplot(3,1,3);
    plot(wtpinput.seg.RElm,wtpinput.seg.AFfile,'s');
    title('Airfoil File Number');
else
    varargout{1} = wtpinput;
end