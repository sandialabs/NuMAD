function wtp = readWTPerfInput(filename,wtperfVersion)
%READWTPERFINPUT  Read a WT_Perf input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   wtp = ReadWTPerfInput('filename',[wtperfVersion])
%      Returns a structure containing the data in a WT_Perf input file.
%      Works with WT_Perf v3.10 and v3.05.00a
%      It will attempt to determine file format version by searching
%      the file header lines.
%
%      filename = filename of WTPerf file to be created.
%      wtperfVersion = (optional) select between 'v310' and 
%                      'v3.05.00' formats


if nargin < 1
    wtp = wtperfStruct('v3.05.00');
    return
elseif isempty(filename)
    wtp = wtperfStruct(wtperfVersion);
    return
end

fr = fileReader(filename);  % create a file reader object
fr.discardLine; % -----  WT_Perf Input File  -----
title1 = fr.getLine;
title2 = fr.getLine;

% determine file version from header lines
wtperfCompatible = regexp([title1,title2],'v[.\d]+','match');
WTPcomp = '';
if ~isempty(wtperfCompatible)
    switch wtperfCompatible{1}
        case {'v3.10','v3.00','v310','v300'}
            WTPcomp = 'v310';
        case {'v3.05','v3.05.00'}
            WTPcomp = 'v3.05.00';
    end
end

% use the found version if wtperfVersion not specified
% default to 'v3.05.00' if version was not found.
if ~exist('wtperfVersion','var') || isempty(wtperfVersion)
    if isempty(WTPcomp)
        wtperfVersion = 'v3.05.00';  % flag to switch the output file format
    else
        wtperfVersion = WTPcomp;
    end
end

% set WTPver variable which controls program flow
switch wtperfVersion
    case {'v310','v300'}
        WTPver = 'v310';
    case {'v3.05','v3.05.00'}
        WTPver = 'v3.05.00';
    otherwise
        error('wtperfVersion "%s" not recognized',wtperfVersion);
end

wtp = wtperfStruct(WTPver); % initialize structure to preserve field order

wtp.title{1} = title1;
wtp.title{2} = title2;
               fr.discardLine; % -----  Input Configuration  -----
wtp.Echo     = fr.str;
wtp.DimenInp = fr.str;
wtp.Metric   = fr.str;
               fr.discardLine; % -----  Model Configuration  -----
wtp.NumSect  = fr.num;
wtp.MaxIter  = fr.num;
switch WTPver
    case 'v310'
wtp.NSplit   = 35; % trying to have an appropriate default
    case 'v3.05.00'
wtp.NSplit   = fr.num;
end
wtp.ATol     = fr.num;
wtp.SWTol    = fr.num;
               fr.discardLine; % -----  Algorithm Configuration  -----
wtp.TipLoss  = fr.str;
wtp.HubLoss  = fr.str;
wtp.Swirl    = fr.str;
wtp.SkewWake = fr.str;
switch WTPver
    case 'v310'
wtp.AdvBrake = fr.str;
wtp.IndProp  = fr.str;
wtp.IndType  = wtp.IndProp; % trying to have an appropriate default
    case 'v3.05.00'
IndType  = fr.str;
wtp.AdvBrake = 'False'; % trying to have an appropriate default
wtp.IndProp  = IndType; % trying to have an appropriate default
wtp.IndType  = IndType;
end
wtp.AIDrag   = fr.str;
wtp.TIDrag   = fr.str;
switch WTPver
    case 'v310'
wtp.TISingularity = 'True';
wtp.DAWT          = 'False';
wtp.Cavitation    = 'False';
wtp.PressAtm      = 101325;
wtp.PressVapor    = 2500;
wtp.CavSF         = 1.0;
wtp.WatDepth      = 33.0;            
    case 'v3.05.00'
wtp.TISingularity = fr.str;
wtp.DAWT     = fr.str;
wtp.Cavitation = fr.str;
               fr.discardLine; % -----  Cavitation Model  -----
wtp.PressAtm = fr.num;
wtp.PressVapor = fr.num;
wtp.CavSF    = fr.num;
wtp.WatDepth = fr.num;
end
               fr.discardLine; % -----  Turbine Data  -----
wtp.NumBlade = fr.num;
wtp.RotorRad = fr.num;
wtp.HubRad   = fr.num;
wtp.PreCone  = fr.num;
wtp.Tilt     = fr.num;
wtp.Yaw      = fr.num;
wtp.HubHt    = fr.num;

wtp.NumSeg   = fr.num;
               fr.discardLine; % RElm   Twist     Chord   AFfile  PrntElem
for n=1:wtp.NumSeg
    wtp.seg.RElm(n)=fscanf(fr.fileID,'%g',[1 1]);
    wtp.seg.Twist(n)=fscanf(fr.fileID,'%g',[1 1]);
    wtp.seg.Chord(n)=fscanf(fr.fileID,'%g',[1 1]);
    wtp.seg.AFfile(n)=fscanf(fr.fileID,'%d',[1 1]);
    wtp.seg.PrntElem{n}=fscanf(fr.fileID,'%s',[1 1]); line=fgetl(fr.fileID);
end

               fr.discardLine; % -----  Aerodynamic Data  -----
wtp.Rho      = fr.num;
wtp.KinVisc  = fr.num;
wtp.ShearExp = fr.num;
wtp.UseCm    = fr.str;
switch WTPver
    case 'v310'
wtp.UseCpmin = 'False'; % trying to have an appropriate default
    case 'v3.05.00'
wtp.UseCpmin = fr.str;
end

wtp.NumAF    = fr.num;
for n=1:wtp.NumAF
    wtp.AF_File{n} = fr.str;
end
               fr.discardLine; % -----  Output Configuration  -----
switch WTPver
    case 'v310'
wtp.UnfPower = 'False'; % trying to have an appropriate default
    case 'v3.05.00'
wtp.UnfPower = fr.str;
end
wtp.TabDel   = fr.str;
switch WTPver
    case 'v310'
wtp.ConvFlag = 1;       % trying to have an appropriate default
wtp.Beep     = 'False'; % trying to have an appropriate default
    case 'v3.05.00'
wtp.ConvFlag = fr.num;
wtp.Beep     = fr.str;
end
wtp.KFact    = fr.str;
wtp.WriteBED = fr.str;
wtp.InputTSR = fr.str;
switch WTPver
    case 'v310'
wtp.OutMaxCp = 'True';
    case 'v3.05.00'
wtp.OutMaxCp = fr.str;
end
wtp.SpdUnits = fr.str;

               fr.discardLine; % -----  Combined-Case Analysis  -----
wtp.NumCases = fr.num;
               fr.discardLine; % WS or TSR   RotSpd   Pitch
if wtp.NumCases>0
   for n=1:wtp.NumCases
       wtp.case{n}.WSorTSR=fscanf(fr.fileID,'%g',[1 1]);
       wtp.case{n}.RotSpd=fscanf(fr.fileID,'%g',[1 1]);
       wtp.case{n}.Pitch=fscanf(fr.fileID,'%g',[1 1]); line=fgetl(fr.fileID);
   end
end

               fr.discardLine; % -----  Parametric Analysis (Ignored if NumCases > 0 )  -----
if (wtp.NumCases == 0)
    wtp.par.ParRow = fr.num;
    wtp.par.ParCol = fr.num;
    wtp.par.ParTab = fr.num;
    wtp.par.OutPwr = fr.str;
    wtp.par.OutCp  = fr.str;
    wtp.par.OutTrq = fr.str;
    wtp.par.OutFlp = fr.str;
    wtp.par.OutThr = fr.str;
    wtp.par.Pit=fscanf(fr.fileID,'%g,',[1 3]); line=fgetl(fr.fileID);
    wtp.par.Omg=fscanf(fr.fileID,'%g,',[1 3]); line=fgetl(fr.fileID);
    wtp.par.Spd=fscanf(fr.fileID,'%g,',[1 3]); line=fgetl(fr.fileID);
end

% close file
fr.delete;

end


function wtp = wtperfStruct(WTPver)
% -----  Input Configuration  -----
wtp.originalVersion = WTPver;
wtp.title{1} = 'WT_Perf input file';
wtp.title{2} = sprintf('Compatible with WT_Perf %s',WTPver);
wtp.Echo     = 'False';
wtp.DimenInp = 'False';
wtp.Metric   = 'True';
% -----  Model Configuration  -----
wtp.NumSect  = 1;
switch WTPver
    case 'v310'
wtp.MaxIter  = 20000;
    case 'v3.05.00'
wtp.MaxIter  = 13;
end
wtp.NSplit   = 35;
wtp.ATol     = 1e-5;
wtp.SWTol    = 1e-5;
% -----  Algorithm Configuration  -----
wtp.TipLoss  = 'True';
wtp.HubLoss  = 'True';
wtp.Swirl    = 'True';
wtp.SkewWake = 'True';
wtp.AdvBrake = 'False';
wtp.IndProp  = 'True';
wtp.IndType  = 'True';
wtp.AIDrag   = 'True';
wtp.TIDrag   = 'True';
wtp.TISingularity = 'True';
wtp.DAWT          = 'False';
wtp.Cavitation    = 'False';
% -----  Cavitation Model  -----
wtp.PressAtm      = 101325;
wtp.PressVapor    = 2500;
wtp.CavSF         = 1.0;
wtp.WatDepth      = 33.0;            
% -----  Turbine Data  -----
wtp.NumBlade = 3;
wtp.RotorRad = 0;
wtp.HubRad   = 0;
wtp.PreCone  = 0;
wtp.Tilt     = 0;
wtp.Yaw      = 0;
wtp.HubHt    = 2;
wtp.NumSeg   = 0;
wtp.seg.RElm = [];
wtp.seg.Twist = [];
wtp.seg.Chord = [];
wtp.seg.AFfile = [];
wtp.seg.PrntElem = cell(0);


% -----  Aerodynamic Data  -----
wtp.Rho      = 1.225;
wtp.KinVisc  = 1.4639e-5;
wtp.ShearExp = 0.0;
wtp.UseCm    = 'False';
wtp.UseCpmin = 'False';
wtp.NumAF    = 0;
wtp.AF_File = cell(0);

% -----  Output Configuration  -----
wtp.UnfPower = 'False';
wtp.TabDel   = 'True';
wtp.ConvFlag = 1;
wtp.Beep     = 'False';
wtp.KFact    = 'True';
wtp.WriteBED = 'True';
wtp.InputTSR = 'True';
wtp.OutMaxCp = 'True';
wtp.SpdUnits = '"mps"';

% -----  Combined-Case Analysis  -----
wtp.NumCases = 0;
wtp.case{1} = struct('WSorTSR',[],'RotSpd',[],'Pitch',[]);

% -----  Parametric Analysis (Ignored if NumCases > 0 )  -----
wtp.par.ParRow = 3;
wtp.par.ParCol = 2;
wtp.par.ParTab = 1;
wtp.par.OutPwr = 'False';
wtp.par.OutCp  = 'True';
wtp.par.OutTrq = 'False';
wtp.par.OutFlp = 'False';
wtp.par.OutThr = 'False';
wtp.par.Pit=[0, 0, 0];
wtp.par.Omg=[20, 20, 0];
wtp.par.Spd=[8, 8, 0];

end
