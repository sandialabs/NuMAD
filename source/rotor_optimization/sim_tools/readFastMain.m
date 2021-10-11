function fast = readFastMain(input_file)
%READFASTMAIN  Read a FAST primary input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   fast = ReadFastInput('input_file_name') 
%   Returns a structure containing the data in a FAST input file.

if ~exist('input_file','var')
    [fn pn] = uigetfile( ...
        {'*.fst','FAST input (*.fst)'; ...
        '*.*','All files (*.*)'});
    if isequal(fn,0)
        return
    else
        input_file = [pn, fn];
    end
end

%% open file
fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"\n',input_file);
    return
end
fast.ffn = input_file;

%% get file description from header
fgetl(fid);
fgetl(fid);
fast.title{1,1}=fgetl(fid);
fast.title{2,1}=fgetl(fid);

%% simulation control
fgetl(fid);
fast.SimCtrl.Echo         =rip(fid,'TF');
fast.SimCtrl.ADAMSPrep    =rip(fid,'str',{'1','2','3'});
fast.SimCtrl.AnalMode     =rip(fid,'str',{'1','2'});
fast.SimCtrl.NumBl        =rip(fid,'num');
fast.SimCtrl.TMax         =rip(fid,'num');
fast.SimCtrl.DT           =rip(fid,'num');

%% turbine control
fgetl(fid);
fast.TurbCtrl.YCMode        =rip(fid,'str',{'0','1','2'});
fast.TurbCtrl.TYCOn         =rip(fid,'num');
fast.TurbCtrl.PCMode        =rip(fid,'str',{'0','1','2'});
fast.TurbCtrl.TPCOn         =rip(fid,'num');
fast.TurbCtrl.VSContrl      =rip(fid,'str',{'0','1','2','3'});
fast.TurbCtrl.VS_RtGnSp     =rip(fid,'num');
fast.TurbCtrl.VS_RtTq       =rip(fid,'num');
fast.TurbCtrl.VS_Rgn2K      =rip(fid,'num');
fast.TurbCtrl.VS_SlPc       =rip(fid,'num');
fast.TurbCtrl.GenModel      =rip(fid,'str',{'1','2','3'});
fast.TurbCtrl.GenTiStr      =rip(fid,'TF');
fast.TurbCtrl.GenTiStp      =rip(fid,'TF');
fast.TurbCtrl.SpdGenOn      =rip(fid,'num');
fast.TurbCtrl.TimGenOn      =rip(fid,'num');
fast.TurbCtrl.TimGenOf      =rip(fid,'num');
fast.TurbCtrl.HSSBrMode     =rip(fid,'str',{'1','2'});
fast.TurbCtrl.THSSBrDp      =rip(fid,'num');
fast.TurbCtrl.TiDynBrk      =rip(fid,'num');
fast.TurbCtrl.TTpBrDp(1,1)  =rip(fid,'num');
fast.TurbCtrl.TTpBrDp(2,1)  =rip(fid,'num');
fast.TurbCtrl.TTpBrDp(3,1)  =rip(fid,'num');
fast.TurbCtrl.TBDepISp(1,1) =rip(fid,'num');
fast.TurbCtrl.TBDepISp(2,1) =rip(fid,'num');
fast.TurbCtrl.TBDepISp(3,1) =rip(fid,'num');
fast.TurbCtrl.TYawManS      =rip(fid,'num');
fast.TurbCtrl.TYawManE      =rip(fid,'num');
fast.TurbCtrl.NacYawF       =rip(fid,'num');
fast.TurbCtrl.TPitManS(1,1) =rip(fid,'num');
fast.TurbCtrl.TPitManS(2,1) =rip(fid,'num');
fast.TurbCtrl.TPitManS(3,1) =rip(fid,'num');
fast.TurbCtrl.TPitManE(1,1) =rip(fid,'num');
fast.TurbCtrl.TPitManE(2,1) =rip(fid,'num');
fast.TurbCtrl.TPitManE(3,1) =rip(fid,'num');
fast.TurbCtrl.BlPitch(1,1)  =rip(fid,'num');
fast.TurbCtrl.BlPitch(2,1)  =rip(fid,'num');
fast.TurbCtrl.BlPitch(3,1)  =rip(fid,'num');
fast.TurbCtrl.BlPitchF(1,1) =rip(fid,'num');
fast.TurbCtrl.BlPitchF(2,1) =rip(fid,'num');
fast.TurbCtrl.BlPitchF(3,1) =rip(fid,'num');

%% environmental conditions
fgetl(fid);
fast.Env.Gravity          =rip(fid,'num');

%% feature flags
fgetl(fid);
fast.Flags.FlapDOF1       =rip(fid,'TF');
fast.Flags.FlapDOF2       =rip(fid,'TF');
fast.Flags.EdgeDOF        =rip(fid,'TF');
fast.Flags.TeetDOF        =rip(fid,'TF');
fast.Flags.DrTrDOF        =rip(fid,'TF');
fast.Flags.GenDOF         =rip(fid,'TF');
fast.Flags.YawDOF         =rip(fid,'TF');
fast.Flags.TwFADOF1       =rip(fid,'TF');
fast.Flags.TwFADOF2       =rip(fid,'TF');
fast.Flags.TwSSDOF1       =rip(fid,'TF');
fast.Flags.TwSSDOF2       =rip(fid,'TF');
fast.Flags.CompAero       =rip(fid,'TF');
fast.Flags.CompNoise      =rip(fid,'TF');

%% initial conditions
fgetl(fid);
fast.Init.OoPDefl         =rip(fid,'num');
fast.Init.IPDefl          =rip(fid,'num');
fast.Init.TeetDefl        =rip(fid,'num');
fast.Init.Azimuth         =rip(fid,'num');
fast.Init.RotSpeed        =rip(fid,'num');
fast.Init.NacYaw          =rip(fid,'num');
fast.Init.TTDspFA         =rip(fid,'num');
fast.Init.TTDspSS         =rip(fid,'num');

%% turbine configuration
fgetl(fid);
fast.TurbConf.TipRad        =rip(fid,'num');
fast.TurbConf.HubRad        =rip(fid,'num');
fast.TurbConf.PSpnElN       =rip(fid,'num');
fast.TurbConf.UndSling      =rip(fid,'num');
fast.TurbConf.HubCM         =rip(fid,'num');
fast.TurbConf.OverHang      =rip(fid,'num');
fast.TurbConf.NacCMxn       =rip(fid,'num');
fast.TurbConf.NacCMyn       =rip(fid,'num');
fast.TurbConf.NacCMzn       =rip(fid,'num');
fast.TurbConf.TowerHt       =rip(fid,'num');
fast.TurbConf.Twr2Shft      =rip(fid,'num');
fast.TurbConf.TwrRBHt       =rip(fid,'num');
fast.TurbConf.ShftTilt      =rip(fid,'num');
fast.TurbConf.Delta3        =rip(fid,'num');
fast.TurbConf.PreCone(1,1)  =rip(fid,'num');
fast.TurbConf.PreCone(2,1)  =rip(fid,'num');
fast.TurbConf.PreCone(3,1)  =rip(fid,'num');
fast.TurbConf.AzimB1Up      =rip(fid,'num');

%% mass and inertia
fgetl(fid);
fast.MassProp.YawBrMass   =rip(fid,'num');
fast.MassProp.NacMass     =rip(fid,'num');
fast.MassProp.HubMass     =rip(fid,'num');
fast.MassProp.TipMass(1,1)=rip(fid,'num');
fast.MassProp.TipMass(2,1)=rip(fid,'num');
fast.MassProp.TipMass(3,1)=rip(fid,'num');
fast.MassProp.NacYIner    =rip(fid,'num');
fast.MassProp.GenIner     =rip(fid,'num');
fast.MassProp.HubIner     =rip(fid,'num');

%% drivetrain
fgetl(fid);
fast.DrvTrn.GBoxEff        =rip(fid,'num');
fast.DrvTrn.GenEff         =rip(fid,'num');
fast.DrvTrn.GBRatio        =rip(fid,'num');
fast.DrvTrn.GBRevers       =rip(fid,'TF');
fast.DrvTrn.HSSBrTqF       =rip(fid,'num');
fast.DrvTrn.HSSBrDT        =rip(fid,'num');
fast.DrvTrn.DynBrkFi       =rip(fid,'qstr');
fast.DrvTrn.DTTorSpr       =rip(fid,'num');
fast.DrvTrn.DTTorDmp       =rip(fid,'num');

%% simple induction generator
fgetl(fid);
fast.SIG.SIG_SlPc         =rip(fid,'num');
fast.SIG.SIG_SySp         =rip(fid,'num');
fast.SIG.SIG_RtTq         =rip(fid,'num');
fast.SIG.SIG_PORt         =rip(fid,'num');

%% thevenin-equivalent induction generator
fgetl(fid);
fast.TEC.TEC_Freq         =rip(fid,'num');
fast.TEC.TEC_NPol         =rip(fid,'num');
fast.TEC.TEC_SRes         =rip(fid,'num');
fast.TEC.TEC_RRes         =rip(fid,'num');
fast.TEC.TEC_VLL          =rip(fid,'num');
fast.TEC.TEC_SLR          =rip(fid,'num');
fast.TEC.TEC_RLR          =rip(fid,'num');
fast.TEC.TEC_MR           =rip(fid,'num');

%% platform model
fgetl(fid);
fast.Ptfm.PtfmModel       =rip(fid,'str',{'0','1','2','3'});
fast.Ptfm.PtfmFile        =rip(fid,'qstr');

%% tower
fgetl(fid);
fast.Twr.TwrNodes         =rip(fid,'num');
fast.Twr.TwrFile          =rip(fid,'qstr');

%% nacelle-yaw
fgetl(fid);
fast.Yaw.YawSpr           =rip(fid,'num');
fast.Yaw.YawDamp          =rip(fid,'num');
fast.Yaw.YawNeut          =rip(fid,'num');

%% furling
fgetl(fid);
fast.Furl.Furling         =rip(fid,'TF');
fast.Furl.FurlFile        =rip(fid,'qstr');

%% rotor-teeter
fgetl(fid);
fast.Teet.TeetMod         =rip(fid,'num');
fast.Teet.TeetDmpP        =rip(fid,'num');
fast.Teet.TeetDmp         =rip(fid,'num');
fast.Teet.TeetCDmp        =rip(fid,'num');
fast.Teet.TeetSStP        =rip(fid,'num');
fast.Teet.TeetHStP        =rip(fid,'num');
fast.Teet.TeetSSSp        =rip(fid,'num');
fast.Teet.TeetHSSp        =rip(fid,'num');

%% tip-brake
fgetl(fid);
fast.TpBr.TBDrConN        =rip(fid,'num');
fast.TpBr.TBDrConD        =rip(fid,'num');
fast.TpBr.TpBrDT          =rip(fid,'num');

%% blade
fgetl(fid);
fast.BldFile{1,1}         =rip(fid,'qstr');
fast.BldFile{2,1}         =rip(fid,'qstr');
fast.BldFile{3,1}         =rip(fid,'qstr');

%% aerodyn
fgetl(fid);
fast.ADFile               =rip(fid,'qstr');

%% noise
fgetl(fid);
fast.NoiseFile            =rip(fid,'qstr');

%% adams
fgetl(fid);
fast.ADAMSFile            =rip(fid,'qstr');

%% linearization file
fgetl(fid);
fast.LinFile              =rip(fid,'qstr');

%% output
fgetl(fid);
fast.Out.SumPrint         =rip(fid,'TF');
try
filePosition = ftell(fid);
fast.Out.OutFileFmt       =rip(fid,'num');
catch ME
    % If reading a 'num' failed, we must have a v7.00 file that does not 
    % have OutFileFmt. Need to rewind file before reading next parameter.
    if isequal(ME.identifier,'readFastMain:ripCheck')
        fseek(fid,filePosition,'bof'); % rewind
    else
        rethrow(ME); % rethrow error if other error occured
    end
end
fast.Out.TabDelim         =rip(fid,'TF');
fast.Out.OutFmt           =rip(fid,'qstr');
fast.Out.TStart           =rip(fid,'num');
fast.Out.DecFact          =rip(fid,'num');
fast.Out.SttsTime         =rip(fid,'num');
fast.Out.NcIMUxn          =rip(fid,'num');
fast.Out.NcIMUyn          =rip(fid,'num');
fast.Out.NcIMUzn          =rip(fid,'num');
fast.Out.ShftGagL         =rip(fid,'num');
fast.Out.NTwGages         =rip(fid,'num');
fast.Out.TwrGagNd         =rip(fid,'csv');
fast.Out.NBlGages         =rip(fid,'num');
fast.Out.BldGagNd         =rip(fid,'csv');
fgetl(fid);
fast.OutList              =ripOutList(fid);

fclose(fid);
end

%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function param = rip(fid,type,varargin)
% read input file parameter
    line = strtrim(fgetl(fid));       % get the next line of text
    options = cell(0,1);              % initialize to empty cell
    if nargin > 2
        options=varargin{1};          % deal out the input arguments
    end
    switch type  % handling of input depends on type
        case 'TF'  % type: true/false
            % use sscanf to get the first non-whitespace 
            % group of characters
            param = sscanf(line,'%s',1);
            % examine only the first character to determine true/false
            if upper(param(1))=='T'
                param = 'True';
            elseif upper(param(1))=='F'
                param = 'False';
            else
                error('readFastMain:ripCheck','Unknown option "%s" for input line:\n%s\n',param,line);
            end
        case 'str'   % type: unquoted string
            % use sscanf to get the first non-whitespace
            % group of characters
            param = sscanf(line,'%s',1);
            % try to find a match from the options list
            k=strmatch(upper(param),upper(options));
            if ~isempty(k)
                param=options{k};
            else
                error('readFastMain:ripCheck','Unknown option "%s" for input line:\n%s\n',param,line);
            end
        case 'qstr'  % type: quoted string
            % a quoted string could contain spaces, 
            % which rules out using sscanf\
            % regexp match: beginning of string ^, quote ", 
            %               anything not quote [^"]*, quote "
            s=regexp(line,'^"[^"]*"','match','once');
            if isempty(s)
                param='';          % a quoted string was not found
            else
                param=s;  % trim off any leading whitespace
            end
        case 'num'  % type: single numeric value OR string option (i.e. 'default')
            param = sscanf(line,'%f',1);
            if isempty(param)  
                % numeric value not found; try to read a string parameter
                if isempty(regexp(line,'^\s*"','once'))  % check for quotes
                    % it's an unquoted string
                    param = sscanf(line,'%s',1);
                else
                    % it's a quoted string
                    s=regexp(line,'^"[^"]*"','match','once');
                    param=regexp(s,'[^"]*','match','once'); % remove quotes         
                end
                % try to find a match from the options list
                k=strmatch(upper(param),upper(options));
                if ~isempty(k)
                    param=options{k};
                else
                    error('readFastMain:ripCheck','Unknown option "%s" for input line:\n%s\n',param,line);
                end
            end
        case 'csv'  % type: comma separated value list
            % get the input field
            % regexp match: beginning of string ^, set of zero or more []* 
            %               digits \d whitespace \s or ,.-+
            csv = regexp(line,'^[\d\s,.-+]*','match','once');
            csv = regexp(csv,'[^\s,]+','match');
            param = zeros(1,length(csv));
            for k=1:length(csv)
                param(k) = str2double(csv{k});
            end
    end
end

function outs = ripOutList(fid)
    outs = textscan(fid,'%s','Delimiter','\n');
    outs = outs{1};
    for k=1:length(outs)
        if strncmp(outs{k},'END',3)
            kk=k;
        end
    end
    outs(kk:end) = [];
    
%     outs = transpose(fread(fid,inf,'uint8=>char'));
%     k = regexp(outs,'END','once');
%     outs = outs(1:k-1);
end