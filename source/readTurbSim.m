function tsim = readTurbSim(input_file)
%READTURBSIM  Read a TurbSim input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   tsim = readTurbSim(input_file_name)
%     Mean for use with TurbSim v1.50
%

if ~exist('input_file','var')
    [fn pn] = uigetfile( ...
        {'*.inp','TurbSim input (*.inp)'; ...
        '*.*','All files (*.*)'});
    if isequal(fn,0)
        return
    else
        input_file = [pn, fn];
    end
end

% open file
fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"\n',input_file);
    return
end
tsim.ffn = input_file;

% get file description from header
tsim.title=fgetl(fid);
% try to determine version of TurbSim input
if strfind(tsim.title,'v1.5')
    tsim.ver = 'v1.5';
elseif strfind(tsim.title,'v1.3')
    tsim.ver = 'v1.3';
else
    tsim.ver = 'unknown';
end
fgetl(fid);

% runtime options
fgetl(fid);  
tsim.RandSeed1       =rip(fid,'num');
tsim.RandSeed2       =rip(fid,'num',{'RanLux','RNSNLW'});
tsim.WrBHHTP         =rip(fid,'TF');
tsim.WrFHHTP         =rip(fid,'TF');
tsim.WrADHH          =rip(fid,'TF');
tsim.WrADFF          =rip(fid,'TF');
tsim.WrBLFF          =rip(fid,'TF');
tsim.WrADTWR         =rip(fid,'TF');
tsim.WrFMTFF         =rip(fid,'TF');
tsim.WrACT           =rip(fid,'TF');
tsim.Clockwise       =rip(fid,'TF');
switch tsim.ver
    case 'v1.3'
    tsim.ScaleIEC    =rip(fid,'TF');
    otherwise
    tsim.ScaleIEC    =rip(fid,'str',{'0','1','2'});
end
fgetl(fid);

% turbine/model specifications
fgetl(fid);
tsim.NumGrid_Z       =rip(fid,'num');
tsim.NumGrid_Y       =rip(fid,'num');
tsim.TimeStep        =rip(fid,'num');
tsim.AnalysisTime    =rip(fid,'num');
tsim.UsableTime      =rip(fid,'num');
tsim.HubHt           =rip(fid,'num');
tsim.GridHeight      =rip(fid,'num');
tsim.GridWidth       =rip(fid,'num');
tsim.VFlowAng        =rip(fid,'num');
tsim.HFlowAng        =rip(fid,'num');
fgetl(fid);

% meteorological boundary conditions
fgetl(fid);
tsim.TurbModel       =rip(fid,'str',{'IECKAI','IECVKM','GP_LLJ','NWTCUP','SMOOTH','WF_UPW','WF_07D','WF_14D','NONE'});
tsim.IECstandard     =rip(fid,'str',{'1','1-ED2','1-ED3','2','3'});
tsim.IECturbc        =rip(fid,'num',{'A','B','C','KHTEST'});
tsim.IEC_WindType    =rip(fid,'str',{'NTM','1ETM','2ETM','3ETM','1EWM1','2EWM1','3EWM1','1EWM50','2EWM50','3EWM50'});
tsim.ETMc            =rip(fid,'num',{'default'});
tsim.WindProfileType =rip(fid,'str',{'JET','LOG','PL','default'});
tsim.RefHt           =rip(fid,'num');
tsim.URef            =rip(fid,'num');
tsim.ZJetMax         =rip(fid,'num',{'default'});
tsim.PLExp           =rip(fid,'num',{'default'});
tsim.Z0              =rip(fid,'num',{'default'});
fgetl(fid);

% non-IEC meteorological boundary conditions
fgetl(fid);
tsim.Latitude        =rip(fid,'num',{'default'});
tsim.RICH_NO         =rip(fid,'num');
tsim.UStar           =rip(fid,'num',{'default'});
tsim.ZI              =rip(fid,'num',{'default'});
tsim.PC_UW           =rip(fid,'num',{'default','none'});
tsim.PC_UV           =rip(fid,'num',{'default','none'});
tsim.PC_VW           =rip(fid,'num',{'default','none'});
tsim.IncDec1         =rip(fid,'num',{'default'});
tsim.IncDec2         =rip(fid,'num',{'default'});
tsim.IncDec3         =rip(fid,'num',{'default'});
tsim.CohExp          =rip(fid,'num',{'default'});
fgetl(fid);

% coherent turbulence scaling parameters
fgetl(fid);
tsim.CTEventPath     =rip(fid,'qstr');
tsim.CTEventFile     =rip(fid,'str',{'random','les','dns'});
tsim.Randomize       =rip(fid,'TF');
tsim.DistScl         =rip(fid,'num');
tsim.CTLy            =rip(fid,'num');
tsim.CTLz            =rip(fid,'num');
tsim.CTStartTime     =rip(fid,'num');
fgetl(fid);

fclose(fid);
end

function param = rip(fid,type,varargin)
% read input file parameter
    %param = fscanf(fid,'%s',[1 1]); descrip=fgetl(fid);
    %param = regexp(param,'[^"]*','match','once'); %remove quotes
    n=regexp(fgetl(fid),'\s*(?<param>"[^"]*"|\S*)\s+(?<descrip>.*)','names');
    n.param = regexp(n.param,'[^"]*','match','once'); %remove quotes
    options = cell(0,1);
    if nargin > 2
        options=varargin{1};
    end
    switch type
        case 'TF'
            if upper(n.param(1))=='T'
                param = 'True';
            elseif upper(n.param(1))=='F'
                param = 'False';
            end
        case 'str'
            k=strmatch(strtrim(upper(n.param)),upper(options));
            if ~isempty(k)
                param=options{k};
            else
                error('Unknown option "%s" for input line:\n%s\n',n.param,strtrim(n.descrip));
            end
        case 'qstr'
            param=sprintf('"%s"',n.param);
        case 'num'
            param = cell2mat(textscan(n.param,'%f'));
            if isempty(param)
                k=strmatch(strtrim(upper(n.param)),upper(options));
                if ~isempty(k)
                    param=options{k};
                else
                    error('Unknown option "%s" for input line:\n%s\n',n.param,strtrim(n.descrip));
                end
            end
    end
end