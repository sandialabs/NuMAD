function cru = readCrunch(input_file)
%READCRUNCH  Read a Crunch input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   cru = ReadCrunch('input_file_name')
%   Returns a structure containing the data in a Crunch input file.
%    Currently does not work with the following Crunch features:
%      -Calculated Channels
%      -Moving Averages
%      -Load Roses
%      -Azimuth Averages
%      -Crosstalk removal
%      -Peak Finding
%      -Peak and valley listing
%      -Probability Mass
%      -Extreme events
%      -Summary files
%      -Statistical extrapolation

if ~exist('input_file','var')
    [fn pn] = uigetfile( ...
        {'*.cru','Crunch input (*.cru)'; ...
        '*.*','All files (*.*)'});
    if isequal(fn,0)
        return
    else
        input_file = [pn, fn];
    end
end

fid=fopen(input_file);
if (fid == -1)
    error('Could not open input "%s"\n',input_file);
    return
end
cru.ffn = input_file;


% begin reading input file
fgetl(fid);
cru.title{1,1}=fgetl(fid);

% Job Options
fgetl(fid);
cru.JobOpt.Echo     =rip(fid,'TF');
cru.JobOpt.Out_Stats=rip(fid,'TF');
cru.JobOpt.Out_Data =rip(fid,'TF');
cru.JobOpt.TabDelim =rip(fid,'TF');
cru.JobOpt.RealFmt  =rip(fid,'qstr');
cru.JobOpt.Aggregate=rip(fid,'TF');
cru.JobOpt.AggRoot  =rip(fid,'qstr');

% Input-Data Layout
fgetl(fid);
cru.IptDataLay.CTRow=rip(fid,'num');
cru.IptDataLay.CURow=rip(fid,'num');
cru.IptDataLay.FDRow=rip(fid,'num');
cru.IptDataLay.NumRecs=rip(fid,'num');
cru.IptDataLay.TStartTEnd=rip(fid,'csv');

% Channel information
fgetl(fid);
cru.ChanInfo.NumInCols=rip(fid,'num');
cru.ChanInfo.NumCols=rip(fid,'num');
fgetl(fid);
Table           =textscan(fid,'%s %s %f %f %f',cru.ChanInfo.NumCols);
cru.ChanInfo.ChanTitle=Table{1};
cru.ChanInfo.ChanUnits=Table{2};
cru.ChanInfo.OrigChan=Table{3};
cru.ChanInfo.Scale=Table{4};
cru.ChanInfo.Offset=Table{5};
fgetl(fid);

% Filtering
fgetl(fid);
cru.Filter.NumFilt=rip(fid,'num');
cru.Filter.FiltCols=rip(fid,'csv');
cru.Filter.FiltType=rip(fid,'num');
cru.Filter.LoCut=rip(fid,'num');
cru.Filter.HiCut=rip(fid,'num');

% Calculated Channels (currently not supported)
fgetl(fid);
cru.CalcChan.NumCChan=rip(fid,'num');
cru.CalcChan.Seed=rip(fid,'num');
fgetl(fid);
Table = textscan(fid,'%s %s %s',cru.CalcChan.NumCChan);
cru.CalcChan.ColTitle = Table{1};
cru.CalcChan.Units    = Table{2};
cru.CalcChan.Equation = Table{3};
if cru.CalcChan.NumCChan > 0
    fgetl(fid);
end
% Moving Averages (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Time and Wind Speed
fgetl(fid);
cru.TimeCol=rip(fid,'num');
cru.WS_Col =rip(fid,'num');

% Load Roses (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Azimuth Averages (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Crosstalk (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Peak finding (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Peak and valley listing (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Probability mass (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Rainflow cycles
fgetl(fid);
cru.NumRFCols            =rip(fid,'num');

fgetl(fid);  % features not used; just using Crunch to do raw cycle counting
fgetl(fid);  % features not used; just using Crunch to do raw cycle counting
fgetl(fid);  % features not used; just using Crunch to do raw cycle counting
fgetl(fid);  % features not used; just using Crunch to do raw cycle counting
fgetl(fid);  % features not used; just using Crunch to do raw cycle counting

fgetl(fid);
Table           =textscan(fid,'%f %f %f %f %f',cru.NumRFCols);
cru.HalfCycMult=Table{2};
fgetl(fid);


% Extreme Events (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Summary Files (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% Statistical Extrapolation (currently not supported)
fgetl(fid);
fgetl(fid);
fgetl(fid);

% input files
fgetl(fid);
cru.NumFiles                    =rip(fid,'num');
for j=1:cru.NumFiles
    cru.InFiles{j}                  =rip(fid,'qstr');
end

fclose(fid);
end

%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function param = rip(fid,type,varargin)
% read input file parameter
try
    line = strtrim(fgetl(fid));       % get the next line of text
catch
    line = '';
    disp('Warning: found blank line while looking for parameter');
end
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
            error('Unknown option "%s" for input line:\n%s\n',param,line);
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
            error('Unknown option "%s" for input line:\n%s\n',param,line);
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
                error('Unknown option "%s" for input line:\n%s\n',param,line);
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
