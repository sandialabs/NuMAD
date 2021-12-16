function twr = readFastTower(input_file)
%READFASTTOWER  Read a FAST tower input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   twr = readFastTower('input_file_name')
%   Returns a structure containing the data in a FAST tower input file.
%

if ~exist('input_file','var')
    [fn pn] = uigetfile( ...
        {'*.dat','Blade/Tower data (*.dat)'; ...
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
twr.ffn = input_file;


% begin reading input file
fgetl(fid);
fgetl(fid);
twr.title{1,1}=fgetl(fid);
fgetl(fid);  % TOWER PARAMETERS
twr.NTwInpSt        =rip(fid,'num');
twr.CalcTMode       =rip(fid,'TF');
twr.TwrFADmp(1)     =rip(fid,'num');
twr.TwrFADmp(2)     =rip(fid,'num');
twr.TwrSSDmp(1)     =rip(fid,'num');
twr.TwrSSDmp(2)     =rip(fid,'num');
fgetl(fid);  % TOWER ADJUSTMENT FACTORS
twr.FAStTunr(1)     =rip(fid,'num');
twr.FAStTunr(2)     =rip(fid,'num');
twr.SSStTunr(1)     =rip(fid,'num');
twr.SSStTunr(2)     =rip(fid,'num');
twr.AdjTwMa         =rip(fid,'num');
twr.AdjFASt         =rip(fid,'num');
twr.AdjSSSt         =rip(fid,'num');
fgetl(fid);  % DISTRIBUTED TOWER PROPERTIES
fgetl(fid);  % column labels
fgetl(fid);  % column units
PropTable           =riptbl(fid,repmat('%f ',1,10));
twr.prop.HtFract    = PropTable{ 1};
twr.prop.TMassDen   = PropTable{ 2};
twr.prop.TwFAStif   = PropTable{ 3};
twr.prop.TwSSStif   = PropTable{ 4};
twr.prop.TwGJStif   = PropTable{ 5};
twr.prop.TwEAStif   = PropTable{ 6};
twr.prop.TwFAIner   = PropTable{ 7};
twr.prop.TwSSIner   = PropTable{ 8};
twr.prop.TwFAcgOf   = PropTable{ 9};
twr.prop.TwSScgOf   = PropTable{10};
fgetl(fid);  % TOWER FORE-AFT MODE SHAPES
twr.TwFAM1Sh(2)     =rip(fid,'num');
twr.TwFAM1Sh(3)     =rip(fid,'num');
twr.TwFAM1Sh(4)     =rip(fid,'num');
twr.TwFAM1Sh(5)     =rip(fid,'num');
twr.TwFAM1Sh(6)     =rip(fid,'num');
twr.TwFAM2Sh(2)     =rip(fid,'num');
twr.TwFAM2Sh(3)     =rip(fid,'num');
twr.TwFAM2Sh(4)     =rip(fid,'num');
twr.TwFAM2Sh(5)     =rip(fid,'num');
twr.TwFAM2Sh(6)     =rip(fid,'num');
fgetl(fid);  % TOWER SIDE-TO-SIDE MODE SHAPES
twr.TwSSM1Sh(2)     =rip(fid,'num');
twr.TwSSM1Sh(3)     =rip(fid,'num');
twr.TwSSM1Sh(4)     =rip(fid,'num');
twr.TwSSM1Sh(5)     =rip(fid,'num');
twr.TwSSM1Sh(6)     =rip(fid,'num');
twr.TwSSM2Sh(2)     =rip(fid,'num');
twr.TwSSM2Sh(3)     =rip(fid,'num');
twr.TwSSM2Sh(4)     =rip(fid,'num');
twr.TwSSM2Sh(5)     =rip(fid,'num');
twr.TwSSM2Sh(6)     =rip(fid,'num');

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

function table = riptbl(fid,frmt)
table = textscan(fid,frmt);
end