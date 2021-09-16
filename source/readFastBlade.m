function blade = readFastBlade(input_file)
%READFASTBLADE  Read a FAST blade input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   blade = readFastBlade('input_file_name') 
%   Returns a structure containing the data in a FAST blade input file.

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
blade.ffn = input_file;


% begin reading input file
fgetl(fid);
fgetl(fid);
blade.title{1,1}=fgetl(fid);
fgetl(fid);  % BLADE PARAMETERS
blade.NBlInpSt      =rip(fid,'num');
blade.CalcBMode     =rip(fid,'TF');
blade.BldFlDmp(1)   =rip(fid,'num');
blade.BldFlDmp(2)   =rip(fid,'num');
blade.BldEdDmp(1)   =rip(fid,'num');
fgetl(fid);  % BLADE ADJUSTMENT FACTORS
blade.FlStTunr(1)   =rip(fid,'num');
blade.FlStTunr(2)   =rip(fid,'num');
blade.AdjBlMs       =rip(fid,'num');
blade.AdjFlSt       =rip(fid,'num');
blade.AdjEdSt       =rip(fid,'num');
fgetl(fid);  % DISTRIBUTED BLADE PROPERTIES
fgetl(fid);  % column labels
fgetl(fid);  % column units
PropTable           =riptbl(fid,repmat('%f ',1,17));
blade.prop.BlFract   = PropTable{ 1};
blade.prop.AeroCent  = PropTable{ 2};
blade.prop.StrcTwst  = PropTable{ 3};
blade.prop.BMassDen  = PropTable{ 4};
blade.prop.FlpStff   = PropTable{ 5};
blade.prop.EdgStff   = PropTable{ 6};
blade.prop.GJStff    = PropTable{ 7};
blade.prop.EAStff    = PropTable{ 8};
blade.prop.Alpha     = PropTable{ 9};
blade.prop.FlpIner   = PropTable{10};
blade.prop.EdgIner   = PropTable{11};
blade.prop.PrecrvRef = PropTable{12};
blade.prop.PreswpRef = PropTable{13};
blade.prop.FlpcgOf   = PropTable{14};
blade.prop.EdgcgOf   = PropTable{15};
blade.prop.FlpEAOf   = PropTable{16};
blade.prop.EdgEAOf   = PropTable{17};
fgetl(fid);  % BLADE MODE SHAPES
blade.BldFl1Sh(2)   =rip(fid,'num');
blade.BldFl1Sh(3)   =rip(fid,'num');
blade.BldFl1Sh(4)   =rip(fid,'num');
blade.BldFl1Sh(5)   =rip(fid,'num');
blade.BldFl1Sh(6)   =rip(fid,'num');
blade.BldFl2Sh(2)   =rip(fid,'num');
blade.BldFl2Sh(3)   =rip(fid,'num');
blade.BldFl2Sh(4)   =rip(fid,'num');
blade.BldFl2Sh(5)   =rip(fid,'num');
blade.BldFl2Sh(6)   =rip(fid,'num');
blade.BldEdgSh(2)   =rip(fid,'num');
blade.BldEdgSh(3)   =rip(fid,'num');
blade.BldEdgSh(4)   =rip(fid,'num');
blade.BldEdgSh(5)   =rip(fid,'num');
blade.BldEdgSh(6)   =rip(fid,'num');

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