function bmodes = readBModesMain(input_file);
%READBMODESMAIN  Read a BMODES primary input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   bmodes = ReadBModesMain('input_file_name')
%   Returns a structure containing the data in a BModes main input file.
%   compatible with BModes v3.00.00

if ~exist('input_file','var')
    [fn pn] = uigetfile( ...
        {'*.bmi','BModes input (*.bmi)'; ...
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
bmodes.ffn = input_file;

%% get file description from header
fgetl(fid);
bmodes.title{1,1}=fgetl(fid);
bmodes.title{2,1}=fgetl(fid);


% bmodes.General.beam_type    =rip(fid,'str',{'1','2'});


%% General parameters
fgetl(fid);
bmodes.Echo                 =rip(fid,'TF');
bmodes.beam_type	        =rip(fid,'num');
if (bmodes.beam_type~=1) && (bmodes.beam_type~=2), warning('BModes beam type must be set to 1 or 2');end
bmodes.rot_rpm              =rip(fid,'num');
bmodes.rpm_mult             =rip(fid,'num');
bmodes.radius               =rip(fid,'num');
bmodes.hub_rad              =rip(fid,'num');
bmodes.precone              =rip(fid,'num');
bmodes.bl_thp               =rip(fid,'num');
bmodes.hub_conn             =rip(fid,'num');
bmodes.modepr               =rip(fid,'num');
bmodes.TabDelim             =rip(fid,'TF');
bmodes.mid_node_tw          =rip(fid,'TF');

%% Blade-tip or tower-top mass properties
fgetl(fid);
fgetl(fid);
bmodes.tip_mass             =rip(fid,'num');
bmodes.cm_loc               =rip(fid,'num');
bmodes.cm_axial             =rip(fid,'num');
bmodes.ixx_tip              =rip(fid,'num');
bmodes.iyy_tip              =rip(fid,'num');
bmodes.izz_tip              =rip(fid,'num');
bmodes.ixy_tip              =rip(fid,'num');
bmodes.izx_tip              =rip(fid,'num');
bmodes.iyz_tip              =rip(fid,'num');

%% Distributed-property identifiers
fgetl(fid);
fgetl(fid);
bmodes.id_mat               =rip(fid,'num');
bmodes.sec_props_file       =rip(fid,'qstr');

%% Property scaling factors
fgetl(fid);
fgetl(fid);
bmodes.sec_mass_mult        =rip(fid,'num');
bmodes.flp_iner_mult        =rip(fid,'num');
bmodes.lag_iner_mult        =rip(fid,'num');
bmodes.flp_stff_mult        =rip(fid,'num');
bmodes.edge_stff_mult       =rip(fid,'num');
bmodes.tor_stff_mult        =rip(fid,'num');
bmodes.axial_stff_mult      =rip(fid,'num');
bmodes.cg_offst_mult        =rip(fid,'num');
bmodes.sc_offst_mult        =rip(fid,'num');
bmodes.tc_offst_mult        =rip(fid,'num');

%% Finite element discretization
fgetl(fid);
fgetl(fid);
fgetl(fid);  % nselt is determined from the length of the array that is read in next...
fgetl(fid);
bmodes.el_loc               =rip(fid,'numhlist');

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
    case 'numhlist'  % type: horizontal list of numbers
        % regexp match: beginning of string ^, set of zero or more []*
        %               digits \d whitespace \s or ,.-+
        csv = regexp(line,'[\d.eE+-]+','match');
        param = zeros(1,length(csv));
        for k=1:length(csv)
            param(k) = str2double(csv{k});
        end
end
end
