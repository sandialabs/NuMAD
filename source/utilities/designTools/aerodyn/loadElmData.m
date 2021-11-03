function out = loadElmData(ffn)
%LOADELMDATA  Read AeroDyn element data into structure
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   Syntax:
%     elm = loadElmData
%     elm = loadElmData('file.elm')
%
%   The output structure has the following fields:
%          ffn - the full file name of the data file
%         info - version and simulation runtime info
%      summary - a table summarizing the column labels and units
%         data - data matrix
%       PrnElm - which elements had 'PrnElm' set in AeroDyn input
%        Alpha - columns that have angle-of-attack data
%         ...
%        ReNum - columns that have Reynolds number data
%
%   Fields 'Alpha' through 'ReNum' are intended for indexing:
%
%     row=1; plot(elm.PrnElm,elm.data(row,elm.ForcN));
%
%   See also loadOutData, plottingUtility

ext = '.elm';  % this script is intended for files with this extension
wext = ['*',ext];  % extension with wildcard
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle the input argument (4 options):
% 1: no input argument
if ~exist('ffn','var')
    [fn,pn] = uigetfile(wext);
    ffn = fullfile(pn,fn);
% 2: input argument is directory
elseif exist(ffn,'dir')
    [fn,pn] = uigetfile(fullfile(ffn,wext));
    ffn = fullfile(pn,fn);
% 3: input argument is file
elseif exist(ffn,'file')
    [pn,fn,ext] = fileparts(ffn);
    if isempty(pn) || pn(1)=='.'
        pn = fullfile(pwd,pn);
        ffn = fullfile(pn,[fn,ext]);
    end
% 4: input argument has wildcard or is bad filename
else
    if regexp(ffn,'[*]')
        [fn,pn] = uigetfile(ffn);
        ffn = fullfile(pn,fn);
    else
        error('Input argument not valid: %s\n',ffn);
    end
end

% fn will be 0 if File Open dialog was cancelled
if ~ischar(fn)
    out=[];
    return;  % exit early if no file selected
end

%%%%%%%%%%%%%%%%%%%%%%
% Open the data file %
fid=fopen(ffn,'r');
if fid==-1
    error('Could not open file: %s\n',ffn);
    return;
end

%%%%%%%%%%%%%%%%%%%%%%
% Read in the header %
info = fgetl(fid); % version of AeroDyn; version of FAST; date-time of simulation
line=fgetl(fid); labels=strtrim(regexp(line,'[^\t]+','match'));
line=fgetl(fid); units =strtrim(regexp(line,'[^\t]+','match'));
nc=numel(labels);  % number of columns
columns = num2cell(1:nc);

%%%%%%%%%%%%%%%%%%%%
% Read in the data %
A=textscan(fid,repmat('%f',1,nc),...  % assumes all rows have nc columns
              'Delimiter','\t');

fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create output structure %
out.ffn = ffn;
out.info = info;
out.summary = [columns' labels' units'];
out.data = cell2mat(A);

% indicate which elements have data 
cols = strmatch('Alpha',labels);
out.PrnElm = str2double(regexp(strcat(labels{cols}),'\d+','match'));

% find columns that correspond to each output (see AeroDyn manual)
out.Alpha   = strmatch('Alpha',labels)';    % angle-of-attack
out.DynPres = strmatch('DynPres',labels)';  % dynamic pressure
out.CLift   = strmatch('CLift',labels)';    % lift coefficient
out.CDrag   = strmatch('CDrag',labels)';    % drag coefficient
out.CNorm   = strmatch('CNorm',labels)';    % normal coeff
out.CTang   = strmatch('CTang',labels)';    % tangent coeff
out.CMom    = strmatch('CMom',labels)';     % pitching moment coeff (not always available)
out.Pitch   = strmatch('Pitch',labels)';    % local pitch (includes structural twist)
out.AxInd   = strmatch('AxInd',labels)';    % axial induction
out.TanInd  = strmatch('TanInd',labels)';   % tangential induction
out.ForcN   = strmatch('ForcN',labels)';    % force normal to rotor plane
out.ForcT   = strmatch('ForcT',labels)';    % force tangent to rotor plane
out.Pmomt   = strmatch('Pmomt',labels)';    % pitching moment (not always available)
out.ReNum   = strmatch('ReNum',labels)';    % Reynolds number
