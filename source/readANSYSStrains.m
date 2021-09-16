function data = readANSYSStrains(filename)
% readANSYSStrains  Read an ANSYS POST1 element table listing.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% data = readANSYSStrains(filename)
% Read an ANSYS POST1 element table listing.
% Usage: data = read_strains(filename)
%  where FILENAME is file name string, default 'Strains.txt'
%        DATA is 7-column matrix [ELEM, EPELX, EPELY, EPELZ, EPELXY, EPELYZ, EPELXZ]
%

defaultfn = 'Strains.txt';

% hard-code filename if not specified
if ~exist('filename','var')
    filename = defaultfn;
end

%     % user select filename if not specified
%     if ~exist('filename','var') || isempty(filename)
%         [fn,pn] = uigetfile( ...
%             {'*.txt','Text files(*.txt)'; ...
%             '*.*','All files (*.*)'},...
%             'Select ANSYS element list',defaultfn);
%         if isequal(fn,0) || isequal(pn,0)
%             disp('Operation canceled by user.')
%             return;
%         end
%         filename = fullfile(pn,fn);
%     end

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);
%assignin('base','filecontents',filecontents);  %debugging

% process the tables
NCOLS = 7;
data = cell(1,NCOLS);
pat = 'ELEM\s*EPELX\s*EPELY\s*EPELZ\s*EPELXY\s*EPELYZ\s*EPELXZ\s*';
tbl_hdrs = regexp(filecontents,pat);  % find location of table headers
tbl_hdrs(end+1) = numel(filecontents);  % append the file end location
%assignin('base','tbl_hdrs',tbl_hdrs);  %debugging
for kTbl = 1:numel(tbl_hdrs)-1
    tbl = filecontents(tbl_hdrs(kTbl):tbl_hdrs(kTbl+1)-1);  % grab table
    tbl = regexprep(tbl,pat,'');  % remove header by replacement
    data = [data; textscan(tbl,repmat(' %f',1,NCOLS))];
end
data = cell2mat(data);

end