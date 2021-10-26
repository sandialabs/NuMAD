function data = readANSYSElem(filename)
% readANSYSElem  Read an ANSYS list of Elements.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Read an ANSYS list of Elements.
% Usage: data = readANSYSElem(FILENAME)
%  where FILENAME is file name string, default 'Elements.txt'
%        DATA is 3-column matrix [ELEM, MAT, SEC]
%

defaultfn = 'Elements.txt';

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
NCOLS = 14; %14 needed when quadratic elements are used. If linear elements are used
            % Then the last four coulumns are NaN. Alternatively NCOLS
            % could be set to 10 for linear elements and no NaNs would
            % apprear
data = cell(1,NCOLS);
pat = 'ELEM\s*MAT\s*TYP\s*REL\s*ESY\s*SEC\s*NODES\s*';
%pat='*** NOTE ***\s*CP =\s*\d*.\d*\s*TIME=\s*\d*:\d*:\d*\r\s*Use the SLIST command to list section data for element\s*\d*.\s*Section\s*\r\s*data overrides the real constant data.|*** NOTE ***\s*CP =\s*\d*.\d*\s*TIME=\s*\d*:\d*:\d*\r\s*Use the SLIST command to list section data for element\s*\d*.\s*Section data\s*\r\s*overrides the real constant data. ';
tbl_hdrs = regexp(filecontents,pat);  % find location of table headers
tbl_hdrs(end+1) = numel(filecontents);  % append the file end location
%assignin('base','tbl_hdrs',tbl_hdrs);  %debugging
for kTbl = 1:numel(tbl_hdrs)-1
    tbl = filecontents(tbl_hdrs(kTbl):tbl_hdrs(kTbl+1)-1);  % grab table
    tbl = regexprep(tbl,pat,'');  % remove header by replacement
    data = [data; textscan(tbl,repmat(' %f',1,NCOLS))];
end
data = cell2mat(data);

% only columns 1,2,6 are needed now
% Note: textscan can skip fields by using %*f in place of %f on the
% columns to skip
%data = data(:,[1 2 6]);

end