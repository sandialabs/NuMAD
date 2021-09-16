function data = readANSYSSections(filename)
% readANSYSSections  Read an ANSYS list of Sections.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Read an ANSYS list of Sections.
% Usage: data = read_sections(filename)
%  where FILENAME is file name string, default 'Sections.txt'
%        DATA is a stucture with fields 'secID' and 'layers'
%        .secID is an array of Section ID Numbers
%        .layers is a cell array of tables with columns:
%          [Layer, Thickness, MatID, Angle]
%  Note: access a table position with syntax data.layers{ksec}(row,col)
%

defaultfn = 'Sections.txt';

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

% seperate the lines of text
filelines = textscan(filecontents,'%s','Delimiter','\n');
filelines = filelines{1};

kline = 1;
ksec = 1;
while true
    t = regexp(filelines{kline},'SECTION ID NUMBER:\s*(\d+)','tokens');
    if isempty(t)
        kline = kline+1;  % go to next line
    else
        % new section found
        SecID = str2double(t{1}{1});  % Section ID Number
        kline = kline+4;  % skip down 4 lines
        t = regexp(filelines{kline},'Number of Layers\s*=\s*(\d+)','tokens');
        NLayers = str2double(t{1}{1}); % Number of Layers in section
        kline = kline+5;  % skip down 5 lines
        LayerTbl = zeros(NLayers,4);
        for k = 1:NLayers
            LayerTbl(k,:) = cell2mat(textscan(filelines{kline},'%f %f %f %f %*f'));
            kline = kline+1;  % go to next line
        end
        data.secID(ksec) = SecID;
        data.layers{ksec} = LayerTbl;
        ksec = ksec + 1;
    end
    
    if kline > numel(filelines)
        break
    end
end
end