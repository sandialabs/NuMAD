function data = readANSYSSectionsTempEC(filename)
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

%% Determine if reading ABD matrix info (ie if Details = FULL in ANSYS SLIST commend
kline = 1;
while true
    t = regexp(filelines{kline},'LIST SECTION ID SETS');
    if t ==1
        t = regexp(filelines{kline+1},'Details = FULL');
        if t==1
            abdFlag=1;
            break
        else
            abdFlag=0;
            break
        end
    else
        kline = kline+1;  % go to next line
    end
end

    kline = 1;
    ksec = 1;
if abdFlag==0
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
elseif abdFlag==1
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
            
            t = regexp(filelines{kline+1},'Total Thickness\s*=\s*(\d+.\d+)','tokens');
            hLaminate = str2double(t{1}{1}); % Total thickness
            
            kline = kline+5;  % skip down 5 lines
            LayerTbl = zeros(NLayers,4);
            for k = 1:NLayers
                LayerTbl(k,:) = cell2mat(textscan(filelines{kline},'%f %f %f %f %*f'));
                kline = kline+1;  % go to next line
            end
            data.totalThick(ksec) = hLaminate;
            data.secID(ksec) = SecID;
            data.layers{ksec} = LayerTbl;
            
            %Plate Stiffness
            if SecID==1
                kline = kline+4;%+3;  % skip down 7 lines
            else
                kline = kline+4;  % skip down 4 lines
            end
            ABD = zeros(6);
            for k = 1:6
                nonzeros=0;
                while ~nonzeros
                    try
                        ABD(k,:) = cell2mat(textscan(filelines{kline},'%f %f %f %f %f %f'));
                        nonzeros=1;
                    catch
                        nonzeros=0;
                    end
                    kline = kline+1;  % go to next line
                end
            end
            data.ABD{ksec} = ABD;
            
            %Transverse shear stiffness matrix
            kline = kline+4;  % skip down 5 lines
            transShearStiff=zeros(2);
            transShearStiff(1,:) = cell2mat(textscan(filelines{kline},'%f %f'));
            transShearStiff(2,:) = cell2mat(textscan(filelines{kline+2},'%f %f'));
            data.transShearStiff{ksec} = transShearStiff;
            
            %Transverse shear stiffness correction factors
            kline = kline+7;  % skip down 7 lines
            transShearCorrection=zeros(2);
            transShearCorrection(1,:) = cell2mat(textscan(filelines{kline},'%f %f'));
            transShearCorrection(2,:) = cell2mat(textscan(filelines{kline+2},'%f %f'));
            data.transShearCorrection{ksec} = transShearCorrection;
            
            %Shell reference surface offset
            kline = kline+4;  % skip down 7 lines
            if contains(filelines{kline},'TOP')
                offset='TOP';
            elseif contains(filelines{kline},'MID')
                offset='MID';
            elseif contains(filelines{kline},'BOT')
                offset='BOT';
            else
                error('Shell Offset not found')
            end
            data.shellOffset{ksec} = offset;
            
            ksec = ksec + 1;
        end

        if kline > numel(filelines)
            break
        end
    end
    
    
    %Overright the thickness information from the shell7 file since SLIST
    %rounds for small layers
    filename='shell7.src';
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
        t = regexp(filelines{kline},'sectype,(\d+),shell','tokens');
        if isempty(t)
            kline = kline+1;  % go to next line
        else
            % new section found
            SecID = str2double(t{1}{1});  % Section ID Number
            index=find(data.secID==SecID);
            if isempty(index)
                disp('Section ID in shell7 not found in model. Exiting loop...')
                break
            end
            [nLayers,~]=size(data.layers{index});
            for iLayer=1:nLayers
                kline = kline+1;  % go to next line
                
                t=regexp(filelines{kline},'secdata,(\d+.\d+),(.*+)','tokens');
                data.layers{index}(iLayer,2)=str2double(t{1}{1});
            end
            kline = kline+1;  % go to next line
            
        end
        if kline > numel(filelines)
            break
        end
    end 
end