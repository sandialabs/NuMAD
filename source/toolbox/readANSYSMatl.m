function mpdata = readANSYSMatl(filename,varargin)
%readANSYSMatl  Read ANSYS list of material properties
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   output = readANSYSMatl(filename)
%
%   Default input file name is Materials.txt
%
%   Output is a structure with material ID's and corresponding properties.
%

if ~exist('filename','var')
    filename = 'Materials.txt';
end

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);
%assignin('base','filecontents',filecontents);

pat = 'MPDATA,(?<prop>\w*)\s*,\s*(?<matnum>\d+),\s*(?<stloc>\d+),\s*(?<value>[^,]+),';
data = regexp(filecontents,pat,'names');
%assignin('base','data',data);

for k = 1:numel(data)
    matnum = str2double(data(k).matnum);
    propname = data(k).prop;
    value = str2double(data(k).value);
    mpdata(matnum).(propname) = value;
end
%assignin('base','mpdata',mpdata);


%Add failure criterial properties

if ~isempty(varargin)
    filename = varargin{1};
    % Open the file and read the entire contents
    fid = fopen(filename);
    if (fid == -1)
        error('Could not open file "%s"',filename);
    end
    filecontents = fread(fid,inf,'uint8=>char')';
    fclose(fid);
    

    % seperate the lines of text
    filelines = textscan(filecontents,'%s','Delimiter','\n');
    filelines = filelines{1};

    kline = 1;
    ksec = 1;
    while true
        t = regexp(filelines{kline},'FC LIMIT .FCLI. Table For Material\s*(\d+)','tokens');
        if isempty(t)
            kline = kline+1;  % go to next line
        else
            % new material found
            matnum = str2double(t{1}{1});  % Material ID Number
            kline = kline+5;  % skip down 4 lines
            pat = '(?<prop>\w*)\s*(?<value>[^,]+)';
            
            while ~isempty(filelines{kline}) %Read in the next few lines untill an empty string is encountered
                data = regexp(filelines{kline},pat,'names');
                value = str2double(data.value);
                if isnan(value) % (?<prop>\w*) in pat does not work for GI/GII returns Nan for value
                    propname = 'g'; %g is the ratio of GI/GII
                    data = regexp(filelines{kline},'(?<prop>\w*\W\w*)\s*(?<value>[^,]+)','names'); 
                    
                else
                    propname = data.prop;
                end
                value = str2double(data.value);
                mpdata(matnum).(propname) = value;
                kline=kline+1;
            end

        end

        if kline > numel(filelines)
            break
        end
        
    end
end

end