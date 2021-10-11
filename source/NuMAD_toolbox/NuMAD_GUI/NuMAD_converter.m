function NuMAD_converter
%NUMAD_CONVERTER  Convert a NuMAD project saved in version 2010.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   NuMAD_converter takes a NuMAD project saved in version 2010 and
%   converts it to the new format.
%
%   After the initial welcome and instruction dialog, the user is prompted
%   to select a Sdata*.nmd from a legacy NuMAD project.
%
%   After conversion is complete, the user is prompted to save the
%   converted NuMAD file.


h=helpdlg('Welcome to the legacy NuMAD converter.  Please select any file (e.g. Sdata1.nmd) in a NuMAD project saved in version 2010.','NuMAD converter');
waitfor(h);  % wait until the user closes the warning dialog

if ispc
    userhome = getenv('USERPROFILE');
elseif isunix
    userhome = getenv('HOME');
else
    userhome = '';
end

[fn,pn] = uigetfile( ...
    {'Sdata*.nmd', 'Legacy NuMAD (Sdata*.nmd)';...
    '*.*',   'All files (*.*)'},...
    'Open a NuMAD blade model',userhome);

if isequal(fn,0)
    % the user canceled file selection
    return
end

file1 = 'Sdata1.nmd';
fid = fopen(fullfile(pn,file1),'rt');
if (fid == -1)
    errordlg(sprintf('Could not open file "%s"',file1),'Error');
    return
end
% read the entire file contents
txt1 = fread(fid,inf,'uint8=>char')';
fclose(fid);

file3 = 'Sdata3.nmd';
fid = fopen(fullfile(pn,file3),'rt');
if (fid == -1)
    errordlg(sprintf('Could not open file "%s"',file3),'Error');
end
% read the entire file contents
txt3 = fread(fid,inf,'uint8=>char')';
fclose(fid);

%============== BEGIN PARSING SData1.nmd ====================
t = regexp(txt1,'UnitSystem\s+(\S*)','tokens');
old.UnitSystem = t{1}{1};

t = regexp(txt1,'BoundaryCondition\s+(\S*)','tokens');
old.BoundaryCondition = t{1}{1};

t = regexp(txt1,'TotalStations\s+(\S*)','tokens');
old.TotalStations = str2double(t{1}{1});

t = regexp(txt1,'MeshDensity\s+(\S*)','tokens');
old.MeshDensity = str2double(t{1}{1});

t = regexp(txt1,'ElementSize\s+(\S*)','tokens');
old.ElementSize = str2double(t{1}{1});

t = regexp(txt1,'BladeRotation\s+(\S*)','tokens');
old.BladeRotation = t{1}{1};

t = regexp(txt1,'ElementSystem\s+(\S*)','tokens');
old.ElementSystem = t{1}{1};

t = regexp(txt1,'Shell7Check\s+(\S*)','tokens');
old.Shell7Check = t{1}{1};

t = regexp(txt1,'DBcheck\s+(\S*)','tokens');
old.DBcheck = t{1}{1};

t = regexp(txt1,'PreCompcheck\s+(\S*)','tokens');
old.PreCompcheck = t{1}{1};

t = regexp(txt1,'AutoGen\s+(\S*)','tokens');
old.AutoGen = t{1}{1};

t = regexp(txt1,'BeamElements\s+(\S*)','tokens');
if ~isempty(t)
    old.BeamElements = str2double(t{1}{1});
    
    old.BeamNodeLocation = zeros(1,old.BeamElements+1);
    for k = 1:numel(old.BeamNodeLocation)
        pat = sprintf('BeamNodeLocation\\(%d\\)\\s+(\\S*)',k);
        t = regexp(txt1,pat,'tokens');
        old.BeamNodeLocation(k) = str2double(t{1}{1});
    end
end

t = regexp(txt1,'SW\(TotalNumber\)\s+(\S*)','tokens');
old.SW_TotalNumber = str2double(t{1}{1});

t = regexp(txt1,'ActiveList\s+[{]*([^{}]*)[}]*','tokens');
str = t{1}{1};
old.ActiveList = regexp(str,'\S+','match');

t = regexp(txt1,'SWcornerFlag\s+(\S*)','tokens');
old.SWcornerFlag = t{1}{1};
%============== END PARSING SData1.nmd ====================




MaterialList = {};
%============== BEGIN PARSING SData3.nmd ====================
for ksw = 1:old.SW_TotalNumber
    pat = sprintf('SW\\(%d,MaterialName\\)\\s+(\\S*)',ksw);
    t = regexp(txt3,pat,'tokens');
    old.SW(ksw).Material = t{1}{1};
    MaterialList{end+1} = t{1}{1};  %#ok
    
    pat = sprintf('SW\\(%d,BeginStation\\)\\s+(\\S*)',ksw);
    t = regexp(txt3,pat,'tokens');
    old.SW(ksw).BeginStation = str2double(t{1}{1});
    
    pat = sprintf('SW\\(%d,EndStation\\)\\s+(\\S*)',ksw);
    t = regexp(txt3,pat,'tokens');
    old.SW(ksw).EndStation = str2double(t{1}{1});
    
    pat = sprintf('SW\\(%d,Corner,(?<index>\\d+)\\)\\s+(?<value>\\S*)',ksw);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.SW(ksw).Corner(str2double(n(k).index)) = str2double(n(k).value);
    end
end

for kstn = 1:old.TotalStations
    pat = sprintf('AirfoilName\\(%d\\)\\s+"([^"]*)"',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).AirfoilName = t{1}{1};
    
    pat = sprintf('TEtype\\(%d\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).TEtype = t{1}{1};
    
    pat = sprintf('DegreesTwist\\(%d\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).DegreesTwist = str2double(t{1}{1});
    
    pat = sprintf('LocationZ\\(%d\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).LocationZ = str2double(t{1}{1});
    
    pat = sprintf('Xoffset\\(%d\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).Xoffset = str2double(t{1}{1});
    
    pat = sprintf('AeroCenter\\(%d\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).AeroCenter = str2double(t{1}{1});
    
    pat = sprintf('BSP\\(%d,LE\\)\\s+(\\S*)',kstn);
    t = regexp(txt3,pat,'tokens');
    old.station(kstn).BSP.LE = str2double(t{1}{1});
    
    pat = sprintf('BSP\\(%d,(?<index>\\d+),x\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.station(kstn).BSP.x(str2double(n(k).index),1) = str2double(n(k).value);
    end
    
    % calculate the chord length
    old.station(kstn).Chord = max(old.station(kstn).BSP.x) - min(old.station(kstn).BSP.x);
    
    pat = sprintf('BSP\\(%d,(?<index>\\d+),y\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.station(kstn).BSP.y(str2double(n(k).index),1) = str2double(n(k).value);
    end
    
    pat = sprintf('BSP\\(%d,(?<index>\\d+),dp\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.station(kstn).BSP.dp{str2double(n(k).index),1} = n(k).value;
    end
    
    N = numel(n);  % size of x, y, and dp arrays
    old.station(kstn).BSP.perChord = cell(N,1);
    pat = sprintf('BSP\\(%d,(?<index>\\d+),perChord\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.station(kstn).BSP.perChord{str2double(n(k).index)} = str2double(n(k).value);
    end
    
    old.station(kstn).BSP.surface = cell(N,1);
    pat = sprintf('BSP\\(%d,(?<index>\\d+),surface\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        old.station(kstn).BSP.surface{str2double(n(k).index)} = n(k).value;
    end
end

for kstn = 1:old.TotalStations-1
    pat = sprintf('SMN\\(%d,(?<index>\\d+)\\)\\s+(?<value>\\S*)',kstn);
    n = regexp(txt3,pat,'names');
    for k=1:numel(n)
        if isequal(n(k).value,'undefined')
            old.station(kstn).SMN{str2double(n(k).index),1} = '**UNSPECIFIED**';
            MaterialList{end+1} = '**UNSPECIFIED**';  %#ok
        else
            old.station(kstn).SMN{str2double(n(k).index),1} = n(k).value;
            MaterialList{end+1} = n(k).value;  %#ok
        end
    end 
end
%============== END PARSING SData3.nmd ====================




%============== BEGIN CONVERTING FORMAT ===================
new.ReleaseVersion = 'v2.0.1';
switch old.BoundaryCondition
    case '1'
        new.ansys.BoundaryCondition = 'cantilevered';
    case '2'
        new.ansys.BoundaryCondition = 'freefree';
end

if old.MeshDensity == 0
    new.ansys.meshing = 'elementsize';
    new.ansys.elementsize = old.ElementSize;
    new.ansys.smartsize = 5;
else
    new.ansys.meshing = 'smartsize';
    new.ansys.smartsize = old.MeshDensity;
    new.ansys.elementsize = 0.05;
end
new.ansys.ElementSystem = old.ElementSystem;

new.ansys.shell7gen = isequal(old.Shell7Check,'1');
new.ansys.dbgen = isequal(old.DBcheck,'1');


new.blade.BladeRotation = old.BladeRotation;
new.blade.PresweepRef = struct('method','normal','table',[0 0 0],'pptype','poly');
new.blade.PrecurveRef = struct('method','shear','table',[0 0 0],'pptype','poly');

activecolors={[255   0   0]/255;  % red
    [  0 102 255]/255;  % blue
    [  0 204   0]/255;  % green
    [204 102   0]/255;  % brown
    [255 255   0]/255;  % yellow
    [  0 153 255]/255;  % sky blue
    [255 153   0]/255;  % orange
    [  0 153   0]/255;  % dark green
    [204 102 255]/255;  % purple
    [204 153   0]/255;  % light brown
    [  0 255 255]/255;  % cyan
    [255  51   0]/255;  % light red
    [  0 255   0]/255;  % bright green
    [204   0   0]/255;  % dark red
    [255 204 102]/255;  % tan
    [  0 102   0]/255;  % forest green
    [255   0 255]/255;  % magenta
    [102 102 102]/255;  % gray - 40%
    [255   0 102]/255;  % name? (pinkish red)
    [102 153   0]/255}; % olive green
%new.active.list = unique([{'**UNSPECIFIED**'}; transpose(old.ActiveList)]);
new.active.list = unique([{'**UNSPECIFIED**'}; transpose(old.ActiveList); transpose(MaterialList)]);
kclr = rem((1:numel(new.active.list))-1,numel(activecolors))+1;
new.active.color = activecolors(kclr);

if old.SW_TotalNumber > 0
    new.shearweb = old.SW;
else
    new.shearweb = struct([]);
end

[afdb aflist] = readAirfoilDB(fullfile(pn,'airfoils'));
new.station = rmfield(old.station,{'BSP','SMN'});
for k = 1:numel(old.station)
    n = strcmp(new.station(k).AirfoilName,aflist);
    new.station(k).coords = afdb(n).coords;
%     new.station(k).LE = afdb(n).LE;
%     new.station(k).s = afdb(n).s;
%     new.station(k).c = afdb(n).c;
    kclr = cellfun(@isempty,old.station(k).BSP.perChord);
    surface = old.station(k).BSP.surface(~kclr);
    onLP = strcmp(surface,'LP');
    perChord = cell2mat(old.station(k).BSP.perChord(~kclr));
    
    new.station(k).dp = [-1.0; perChord/100 .* (-1 + 2*onLP); 1.0];
    new.station(k).dptype = [{'single'}; lower(old.station(k).BSP.dp(~kclr)); {'single'}];
    new.station(k).dpmaterial = cell(numel(new.station(k).dp),1);
    new.station(k).sm = old.station(k).SMN;
    for kdp = 1:numel(new.station(k).dp)
        switch new.station(k).dptype{kdp}
            case {'flare','hourglass'}
                new.station(k).dpmaterial(kdp) = new.station(k).sm(kdp);
                new.station(k).sm(kdp) = [];
            otherwise
                new.station(k).dpmaterial(kdp) = {''};
        end
    end
% JCB:  this was a possible fix to the problem with DPs at the tip
%       station, however, it seems better to fix in the read/write
%       NuMADinput scipts
%     if k == numel(old.station)
%         for kdp = 1:numel(new.station(k).dp)-1
%             new.station(k).sm{kdp,1} = '**UNSPECIFIED**';
%         end
%     end
end

% DEBUGGING
%assignin('base','new',new);
%assignin('base','old',old);

% prepare to save the new file
t = regexp(pn,'\\([^\\]*)\\$','tokens');  % grab the project directory name
fn = [t{1}{1} '.nmd'];
filename = fullfile(pn,fn);

% have the user select a filename
h=helpdlg('Please enter the new project/file name.  Overwriting the original Sdata*.nmd is NOT recommended.','NuMAD converter');
waitfor(h);  % wait until the user closes the warning dialog
[fn,pn] = uiputfile( ...
    {'*.nmd', 'NuMAD blade (*.nmd)';...
    '*.*',   'All files (*.*)'},...
    'Save As',filename);
if isequal(fn,0)
    % the user canceled file selection
    return
end
[~,~,ext] = fileparts(fn);
if isempty(ext)
    % add the NuMAD extension if it is missing
    fn = [fn '.nmd'];
end
filename = fullfile(pn,fn);

% write the new NuMAD input file
writeNuMADinput(new,filename);
msgbox(sprintf('%s written to %s',fn,pn),'Notification');


end %end NuMAD_converter()

