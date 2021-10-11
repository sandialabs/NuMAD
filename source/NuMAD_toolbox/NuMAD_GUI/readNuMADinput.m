function [station, shearweb, active, ansys, BladeRotation, blade, plot3d, flutterInput] = readNuMADinput(filename)
%READNUMADINPUT  Read a NuMAD project file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   [station, shearweb, active, ansys, BladeRotation, ...
%       blade, plot3d, flutterInput] = readNuMADinput(FILENAME)
%     Read the NuMAD project FILENAME (e.g. 'my_blade.nmd')
%     and return the data structures
%

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2010.12.17  JCB: first draft, based on readMatDB()

%filename = 'pseudoseri8.txt'; clear station;

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

% use find-and-replace to convert certain tokens
filecontents = strrep(filecontents,'genlinex','presweep');
filecontents = strrep(filecontents,'genliney','precurve');

% %%%%%%%%%%%%  PROGRAMMER'S NOTE %%%%%%%%%%%%%
% The following regular expression pattern matches any number of characters
% found between the opening and closing "reference" tags
%pattern = '<reference>(.*)</reference>';
%t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
% try
%     % jcb: is there a better way to extract the contents of t?
%     af.reference = cell2mat(t{1});  
% catch me
%     af.reference = '';
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% get the blade rotation direction (cw or ccw)
pattern = '<rotation>(.*)</rotation>';
t = regexp(filecontents, pattern, 'tokens');
if isempty(t)
    fprintf('assuming ccw blade\n');
    BladeRotation = 'ccw';
else
    BladeRotation = cell2mat(t{1});
end
blade.BladeRotation = BladeRotation;

blade.PresweepRef = struct('method','normal','table',[0 0 0],'pptype','poly');
blade.PrecurveRef = struct('method','shear','table',[0 0 0],'pptype','poly');
% get the Generating Line parameters
% find the start and end of PresweepRef block
blkS = regexp(filecontents, '<presweep>');
blkE = regexp(filecontents, '</presweep>');
if numel(blkS)==1
    % adjust start and end indices to encompass only child elements
    blkS = blkS + length('<presweep>');
    blkE = blkE - 1;
    blkC = filecontents(blkS:blkE);
    
    pattern = '<method>(.*)</method>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PresweepRef.method = cell2mat(t{1});
    
    pattern = '<table>(.*)</table>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PresweepRef.table = cell2mat(textscan(cell2mat(t{1}),'%f %f %f'));
    
    pattern = '<pptype>(.*)</pptype>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PresweepRef.pptype = cell2mat(t{1});
end
% find the start and end of PrecurveRef block
blkS = regexp(filecontents, '<precurve>');
blkE = regexp(filecontents, '</precurve>');
if numel(blkS)==1
    % adjust start and end indices to encompass only child elements
    blkS = blkS + length('<precurve>');
    blkE = blkE - 1;
    blkC = filecontents(blkS:blkE);
    
    pattern = '<method>(.*)</method>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PrecurveRef.method = cell2mat(t{1});
    
    pattern = '<table>(.*)</table>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PrecurveRef.table = cell2mat(textscan(cell2mat(t{1}),'%f %f %f'));
    
    pattern = '<pptype>(.*)</pptype>';
    t = regexp(blkC, pattern, 'tokens');
    blade.PrecurveRef.pptype = cell2mat(t{1});
end

% find the start and end of each station block
blkS = regexp(filecontents, '<station>');
blkE = regexp(filecontents, '</station>');
% adjust start and end indices to encompass only child elements
blkS = blkS + length('<station>');
blkE = blkE - 1;

% process each station block
Nblocks = length(blkS);
for kblk = 1:Nblocks
    blkC = filecontents(blkS(kblk):blkE(kblk));
    
    pattern = '<airfoilname>(.*)</airfoilname>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).AirfoilName = cell2mat(t{1});
    
    pattern = '<tetype>(.*)</tetype>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).TEtype = cell2mat(t{1});
    
    pattern = '<degreestwist>(.*)</degreestwist>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).DegreesTwist = str2double(cell2mat(t{1}));
    
    pattern = '<locationz>(.*)</locationz>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).LocationZ = str2double(cell2mat(t{1}));
    
    pattern = '<xoffset>(.*)</xoffset>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).Xoffset = str2double(cell2mat(t{1}));
    
    pattern = '<aerocenter>(.*)</aerocenter>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).AeroCenter = str2double(cell2mat(t{1}));
    
    pattern = '<chord>(.*)</chord>';
    t = regexp(blkC, pattern, 'tokens');
    station(kblk).Chord = str2double(cell2mat(t{1}));
    
    pattern = '<coords>(.*)</coords>';
    t = regexp(blkC, pattern, 'tokens');
    coord_text = cell2mat(t{1});
    station(kblk).coords = cell2mat(textscan(coord_text,'%f %f'));
    
    pattern = '<delineationpoint>(.*)</delineationpoint>';
    t = regexp(blkC, pattern, 'tokens');
    num_text = cell2mat(t{1});
    values = textscan(num_text,'%f %s %s');
    station(kblk).dp = values{1};
    station(kblk).dptype = values{2};
    station(kblk).dpmaterial = values{3};
    
    pattern = '<surfacematerial>(.*)</surfacematerial>';
    t = regexp(blkC, pattern, 'tokens');
    num_text = cell2mat(t{1});
    values = textscan(num_text,'%s');
    station(kblk).sm = values{1};
end

% The tip station does not need any surfacematerial definition
% (and thus is not in the input file) but not having anything
% in the data structure can create bugs when creating DPs or
% creating a station outboard of the tip station.
%
% This FOR loop populates the surface material list for the tip
% station
for k = 1:numel(station(end).dp)-1
    station(end).sm{k,1} = '**UNSPECIFIED**';
end

% find the start and end of each shear web block
blkS = regexp(filecontents, '<shearweb>');
blkE = regexp(filecontents, '</shearweb>');
% adjust start and end indices to encompass only child elements
blkS = blkS + length('<shearweb>');
blkE = blkE - 1;

% process each shear web block
%shearweb = struct('Material',[],'BeginStation',[],'EndStation',[],'Corner',[]);
shearweb = struct([]);
Nblocks = length(blkS);
for kblk = 1:Nblocks
    blkC = filecontents(blkS(kblk):blkE(kblk));
    
    pattern = '<material>(.*)</material>';
    t = regexp(blkC, pattern, 'tokens');
    shearweb(kblk).Material = cell2mat(t{1});
    
    pattern = '<beginstation>(.*)</beginstation>';
    t = regexp(blkC, pattern, 'tokens');
    shearweb(kblk).BeginStation = str2double(cell2mat(t{1}));
    
    pattern = '<endstation>(.*)</endstation>';
    t = regexp(blkC, pattern, 'tokens');
    shearweb(kblk).EndStation = str2double(cell2mat(t{1}));
    
    pattern = '<corner>(.*)</corner>';
    t = regexp(blkC, pattern, 'tokens');
    num_text = cell2mat(t{1});
    values = textscan(num_text,'%f %f %f %f');
    shearweb(kblk).Corner = cell2mat(values);
end

% find the start and end of the active material block
blkS = regexp(filecontents, '<activematerials>');
blkE = regexp(filecontents, '</activematerials>');
% adjust start and end indices to encompass only child elements
blkS = blkS + length('<activematerials>');
blkE = blkE - 1;

pattern = '<list>(.*)</list>';
blkC = filecontents(blkS(1):blkE(1));
t = regexp(blkC, pattern, 'tokens');
values = textscan(cell2mat(t{1}),'%s');
active.list = values{1};

% pattern = '<colors>(.*)</colors>';
% blkC = filecontents(blkS(1):blkE(1));
% t = regexp(blkC, pattern, 'tokens');
% values = textscan(cell2mat(t{1}),'%f %f %f');
% active.color = cell(size(active.list));
% for k = 1:numel(active.list)
%     active.color{k} = [values{1}(k) values{2}(k) values{3}(k)];
% end



% find the start and end of the ANSYS block
blkS = regexp(filecontents, '<ansys>');
blkE = regexp(filecontents, '</ansys>');
% adjust start and end indices to encompass only child elements
blkS = blkS + length('<ansys>');
blkE = blkE - 1;
blkC = filecontents(blkS(1):blkE(1));

ansys.BoundaryCondition = regexp_str(blkC,'<boundarycondition>(.*)</boundarycondition>');
ansys.ElementSystem = regexp_str(blkC,'<elementsystem>\D*(.*)\D*</elementsystem>');
ansys.MultipleLayerBehavior = regexp_str(blkC,'<multiplelayerbehavior>(.*)</multiplelayerbehavior>','distinct');
ansys.meshing = regexp_str(blkC,'<meshing>(.*)</meshing>');
ansys.smartsize = regexp_dbl(blkC,'<smartsize>(.*)</smartsize>');
ansys.elementsize = regexp_dbl(blkC,'<elementsize>(.*)</elementsize>');
ansys.shell7gen = regexp_dbl(blkC,'<shell7gen>(.*)</shell7gen>');
ansys.dbgen = regexp_dbl(blkC,'<dbgen>(.*)</dbgen>');
fclist = regexp_str(blkC,'<failurecriteria>(.*)</failurecriteria>');
fcopts = {'EMAX','SMAX','TWSI','TWSR','HFIB','HMAT','PFIB','PMAT',...
    'L3FB','L3MT','L4FB','L4MT','USR1','USR2','USR3','USR4',...
    'USR5','USR6','USR7','USR8','USR9'};
ansys.FailureCriteria = cell(numel(fcopts),2);
ansys.FailureCriteria(:,1) = fcopts';
for kfc = 1:numel(fcopts)
    s = regexpi(fclist,fcopts{kfc});
    ansys.FailureCriteria{kfc,2} = ~isempty(s);
end

% get the Plot3D parameters
% find the start and end of plot3d block
blkS = regexp(filecontents, '<plot3d>');
blkE = regexp(filecontents, '</plot3d>');
if numel(blkS)==1
    % adjust start and end indices to encompass only child elements
    blkS = blkS + length('<plot3d>');
    blkE = blkE - 1;
    blkC = filecontents(blkS:blkE);
    
    pattern = '<n_panels>(.*)</n_panels>';
    t = regexp(blkC, pattern, 'tokens');
    plot3d.n_panels = str2double(cell2mat(t{1}));
    
    pattern = '<spacing>(.*)</spacing>';
    t = regexp(blkC, pattern, 'tokens');
    plot3d.spacing = cell2mat(t{1});
    
    pattern = '<interpmethod>(.*)</interpmethod>';
    t = regexp(blkC, pattern, 'tokens');
    plot3d.interpmethod = cell2mat(t{1});
    
    pattern = '<breakpoints>(.*)</breakpoints>';
    t = regexp(blkC, pattern, 'tokens');
    if isempty(t)
        plot3d.breakpoints = [];
    else
        plot3d.breakpoints = cell2mat(textscan(cell2mat(t{1}),'%f'));
    end
    
    pattern = '<newspanloc>(.*)</newspanloc>';
    t = regexp(blkC, pattern, 'tokens');
    if isempty(t)
        plot3d.newspanloc = [];
    else
        plot3d.newspanloc = cell2mat(textscan(cell2mat(t{1}),'%f'));
    end
else
    plot3d.n_panels = 100;
    plot3d.spacing = 'cosine';
    plot3d.interpmethod = 'linear';
    plot3d.breakpoints = [];
    plot3d.newspanloc = [];
    
end

% get the Flutter Tool (BLAST) Input parameters
% find the start and end of flutterinput block
blkS = regexp(filecontents, '<flutterinput>');
blkE = regexp(filecontents, '</flutterinput>');

% adjust start and end indices to encompass only child elements
blkS = blkS + length('<flutterinput>');
blkE = blkE - 1;
blkC = filecontents(blkS:blkE);

flutterInput.fstFile = regexp_str(blkC,'<fstfile>(.*)</fstfile>');
flutterInput.bladeFile = regexp_str(blkC,'<bladefile>(.*)</bladefile>');
flutterInput.aeroFile = regexp_str(blkC,'<aerofile>(.*)</aerofile>');
flutterInput.outFile = regexp_str(blkC,'<outfile>(.*)</outfile>');
if isempty(flutterInput.outFile)
    flutterInput.outFile = 'name.out';
end
pattern = '<lcs>(.*)</lcs>';
t = regexp(blkC, pattern, 'tokens');
if isempty(t) || isempty(t{1}) || all(isspace(t{1}{1}))
    flutterInput.LCS = [];
else
    flutterInput.LCS = cell2mat(textscan(cell2mat(t{1}),'%f'));
end
flutterInput.OmegaStart = regexp_dbl(blkC,'<omegastart>(.*)</omegastart>');
flutterInput.OmegaIncrement = regexp_dbl(blkC,'<omegaincrement>(.*)</omegaincrement>');
flutterInput.OmegaEnd = regexp_dbl(blkC,'<omegaend>(.*)</omegaend>');

end  % end readNuMADinput

function str = regexp_str(input_str,pattern,default_str)
    if ~exist('default_str','var')
        default_str = '';
    end
    t = regexp(input_str, pattern, 'tokens');
    if isempty(t)
        str = default_str;
    else
        str = cell2mat(t{1});
    end
end

function dbl = regexp_dbl(input_str,pattern,default_val)
    if ~exist('default_val','var')
        default_val = [];
    end
    t = regexp(input_str, pattern, 'tokens');
    if isempty(t) || isempty(t{1}) || isempty(t{1}{1});
        dbl = default_val;
    else
        dbl = str2double(cell2mat(t{1}));
    end
end
