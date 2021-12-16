 function BladeDef_to_NuMADfile(blade,numad_file,matdb_file,airfoils_path)
%BladeDef_to_NuMADfile   Construct NuMAD files from BladeDef.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%
%   Usage:
%     BladeDef_to_NuMADfile(blade,'BladeDef_test.nmd','MatDBsi.txt');
%
%     BladeDef_to_NuMADfile(blade,'folder/BladeDef_test.nmd',...
%                                 'folder/MatDBsi.txt','folder/airfoils/');
%
%   See also xlsBlade, BladeDef

[nmd_pathstr, nmd_filename] = fileparts(numad_file);
matdb_pathstr = fileparts(matdb_file);

if ~isempty(nmd_pathstr) && exist(nmd_pathstr,'dir')~=7
    error('Path does not exist: %s\n',nmd_pathstr);
end

if ~isempty(matdb_pathstr) && exist(matdb_pathstr,'dir')~=7
    error('Path does not exist: %s\n',nmd_pathstr);
end

if ~exist('airfoils_path','var') || isempty(airfoils_path)
    AIRFOIL_DIR_EXISTS = false;
elseif exist(airfoils_path,'dir')==7
    AIRFOIL_DIR_EXISTS = true;
else
    error('Path does not exist: %s\n',airfoils_path);
end

mm_to_m = 1e-3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Build the required data structures   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data.ReleaseVersion = 'v2.0.1';
% blade rotation =====================================================
if blade.rotorspin == 1
    data.blade.BladeRotation = 'cw';
else
    data.blade.BladeRotation = 'ccw';
end
% presweep ===========================================================
if all(blade.isweep==0)
    table = [0, 0, nan];
else
    N = length(blade.ispan);
    table = [blade.ispan(:), blade.isweep(:), nan(N,1)];
end
data.blade.PresweepRef.method = 'shear';
data.blade.PresweepRef.table = table;
data.blade.PresweepRef.pptype = 'spline';
% precurve ===========================================================
if all(blade.iprebend==0)
    table = [0, 0, nan];
else
    N = length(blade.ispan);
    table = [blade.ispan(:), blade.iprebend(:), nan(N,1)];
end
data.blade.PrecurveRef.method = 'shear';
data.blade.PrecurveRef.table = table;
data.blade.PrecurveRef.pptype = 'spline';

% stations ===========================================================
data.station = struct('AirfoilName' ,'',...
    'TEtype'      ,'',...
    'DegreesTwist',[],...
    'LocationZ'   ,[],...
    'Xoffset'     ,[],...
    'AeroCenter'  ,[],...
    'Chord'       ,[],...
    'coords'      ,[],...
    'dp'          ,[],...
    'dptype'      ,cell(0,1),...
    'dpmaterial'  ,cell(0,1),...
    'sm'          ,cell(0,1));
af = AirfoilDef();
time_stamp = datestr(now);
N = length(blade.ispan);
M = size(blade.keycpos,1);
for k=1:N
    data.station(k).AirfoilName  = sprintf('Interp_%06d',fix(1e3*blade.ispan(k)));
    data.station(k).LocationZ    = blade.ispan(k);
    data.station(k).DegreesTwist = blade.idegreestwist(k);
    data.station(k).Xoffset      = blade.ichordoffset(k) + blade.naturaloffset*blade.xoffset(k);
    data.station(k).AeroCenter   = blade.iaerocenter(k);
    data.station(k).Chord        = blade.ichord(k);
    data.station(k).TEtype       = cell2mat(blade.getprofileTEtype(k));
    data.station(k).coords       = blade.downsampleProfile(k,160); % reduce to about 160 points

    for j=1:M
        data.station(k).dp(j,1)         = blade.keycpos(j,k);
        data.station(k).dptype{j,1}     = 'single';
        data.station(k).dpmaterial{j,1} = '';
    end
    if AIRFOIL_DIR_EXISTS
        af.name = data.station(k).AirfoilName;
        af.reference = sprintf('Airfoil "%s" created on %s for blade project %s',...
            af.name, time_stamp, nmd_filename);
        af.coordinates = data.station(k).coords;
        af_filename = sprintf('%s.txt',af.name);
        af.writeAirfoil(fullfile(airfoils_path,af_filename),'numad');
    end
end
for k=1:N-1
    for j=1:M-1
        data.station(k).sm{j,1} = blade.stacks(j,k).name;
    end
end
% add DPs needed for shear webs
N = length(blade.ispan);
webindices = blade.webindices;
DP_ADDED = false;
for kw=1:numel(blade.webcpos)
    if isnan(webindices{kw}(1))   % HP side
        for k=1:N
            data.station(k).dp(end+1,1) = blade.webcpos{kw}(1,k);
            data.station(k).dptype{end+1,1} = 'single';
            data.station(k).dpmaterial{end+1,1} = '';
        end
        webindices{kw}(1) = size(data.station(1).dp,1);
        DP_ADDED = true;
    end
    if isnan(webindices{kw}(2))   % LP side
        for k=1:N
            data.station(k).dp(end+1,1) = blade.webcpos{kw}(2,k);
            data.station(k).dptype{end+1,1} = 'single';
            data.station(k).dpmaterial{end+1,1} = '';
        end
        webindices{kw}(2) = size(data.station(1).dp,1);
        DP_ADDED = true;
    end
end

if DP_ADDED
    for k=1:N
        [data.station(k).dp, sortorder] = sort(data.station(k).dp);
        if k>1 && ~isequal(prev_sortorder,sortorder)
%             keyboard%ble
% ble            error('Error placing shear web.  Must stay within same region of blade (e.g., TE-Panel)');
        end
        data.station(k).dptype     = data.station(k).dptype(sortorder);
        data.station(k).dpmaterial = data.station(k).dpmaterial(sortorder);
        prev_dp = 0; ksm = 0;
        sm = data.station(k).sm;
        if (k ~= N)
            for j=1:length(sortorder)-1  % for each DP except last
                this_dp = sortorder(j);  % get the DP index
                if (this_dp-prev_dp) == 1
                    % continue as normal
                    prev_dp = this_dp;
                    ksm = ksm + 1;
                else
                    % DP added here: duplicate previous sm
                    % (don't change ksm)
                end
                % ble <<<<<<<<<<<<<<<<<
                if ksm == 0
                    keyboard
                end
                % ble >>>>>>>>>>>>>>>>>>
                assert(ksm~=0,'Unexpected error: shear web before first DP?');
                data.station(k).sm{j,1} = sm{ksm,1};
            end
        end
        if k<4 % ble
            prev_sortorder = sortorder;
        end
    end
    for kw=1:numel(blade.webcpos)
%         this_dp=find(webindices{kw}(1)==sortorder,1);
%         webindices{kw}(1) = this_dp;
%         this_dp=find(webindices{kw}(2)==sortorder,1);
%         webindices{kw}(2) = this_dp;
        % ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        this_dp=find(webindices{kw}(1)==prev_sortorder,1);
        webindices{kw}(1) = this_dp;
        this_dp=find(webindices{kw}(2)==prev_sortorder,1);
        webindices{kw}(2) = this_dp;
        % ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    end
% keyboard%ble
end
% shearwebs ==========================================================
if numel(blade.swstacks) > 0
    data.shearweb = struct('Material','','BeginStation',[],'EndStation',[],'Corner',[]);
else
    data.shearweb = [];
end
N = length(blade.ispan);
j = 1;
for kw=1:numel(blade.swstacks)
    % The "corner" indices of the shearweb should be the same
    % for every ply in the shearweb.
    ind = webindices{kw};
    for k=1:N-1
        if ~isempty(blade.swstacks{kw}(k).plygroups)
            data.shearweb(j).Material     = blade.swstacks{kw}(k).name;
            data.shearweb(j).BeginStation = k;
            data.shearweb(j).EndStation   = k+1;
            data.shearweb(j).Corner       = ind([2,1,1,2]) - 1;  % dp number is offset by 1 in NuMAD
            j = j + 1;
        end
    end
end

% active materials list - leave empty
data.active.list = {'**UNSPECIFIED**'};
% ansys ==============================================================
data.ansys.BoundaryCondition     = 'cantilevered';
data.ansys.ElementSystem         = '181';
data.ansys.MultipleLayerBehavior = 'distinct';
data.ansys.meshing               = 'elementsize';
data.ansys.smartsize             = 5;
% make the element size 5 percent of max chord
data.ansys.elementsize           = 0.05*max(blade.chord);
data.ansys.shell7gen             = 1;
data.ansys.dbgen                 = 1;
% failure criteria struct-field not required
%     data.ansys.FailureCriteria       = {};

% plot3d struct-field not required ===================================

% flutterInput struct-field not required =============================

% prepare material database ==========================================
matdb = struct('type'     ,'',...
    'name'     ,'',...
    'reference','',...
    'dens'     ,[],...
    'nuxy'     ,[],...
    'ex'       ,[],...
    'ey'       ,[],...
    'ez'       ,[],...
    'gxy'      ,[],...
    'gyz'      ,[],...
    'gxz'      ,[],...
    'prxy'     ,[],...
    'pryz'     ,[],...
    'prxz'     ,[],...
    'xten'     ,[],...
    'xcmp'     ,[],...
    'yten'     ,[],...
    'ycmp'     ,[],...
    'zten'     ,[],...
    'zcmp'     ,[],...
    'xy'       ,[],...
    'yz'       ,[],...
    'xz'       ,[],...
    'xycp'     ,[],...
    'yzcp'     ,[],...
    'xzcp'     ,[],...
    'xzit'     ,[],...
    'xzic'     ,[],...
    'yzit'     ,[],...
    'yzic'     ,[],...
    'g1g2'     ,[],...
    'etal'     ,[],...
    'etat'     ,[],...
    'alp0'     ,[],...
    'thicknessType',[],...
    'uniqueLayers' ,[],...
    'symmetryType' ,[],...
    'layer'        ,[]);
for k=1:length(blade.materials)
    matdb(k).name = blade.materials(k).name;
    matdb(k).type = blade.materials(k).type;
    matdb(k).ex   = blade.materials(k).ex;
    matdb(k).ey   = blade.materials(k).ey;
    matdb(k).ez   = blade.materials(k).ez;
    matdb(k).gxy  = blade.materials(k).gxy;
    matdb(k).gyz  = blade.materials(k).gyz;
    matdb(k).gxz  = blade.materials(k).gxz;
    if strcmpi(matdb(k).type,'isotropic')
        matdb(k).nuxy = blade.materials(k).prxy;
    else
        matdb(k).prxy = blade.materials(k).prxy;
        matdb(k).pryz = blade.materials(k).pryz;
        matdb(k).prxz = blade.materials(k).prxz;
    end
    matdb(k).dens = blade.materials(k).density;
    
    if length(blade.materials(k).uts)==1 % isotropic
    matdb(k).xten   = blade.materials(k).uts(1);
    matdb(k).yten   = blade.materials(k).uts(1);
    matdb(k).zten   = blade.materials(k).uts(1);
    matdb(k).xcmp   = blade.materials(k).ucs(1);
    matdb(k).ycmp   = blade.materials(k).ucs(1);
    matdb(k).zcmp   = blade.materials(k).ucs(1);
    try 
        matdb(k).xy     = blade.materials(k).uss(1);
        matdb(k).yz     = blade.materials(k).uss(1);
        matdb(k).xz     = blade.materials(k).uss(1);
    catch
        %Do Nothing
    end
    else % orthotropic
    matdb(k).xten   = blade.materials(k).uts(1);
    matdb(k).yten   = blade.materials(k).uts(2);
    matdb(k).zten   = blade.materials(k).uts(3);
    matdb(k).xcmp   = blade.materials(k).ucs(1);
    matdb(k).ycmp   = blade.materials(k).ucs(2);
    matdb(k).zcmp   = blade.materials(k).ucs(3);
    matdb(k).xy     = blade.materials(k).uss(1);
    matdb(k).yz     = blade.materials(k).uss(2);
    matdb(k).xz     = blade.materials(k).uss(3);   
    end
    matdb(k).xzit   = blade.materials(k).xzit;
    matdb(k).xzic   = blade.materials(k).xzic;
    matdb(k).yzit   = blade.materials(k).yzit;
    matdb(k).yzic   = blade.materials(k).yzic;
    matdb(k).g1g2   = blade.materials(k).g1g2;
    matdb(k).alp0   = blade.materials(k).alp0;
    matdb(k).etat   = blade.materials(k).etat;
    matdb(k).etal   = blade.materials(k).etal;
    matdb(k).reference = blade.materials(k).reference;
end
N = numel(matdb);
for k=1:numel(blade.stacks)
    matdb(N+k).name          = blade.stacks(k).name;
    matdb(N+k).type          = 'composite';
    matdb(N+k).reference     = 'Reference text';
    matdb(N+k).thicknessType = 'Constant';
    matdb(N+k).uniqueLayers  = numel(blade.stacks(k).plygroups);
    matdb(N+k).symmetryType  = 'none';
    for j=1:matdb(N+k).uniqueLayers
        matid = blade.stacks(k).plygroups(j).materialid;
        matdb(N+k).layer(j).layerName  = matdb(matid).name;
        matdb(N+k).layer(j).thicknessA = mm_to_m * blade.stacks(k).plygroups(j).thickness;
        matdb(N+k).layer(j).thicknessB = matdb(N+k).layer(j).thicknessA;
        matdb(N+k).layer(j).quantity   = blade.stacks(k).plygroups(j).nPlies;
        matdb(N+k).layer(j).theta      = blade.stacks(k).plygroups(j).angle;
    end
end
for kw=1:numel(blade.swstacks)
    N = numel(matdb);
    for k=1:numel(blade.swstacks{kw})
        matdb(N+k).name          = blade.swstacks{kw}(k).name;
        matdb(N+k).type          = 'composite';
        matdb(N+k).reference     = 'Reference text';
        matdb(N+k).thicknessType = 'Constant';
        matdb(N+k).uniqueLayers  = numel(blade.swstacks{kw}(k).plygroups);
        matdb(N+k).symmetryType  = 'none';
        for j=1:matdb(N+k).uniqueLayers
            matid = blade.swstacks{kw}(k).plygroups(j).materialid;
            matdb(N+k).layer(j).layerName  = matdb(matid).name;
            matdb(N+k).layer(j).thicknessA = mm_to_m * blade.swstacks{kw}(k).plygroups(j).thickness;
            matdb(N+k).layer(j).thicknessB = matdb(N+k).layer(j).thicknessA;
            matdb(N+k).layer(j).quantity   = blade.swstacks{kw}(k).plygroups(j).nPlies;
            matdb(N+k).layer(j).theta      = blade.swstacks{kw}(k).plygroups(j).angle;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Write the NuMAD files          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     assignin('base','data',data);
%     assignin('base','matdb',matdb);
writeNuMADinput(data,numad_file);
writeMatDB(matdb,matdb_file);
end