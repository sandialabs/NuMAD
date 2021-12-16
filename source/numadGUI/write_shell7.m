function varargout = write_shell7(app,blade,filename)
%WRITE_SHELL7  Generate the ANSYS input file that creates the blade 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   write_shell7(app,blade,filename)
%     app = app data structure in NuMAD
%     blade = blade data structure in NuMAD
%     filename = ANSYS input file to be written (typically 'shell7.src')
%

% 2011.06.24  jcb: Fixed a logical error in which DPs could be skipped
%    while forming the KeyPoints structure.

%jcb: passing both "app" and "blade" is temporary until the data
%     re-organization is finalized

% original TCL code comments are prefaced with "tcl:"

%tcl: DO NOT USE DOUBLE QUOTES IN THE shell7 CONTENTS
   
%tcl:         Created by Daniel Laird for Sandia Labs
%tcl:         3-D shell analysis
%tcl:         Modified on April 26, 2000
%tcl:
%tcl:         THIS FILE IS FOR USE WITH INTERACTIVE SHELL MODEL INPUT
%tcl:         IT USES SPLINE-ING CAPABILITIES WITHIN ANSYS
   
%tcl:         ANSYS source file for use with NuMAD

%tcl: Copy macros into the working directory
%jcb:  need to decide when (in what script) to copy macros over
if isequal(0,nargout)
    parent_pn = app.numadpath;
    [success,message,~] = copyfile(fullfile(parent_pn,'macros'),app.settings.job_path);
    if ~success
        errordlg(message,'write_shell7: error copying macros');
        return;
    end
end

% % Attempt to get version number of ANSYS
% ansys_path = app.settings.ansys_path;
% t = regexp(ansys_path,'v(\d+)|ANSYS(\d+)','tokens');
% if ~isempty(t)
%     ansys_version = t{1}{1};  % expect to get '140' for Version 14.0
% else
%     ansys_version = 'unknown';
% end

TotalStations = numel(app.station);
TotalShearwebs = numel(app.shearweb);

% Insert the DPs into the airfoil coordinates and remember their indices
InterpMethod = 'spline';  % interpolation method for DPs
KeyPoints = struct('x',[],'y',[],'LE',[],'DP',[]);
for kStation = 1:TotalStations
    station = app.station(kStation);
    nCoordPairs = size(station.coords,1);
    %ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    % ANSYS will error out when station length is greater than ~175
    Nansys = 150;
% %     if nCoordPairs > Nansys
        % resample airfoil coordinates to reduce size        
        af_out = resampleAirfoil(station.coords, round(Nansys/2-1), 'arc');
                
% %         figure
% %         hold on
% %         plot(station.coords(:,1),station.coords(:,2),'k--')
% %         plot(af_out(:,1),af_out(:,2),'^')
% %         grid on; axis equal
        
        % redefine station coordinates for lower point count
        station.coords = af_out;
        app.station(kStation).coords = af_out; % will this affect other uses? Does app get saved back to NuMAD?? Is this necessary??
        nCoordPairs = size(station.coords,1);
% %     end
    %ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    nDPs = numel(station.dp) - 2;
    [KeyPoints(kStation).x, KeyPoints(kStation).y] = deal(zeros(nCoordPairs+nDPs,1));
    station.xn = station.coords(:,1);
    [~,j] = min(station.xn);
    station.coordsLE = j;
    KeyPoints(kStation).LE = station.coordsLE;
    if station.xn(j) >= 0  % jcb: what was my logic here?
        j = j-1;
    end
    station.xn(1:j) = -1*station.xn(1:j);  % give the HP points a negative xn value
    % jcb: The following loop inserts the DPs into the coordinate pairs.  
    %    If the DP matches a CoordPair, a point is not inserted and we must 
    % reduce the length of KeyPoints.x and .y accordingly (the 'cp_delete'
    % variable keeps track of how many DPs matched CoordPairs.
    %    A DP is inserted when the xn (signed chord value) of the current
    % CoordPair is greater than the current DP's signed chord position.
    %    After insertion of a new point, the offset between indices of 
    % station.coords and KeyPoints.xy must increase by one (the 'cp_offset' 
    % variable keeps track of how many DPs have been inserted).
    kcp = 1; kdp = 1;
    cp_offset = 0;
    cp_delete = 0;
    
    while (kcp <= nCoordPairs)
        if (abs(station.xn(kcp) - station.dp(kdp+1)) <= 1e-4) && (kdp <= nDPs)
            % the current coord pair matches the current DP (within the tolerance)  
            % jcb: should the user be able to set the tolerance?
            KeyPoints(kStation).DP(kdp) = kcp+cp_offset;
            kdp = kdp + 1;  % next DP;
            cp_delete = cp_delete + 1;  % KeyPoints length will be reduced by 1
            KeyPoints(kStation).x(kcp+cp_offset) = station.coords(kcp,1);
            KeyPoints(kStation).y(kcp+cp_offset) = station.coords(kcp,2);
            kcp = kcp + 1;  % next CoordPair
            % jcb: ToDo - in model check, need to ensure two DPs are not at
            % the same point (note: equivalent to Hourglass)
        elseif (station.xn(kcp) > station.dp(kdp+1)) && (kdp <= nDPs)
            KeyPoints(kStation).x(kcp+cp_offset) = abs(station.dp(kdp+1));
            arc = interp1(station.c(2:end-1),station.s(2:end-1),station.dp(kdp+1),'linear');
            KeyPoints(kStation).y(kcp+cp_offset) = interp1(station.s(2:end-1),station.xy(2:end-1,2),arc,InterpMethod);
            KeyPoints(kStation).DP(kdp) = kcp+cp_offset;
            kdp = kdp + 1;  % next DP
            cp_offset = cp_offset + 1;  % offset next coord pair by one more
            % jcb: 2012-08-28 comparison below needs to be '<=' because DP
            % can come immediately before LE coordinate point (but the DP 
            % isn't exaclty at the LE because the first IF condition above 
            % failed)
            if kcp <= station.coordsLE  
                % LE index has moved if we inserted a HP dp
                KeyPoints(kStation).LE = KeyPoints(kStation).LE + 1;
            end
        else
            KeyPoints(kStation).x(kcp+cp_offset) = station.coords(kcp,1);
            KeyPoints(kStation).y(kcp+cp_offset) = station.coords(kcp,2);
            kcp = kcp + 1;  % next CoordPair
        end
    end
    % need to add any DPs that fall after the last CoordPair
    while kdp <= nDPs
        KeyPoints(kStation).x(kcp+cp_offset) = abs(station.dp(kdp+1));
        arc = interp1(station.c(2:end-1),station.s(2:end-1),station.dp(kdp+1),'linear');
        KeyPoints(kStation).y(kcp+cp_offset) = interp1(station.s(2:end-1),station.xy(2:end-1,2),arc,InterpMethod);
        KeyPoints(kStation).DP(kdp) = kcp+cp_offset;
        kdp = kdp + 1;  % next DP
        cp_offset = cp_offset + 1;  % offset next coord pair by one more
    end
    if cp_delete > 0
        KeyPoints(kStation).x(end-cp_delete+1:end) = [];
        KeyPoints(kStation).y(end-cp_delete+1:end) = [];
    end
    switch station.TEtype
        case {'round','sharp'}
            KeyPoints(kStation).DP = [1 KeyPoints(kStation).DP 1];
        case {'flat'}
            KeyPoints(kStation).DP = [1 2 KeyPoints(kStation).DP(2:end) 1];
    end
end

SkinAreas(1:TotalStations-1) = struct('startIB',[],'endIB',[],'startOB',[],'endOB',[],'Material',{''});
for kStation = 1:numel(SkinAreas)
    stationIB = app.station(kStation);
    stationOB = app.station(kStation+1);
    
    for kdp = 1:numel(stationIB.dp)-1
        switch stationIB.dptype{kdp}
            case {'single','double'}
                % single and double are equivalent on the area inboard edge
                % start and end are current and next DP
                SkinAreas(kStation).startIB(end+1) = kdp;
                SkinAreas(kStation).endIB(end+1)   = kdp+1;
                SkinAreas(kStation).Material{end+1} = stationIB.sm{kdp};
            case {'flare','hourglass'}
                % flare and hourglass are equivalent on the area inboard edge
                % start and end of first area is current DP
                % start and end of next area is current and next DP
                SkinAreas(kStation).startIB(end+(1:2)) = [kdp kdp];
                SkinAreas(kStation).endIB(end+(1:2))   = [kdp kdp+1];
                SkinAreas(kStation).Material{end+1} = stationIB.dpmaterial{kdp};  %jcb:  append (end+1) dp material
                SkinAreas(kStation).Material{end+1} = stationIB.sm{kdp};  %jcb:  append (end+1) segment material
        end     
    end
    
    for kdp = 1:numel(stationOB.dp)-1
        switch stationOB.dptype{kdp}
            case {'single','flare'}
                % single and flare are equivalent on the area outboard edge
                % start and end are current and next DP
                SkinAreas(kStation).startOB(end+1) = kdp;
                SkinAreas(kStation).endOB(end+1)   = kdp+1;
            case {'double','hourglass'}
                % double and hourglass are equivalent on the area outboard edge
                % start and end of first area is current DP
                % start and end of next area is current and next DP
                SkinAreas(kStation).startOB(end+(1:2)) = [kdp kdp];
                SkinAreas(kStation).endOB(end+(1:2))   = [kdp kdp+1];
        end     
    end
end

%tcl: Determine which composite materials are used in the model
compsInModel = {};
%tcl: search shear web materials
for k = 1:TotalShearwebs
    compsInModel = [compsInModel; {app.shearweb(k).Material}];  %#ok
end
%tcl: search skin materials
for k = 1:numel(SkinAreas)
    compsInModel = [compsInModel; transpose(SkinAreas(k).Material)];  %#ok
end
compsInModel = unique(compsInModel);

% load the material database and create searchable list
if ~isempty(app.settings.job_name)
    % job name exists, use the local material database
    app.matdb_path = fullfile(app.settings.job_path,'MatDBsi.txt');
else
    % if no job name, use the master material database
    % jcb: we shouldn't arrive here because a file save is required first
    app.matdb_path = fullfile(app.numadpath,'MatDBsi.txt');
end
app.matdb = readMatDB(app.matdb_path);
for k=1:numel(app.matdb)
    app.matlist{k} = app.matdb(k).name;
    mattype{k} = app.matdb(k).type;  %#ok
end
app.isotropic =  strcmp('isotropic',mattype);
app.orthotropic = strcmp('orthotropic',mattype);
app.composite = strcmp('composite',mattype);

% Determine which isotropic and orthotropic materials are used in the model
isoorthoInModel = {};
for kcomp = 1:numel(compsInModel)
    n = strcmp(compsInModel(kcomp),app.matlist);
    if ~any(n)
        errordlg(sprintf('Material "%s" not found in database.',compsInModel{kcomp}),'Error');
        error('Material "%s" not found in database.',compsInModel{kcomp});
    end
    mat = app.matdb(n);
    layerNames = cell(numel(mat.layer),1);
    for klay = 1:numel(mat.layer)
        layerNames(klay) = {mat.layer(klay).layerName};
    end
    isoorthoInModel = unique([isoorthoInModel; layerNames]);
end

if isequal(3,nargout)
    varargout{1} = KeyPoints;
    varargout{2} = SkinAreas;
    varargout{3} = isoorthoInModel;
    return   % stop here and do not write output
end

% % DEBUGGING
% assignin('base','KeyPoints',KeyPoints);
% assignin('base','SkinAreas',SkinAreas);

% % DEBUGGING
% assignin('base','app_shell7',app);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write the shell7.src file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fid = fopen('shell7.src','wt');
fid = fopen(filename,'wt');
try
    fprintf(fid,'hide_warndlg_keyopt=1\n');
    fprintf(fid,'hide_warndlg_areas=1\n');
    
    fprintf(fid,'\n/nerr,500,50000\n');
    fprintf(fid,'\n/filename,master\n');
    
    %tcl: define location of ANSYS macro library
    fprintf(fid,'\n/prep7\n');
    
    %tcl: DEFINE ELEMENT TYPES
    fprintf(fid,'\n! DEFINE ELEMENT TYPES =================================\n');
    
    %tcl: shell99, constant thickness layers, store data for all layers, nodes at mid surface
    fprintf(fid,'\n   et,11,shell99,,0');
    fprintf(fid,'\n   keyopt,11,8,1\n');
    
    %tcl: shell99, constant thickness layers, store data for all layers, nodes at bottom surface
    fprintf(fid,'\n   et,12,shell99,,0');
    fprintf(fid,'\n   keyopt,12,8,1');
    fprintf(fid,'\n   keyopt,12,11,1\n');
    
    %tcl: shell99, constant thickness layers, store data for all layers, nodes at top surface
    fprintf(fid,'\n   et,14,shell99,,0');
    fprintf(fid,'\n   keyopt,14,8,1');
    fprintf(fid,'\n   keyopt,14,11,2\n');

    %tcl: 50 layer shell, tapered layers, store data for all layers, nodes at top surface
    fprintf(fid,'\n   et,13,shell99,,1');
    fprintf(fid,'\n   keyopt,13,8,1');
    fprintf(fid,'\n   keyopt,13,11,1\n');

    %tcl: shell99, tapered layers, store data for all layers, nodes at mid surface
    fprintf(fid,'\n   et,15,shell99,,1');
    fprintf(fid,'\n   keyopt,15,8,1\n');

    %tcl: shell91, store data for all layers, nodes at mid surface
    fprintf(fid,'\n   et,16,shell91,100');
    fprintf(fid,'\n   keyopt,16,8,1\n');

    %tcl: shell91, store data for all layers, nodes at bottom surface
    fprintf(fid,'\n   et,17,shell91,100');
    fprintf(fid,'\n   keyopt,17,8,1');
    fprintf(fid,'\n   keyopt,17,11,1\n');

    %tcl: shell91, store data for all layers, nodes at top surface
    fprintf(fid,'\n   et,18,shell91,100');
    fprintf(fid,'\n   keyopt,18,8,1');
    fprintf(fid,'\n   keyopt,18,11,2\n');

    %tcl: shell91, sandwich option, store data for all layers, nodes at mid surface
    fprintf(fid,'\n   et,19,shell91,100');
    fprintf(fid,'\n   keyopt,19,5,1');
    fprintf(fid,'\n   keyopt,19,8,1');
    fprintf(fid,'\n   keyopt,19,9,1\n');

    %tcl: shell91, sandwich option store data for all layers, nodes at bottom surface
    fprintf(fid,'\n   et,20,shell91,100');
    fprintf(fid,'\n   keyopt,20,5,1');
    fprintf(fid,'\n   keyopt,20,8,1');
    fprintf(fid,'\n   keyopt,20,11,1');
    fprintf(fid,'\n   keyopt,20,9,1\n');

    %tcl: shell91, sandwich option store data for all layers, nodes at top surface
    fprintf(fid,'\n   et,22,shell91,100');
    fprintf(fid,'\n   keyopt,22,5,1');
    fprintf(fid,'\n   keyopt,22,8,1');
    fprintf(fid,'\n   keyopt,22,11,2');
    fprintf(fid,'\n   keyopt,22,9,1\n');

    %tcl: structural mass
    fprintf(fid,'\n   et,21,mass21,,,0');
    fprintf(fid,'\n   r,999,0.0,0.0,0.00001,0.0,0.0,0.0\n');

    %tcl: shell281, 8-node structural shell, store data for TOP, BOTTOM,
    %     and MID for all layers
    fprintf(fid,'\n   et,31,shell281');
    fprintf(fid,'\n   keyopt,31,8,2');
    fprintf(fid,'\n*if,hide_warndlg_keyopt,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 3'); % turn warning dialog off
    fprintf(fid,'\n*endif');
    fprintf(fid,'\n   !Set keyopt(2)=1 for improved formulation in R12 & R13');
    fprintf(fid,'\n   keyopt,31,2,1');
    fprintf(fid,'\n*if,hide_warndlg_keyopt,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 2'); % turn warning dialog on
    fprintf(fid,'\n*endif\n');
    %jcb: I thought about checking the ansys version and making conditional
    %statements, but someone could share the shell7 with someone using a
    %different version of ansys
%     if strncmp(ansys_version,'12',2) || strncmp(ansys_version,'13',2)
%         % Releases 12 & 13 of ANSYS require keyopt(2)=1 with shell281
%         %  to make use of the improved formulation
%         fprintf(fid,'\n   keyopt,31,2,1\n');
%     end
  
    % shell181, 4-node structural shell, store data for TOP, BOTTOM, and 
    %    MID for all layers
    fprintf(fid,'\n   et,32,shell181');
    fprintf(fid,'\n   keyopt,32,8,2\n');
    fprintf(fid,'\n   keyopt,32,3,2\n');  % added by brr 10/25/2012

    %tcl: Write material properties
    %tcl:    This changed dramatically on 2001 November 08
    %tcl:    Now only materials used in the model are written to shell7.src
    %tcl:    Also, material numbers are no longer recorded until write_shell7
    %tcl:    Two new local arrays ansysMPnumber and ansysRnumber are used within write_shell7
    
    fcfields = {'xten','xcmp','yten','ycmp','zten','zcmp','xy','yz','xz',...
        'xycp','yzcp','xzcp','xzit','xzic','yzit','yzic','g1g2',...
        'etal','etat','alp0'};
    fcvalues = cell(numel(fcfields),1);
    
    fprintf(fid,'\n! WRITE MATERIAL PROPERTIES ============================\n');
    fprintf(fid,'\n  ! FAILURE CRITERIA LIMIT TABLE LEGEND:');
    fprintf(fid,'\n  ! tb,fcli,<mat>,ntemp,npts,tbopt');
    fprintf(fid,'\n  !   (tbopt=1 for stress limits; default ntemp=1,npts=20)');
    fprintf(fid,'\n  ! tbdata,1,xten,xcmp,yten,ycmp,zten,zcmp');
    fprintf(fid,'\n  ! tbdata,7,xy,yz,xz,xycp,yzcp,xzcp');
    fprintf(fid,'\n  ! tbdata,13,xzit,xzic,yzit,yzic');
    fprintf(fid,'\n  ! tbdata,17,g1g2,etal,etat,alp0\n');
    for kmp = 1:numel(isoorthoInModel)
        n = strcmp(isoorthoInModel(kmp),app.matlist);
        mat = app.matdb(n);
        switch mat.type
            case 'isotropic'
                fprintf(fid,'\n   ! %s'           ,mat.name);
                fprintf(fid,'\n   mp,ex,%d,%g'    ,kmp,mat.ex);
                fprintf(fid,'\n   mp,dens,%d,%g'  ,kmp,mat.dens);
                fprintf(fid,'\n   mp,nuxy,%d,%g',kmp,mat.nuxy);
            case 'orthotropic'
                fprintf(fid,'\n   ! %s'         ,mat.name);
                fprintf(fid,'\n   mp,ex,%d,%g'  ,kmp,mat.ex);
                fprintf(fid,'\n   mp,ey,%d,%g'  ,kmp,mat.ey);
                fprintf(fid,'\n   mp,ez,%d,%g'  ,kmp,mat.ez);
                fprintf(fid,'\n   mp,prxy,%d,%g',kmp,mat.prxy);
                fprintf(fid,'\n   mp,pryz,%d,%g',kmp,mat.pryz);
                fprintf(fid,'\n   mp,prxz,%d,%g',kmp,mat.prxz);
                fprintf(fid,'\n   mp,gxy,%d,%g' ,kmp,mat.gxy);
                fprintf(fid,'\n   mp,gyz,%d,%g' ,kmp,mat.gyz);
                fprintf(fid,'\n   mp,gxz,%d,%g' ,kmp,mat.gxz);
                fprintf(fid,'\n   mp,dens,%d,%g',kmp,mat.dens);
            otherwise
                error('Unknown material type in database');
        end
        switch mat.type
            case {'isotropic','orthotropic'}
                % Note that entering a blank or a zero for XYCP,YZCP, or XZCP
                % triggers the default value of -1.0. To specify an effective zero,
                % use a small, nonzero value (such as 1E-14).
                if isequal(0,mat.xycp), mat.xycp=1e-14; end
                if isequal(0,mat.yzcp), mat.yzcp=1e-14; end
                if isequal(0,mat.xzcp), mat.xzcp=1e-14; end
                % convert degrees to radians
                %mat.alp0 = mat.alp0 * pi/180;
                % read all of the failure criteria values
                for kfc = 1:numel(fcfields)
                    fcname = fcfields{kfc};
                    fcvalues{kfc} = mat.(fcname);
                end
                if all(cellfun('isempty',fcvalues))
                    % do not print anything if all failure criteria
                    % properties are empty
                else
                    fprintf(fid,'\n   tb,fcli,%d,1,20,1',kmp);
                    fprintf(fid,'\n   tbdata,1,%g,%g,%g,%g,%g,%g',mat.xten,mat.xcmp,mat.yten,mat.ycmp,mat.zten,mat.zcmp);
                    fprintf(fid,'\n   tbdata,7,%g,%g,%g,%g,%g,%g',mat.xy,mat.yz,mat.xz,mat.xycp,mat.yzcp,mat.xzcp);
                    fprintf(fid,'\n   tbdata,13,%g,%g,%g,%g',mat.xzit,mat.xzic,mat.yzit,mat.yzic);
                    fprintf(fid,'\n   tbdata,17,%g,%g,%g,%g',mat.g1g2,mat.etal,mat.etat,mat.alp0);
%                     for kf = 1:numel(fcvalues)
%                         if ~isempty(fcvalues{kf})
%                             fprintf(fid,'\n   tbdata,%d,%g',kf,fcvalues{kf});
%                         end
%                     end
                end
            otherwise
                error('Unknown material type in database');
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'\n! WRITE THE COMPOSITES =================================\n');
    rCounter=1;
    switch app.ansys.ElementSystem
        %tcl:  first few lines are same for shell91 and shell99
        case {'91','99'}
            for rCounter = 1:numel(compsInModel)
                n = strcmp(compsInModel(rCounter),app.matlist);
                mat = app.matdb(n);
                for klay = 1:numel(mat.layer)
                    % shell91/99: layer.quantity only multiplies layer
                    % thickness and does not repeat layers
                    mat.layer(klay).thicknessA = mat.layer(klay).thicknessA * mat.layer(klay).quantity;
                    mat.layer(klay).thicknessB = mat.layer(klay).thicknessB * mat.layer(klay).quantity;
                end
                fprintf(fid,'\n   ! %s',mat.name);
                LSYM = 0;
                if ~isequal(mat.symmetryType,'none')
                    LSYM = 1;
                end
                fprintf(fid,'\n   r,%d,%d,%d',rCounter,mat.uniqueLayers,LSYM);
                fprintf(fid,'\n   rmore');
                
                if isequal(app.ansys.ElementSystem,'91')
                    for klay = 1:numel(mat.layer)
                        ansysMPnumber = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                        fprintf(fid,'\n   rmore,%d,%g,%g,%g,%g,%g',...
                            ansysMPnumber,mat.layer(klay).theta,...
                            mat.layer(klay).thicknessA,mat.layer(klay).thicknessB,...
                            mat.layer(klay).thicknessB,mat.layer(klay).thicknessA);
                    end
                elseif isequal(app.ansys.ElementSystem,'99')
                    if isequal(mat.thicknessType,'Constant')
                        for klay = 1:2:floor(numel(mat.layer)/2)*2
                            ansysMPnumber(1) = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                            ansysMPnumber(2) = find(strcmp(mat.layer(klay+1).layerName,isoorthoInModel)==1);
                            fprintf(fid,'\n   rmore,%d,%g,%g,%d,%g,%g',...
                                ansysMPnumber(1),mat.layer(klay).theta,...
                                mat.layer(klay).thicknessA,...
                                ansysMPnumber(2),mat.layer(klay+1).theta,...
                                mat.layer(klay+1).thicknessA); 
                        end
                        %tcl: ### WRITE LAST LAYER IF TOTAL NUMBER OF LAYERS IS ODD
                        if klay == numel(mat.layer)
                            ansysMPnumber = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                            fprintf(fid,'\n   rmore,%d,%g,%g',...
                                ansysMPnumber,mat.layer(klay).theta,...
                                mat.layer(klay).thicknessA);
                        end
                        
                    elseif isequal(mat.thicknessType,'Tapered')
                        for klay = 1:numel(mat.layer)
                            ansysMPnumber = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                            fprintf(fid,'\n   rmore,%d,%g,%g,%g,%g,%g',...
                                ansysMPnumber,mat.layer(klay).theta,...
                                mat.layer(klay).thicknessA,mat.layer(klay).thicknessB,...
                                mat.layer(klay).thicknessB,mat.layer(klay).thicknessA);
                        end 
                    end
                end
                fprintf(fid,'\n');
            end
            
        %case '191'
            
        case {'281','181'}
            for secCounter = 1:numel(compsInModel)
                n = strcmp(compsInModel(secCounter),app.matlist);
                mat = app.matdb(n);
                if isequal('multiply',app.ansys.MultipleLayerBehavior)
                    for klay = 1:numel(mat.layer)
                        mat.layer(klay).thicknessA = mat.layer(klay).thicknessA * mat.layer(klay).quantity;
                        mat.layer(klay).thicknessB = mat.layer(klay).thicknessB * mat.layer(klay).quantity;
                        mat.layer(klay).quantity = 1;
                    end
                end
                fprintf(fid,'\n   ! %s',mat.name);
                %tcl: define section type
                fprintf(fid,'\n   sectype,%d,shell'  ,secCounter);
                for klay = 1:numel(mat.layer)
                    % jcb: might want to change the following line to use
                    % key/value approach rather than searching every time
                    ansysMPnumber = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                    %tcl: define layers corresponding to section type
                    for kqty = 1:mat.layer(klay).quantity
                    fprintf(fid,'\n      secdata,%g,%d,%g,,%s',...
                        mat.layer(klay).thicknessA,...
                        ansysMPnumber,...
                        mat.layer(klay).theta,...
                        mat.layer(klay).layerName);
                    end
                end
                %tcl: define section offset
                fprintf(fid,'\n   secoffset,bot\n');
                
                %tcl: duplicate section type information for possible use in shear web
                %tcl:   This is necessary because the nodal location is part of the material section definition
                fprintf(fid,'\n   ! %s',mat.name);
                fprintf(fid,'\n   sectype,%d,shell'  ,1000+secCounter);
                for klay = 1:numel(mat.layer)
                    % jcb: might want to change the following line to use
                    % key/value approach rather than searching every time
                    ansysMPnumber = find(strcmp(mat.layer(klay).layerName,isoorthoInModel)==1);
                    %tcl: define layers corresponding to section type
                    for kqty = 1:mat.layer(klay).quantity
                    fprintf(fid,'\n      secdata,%g,%d,%g,,%s',...
                        mat.layer(klay).thicknessA,...
                        ansysMPnumber,...
                        mat.layer(klay).theta,...
                        mat.layer(klay).layerName);
                    end
                end
                %tcl: define section offset
                fprintf(fid,'\n   secoffset,mid\n');
            end
        otherwise
            errordlg(sprintf('Element System %s not yet available',app.ansys.ElementSystem),'write_shell7 error')
            error('Element System %s not yet available',app.ansys.ElementSystem);
            
    end
    
    [~,jobtitle,~] = fileparts(app.settings.job_name);
    fprintf(fid,'\n/title,%s',jobtitle);
    fprintf(fid,'\nZrCount=%d\n',rCounter);
    
%jcb:  it appears these "Mass*" variables are not used anywhere
% tcl % tcl % tcl % tcl % tcl % tcl % tcl % tcl % tcl 
%    set numZelem 100
%    #        for {set i 2} {$i <= $TotalStations} {incr i} {
%    #           set numZelem [expr $numZelem + $SelemDiv($i)]
%    #        }
%    
%    
%    if {$MassDistributionCheckStatus == 1} {
%       append shell7Contents "
%          *dim,MassVal,,$numZelem
%          *dim,MassRmin,,$numZelem
%          *dim,MassRmax,,$numZelem
%       "
%    }
% tcl % tcl % tcl % tcl % tcl % tcl % tcl % tcl % tcl 
    

    %tcl: DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES
    %tcl:    THE LINES ARE PRODUCED WITH THREE DIFFERENT SPLINING MACROS    
    fprintf(fid,'\n! DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES\n');
    
    % Create a coordinate system roughly in the fiber direction (+X down blade, +Z up toward LP side)
    % --> beginning with global csys, rotate -90 about local Z, then -90 about local Y
    fprintf(fid,'\nlocal,1000,CART,0,0,0, -90,0,-90\n');
    twistFlag = 1;  % ccw rotor rotation
    if strcmp(app.BladeRotation,'cw')
        twistFlag = -1;  % cw rotor rotation
    end
    for kStation = 1:TotalStations
        station = app.station(kStation);
        keypts = KeyPoints(kStation);
        x = (keypts.x - station.Xoffset) * station.Chord * twistFlag;
        y = (keypts.y              ) * station.Chord;
        %jcb: as of 2011-05-26, the keypoints are transformed directly
        %     rather than with ansys commands
        twist = twistFlag * station.DegreesTwist * pi/180;
        coords = zeros(length(x),4);
        coords(:,1) = cos(twist) * x - sin(twist) * y;
        coords(:,2) = sin(twist) * x + cos(twist) * y;
        %coords(:,3) = zeros(size(x));
        coords(:,4) = ones(size(x));
%         xt = cos(twist) * x - sin(twist) * y;
%         yt = sin(twist) * x + cos(twist) * y;
%         zt = ones(size(x)) * station.LocationZ;
%         xn = station.xn;

        % use the generating line to translate and rotate the coordinates
        presweep_slope = ppval(blade.PresweepRef.dpp,station.LocationZ);
        precurve_slope = ppval(blade.PrecurveRef.dpp,station.LocationZ);
        presweepDeg = 180/pi*atan(presweep_slope*twistFlag);
        precurveDeg = 180/pi*atan(-precurve_slope);
        [presweep_rot, precurve_rot] = deal(0);
        if isequal(blade.PresweepRef.method,'normal')
            presweep_rot = atan(presweep_slope*twistFlag);
        end
        if isequal(blade.PrecurveRef.method,'normal')
            precurve_rot = atan(-precurve_slope);
        end
        transX = twistFlag*ppval(blade.PresweepRef.pp,station.LocationZ);
        transY = ppval(blade.PrecurveRef.pp,station.LocationZ);
        R = makehgtform('yrotate',presweep_rot,'xrotate',precurve_rot);
        T = makehgtform('translate',transX,transY,station.LocationZ);
        coords = coords * R' * T';
        xt = coords(:,1);
        yt = coords(:,2);
        zt = coords(:,3);

        % ensure we are in csys0 and no keypoints are selected
        fprintf(fid,'\ncsys,0');
        fprintf(fid,'\n   ksel,none');
        for kcp = 1:numel(x)
            %jcb:  old method - unrotated points
            %fprintf(fid,'\n      k,%d,%g,%g,%g',((kStation-1)*1000+kcp),x(kcp),y(kcp),station.LocationZ);
            %jcb:  new method - points transformed by twist and generating line
            fprintf(fid,'\n      k,%d,%g,%g,%g',((kStation-1)*1000+kcp),xt(kcp),yt(kcp),zt(kcp));
        end
        switch station.TEtype
            case 'round'
                fprintf(fid,'\n   zSmoothe');
            case 'sharp'
                fprintf(fid,'\n   zAirfoil');
            case 'flat'
                fprintf(fid,'\n   zFlatback');
        end
        % Create a coordinate system to be used later for aligning the fiber direction.
        % First, load the csys defined earlier (+X down blade, +Z up toward LP side)
        fprintf(fid,'\n   csys,1000');
        % Next, translate & rotate relative to this active csys (use CLOCAL, not LOCAL)
        % translation: global X,Y,Z => local y,z,x
        % rotation: presweep is local z rotation & precurve is local y rotation
        fprintf(fid,'\n   clocal,%d,CART,%g,%g,%g, %g,%g,%g\n',(1000+kStation),...
            station.LocationZ,transX,transY, presweepDeg,0,precurveDeg);
        
        % Create coordinate system at the tip
        fprintf(fid,'\nlocal,12,CART,%g,%g,%g, %g,%g,%g\n',...
            transX,transY,station.LocationZ, 0,precurve_rot*180/pi,presweep_rot*180/pi);
    end
    fprintf(fid,'\n   csys,0');
    fprintf(fid,'\nksel,all\n');
    
    %jcb: as of 2011-05-26, the keypoints are transformed directly
    %     rather than with ansys commands
%     fprintf(fid,'\n! ROTATE SECTIONS ======================================\n');
%     fprintf(fid,'\ncsys,1');
%     for kStation = 1:TotalStations
%         fprintf(fid,'\n   lsel,s,loc,z,%g',app.station(kStation).LocationZ);
%         fprintf(fid,'\n   lgen,2,all,,,,%g,,,,1\n',twistFlag*app.station(kStation).DegreesTwist);
%     end
%     fprintf(fid,'\ncsys,0');
    fprintf(fid,'\nallsel\n');
    
    %tcl: GENERATE SPANWISE LINES FOR LATER AREA CREATION
    fprintf(fid,'\n! GENERATE SPANWISE LINES FOR LATER AREA CREATION ======\n');
    %tcl:    THIS IS TO AVOID STRANGE BLADE ENVELOPES RESULTING FROM VERY CHORDWISE TAPERED SKIN MATERIALS
    %tcl:    added Sept 2008
   
    %tcl:    HP and LP areas modified January 2010 to no longer wrap around the knee of a flatback airfoil
   
    %tcl:    determine maximum line number before starting this process
    fprintf(fid,'\n*get,z_lines_after_splining,line,,num,max\n');   %jcb:  is "z_lines_after_splining" used anywhere?
    
    %tcl:    dimension arrays
    fprintf(fid,'\n*dim,z_HP_line,,%d',TotalStations);
    fprintf(fid,'\n*dim,z_LP_line,,%d',TotalStations);
    fprintf(fid,'\n*dim,z_HP_area,,%d',TotalStations-1);
    fprintf(fid,'\n*dim,z_LP_area,,%d',TotalStations-1);
    fprintf(fid,'\n*dim,z_TE_line,,%d',TotalStations-1);
    fprintf(fid,'\n*dim,z_LE_line,,%d',TotalStations-1);
    fprintf(fid,'\n*dim,z_LP_knee_line,,%d',TotalStations-1);
    fprintf(fid,'\n*dim,z_HP_knee_line,,%d\n',TotalStations-1);
    
    
    %tcl:    Generate HP and LP lines at each station
    fprintf(fid,'\n! Generate HP and LP lines at each station =============\n');
    for kStation = 1:TotalStations
        isFlatback = strcmp(app.station(kStation).TEtype,'flat');
        nCoordPairs = numel(KeyPoints(kStation).x);
        hpTE = (1000*(kStation-1)+1);
        LE   = (1000*(kStation-1)+KeyPoints(kStation).LE);
        lpTE = (1000*(kStation-1)+nCoordPairs);
        
        %tcl:  select KPs for HP
        fprintf(fid,'\n   ksel,s,kp,,%d,%d',hpTE+isFlatback,LE);
        %tcl:  select lines completely within KP set
        fprintf(fid,'\n      lslk,s,1');
        %tcl:  combine lines
        fprintf(fid,'\n      lcomb,all,,1');
        %tcl:  determine the new line number for HP
        fprintf(fid,'\n      *get,z_HP_line(%d),line,,num,max',kStation);
        
        %tcl:  select most KPs for LP (can't do them all at once or the HP line is inadvertantly selected)
        fprintf(fid,'\n   ksel,s,kp,,%d,%d',LE,lpTE);
        %tcl:  select lines completely within KP set
        fprintf(fid,'\n      lslk,s,1');
        %tcl: if flatback, all necessary lines have been selected
        %tcl: if not, final lines needs to be added to set
        if isFlatback == 0
            %tcl:  select first and last KPs for LP (this will allow last line to be selected without selecting HP)
            fprintf(fid,'\n         ksel,s,kp,,%d,%d,%d',hpTE,lpTE,nCoordPairs-1);
            %tcl:  add lines completely within KP set (final line in station)
            fprintf(fid,'\n         lslk,a,1');
        end
        %tcl:  combine lines
        fprintf(fid,'\n      lcomb,all,,1');
        %tcl:  determine the new line number for LP
        fprintf(fid,'\n      *get,z_LP_line(%d),line,,num,max',kStation);
        fprintf(fid,'\n');
    end
    
    
    %tcl:    Generate TE and LE spanwise lines from station to station
    %tcl:    Also generate spanwise lines for flatback knees
    fprintf(fid,'\n! Generate TE and LE spanwise lines from station to station\n');
    fprintf(fid,'\ncsys,0');
    fprintf(fid,'\nallsel');
    for kStation = 1:TotalStations-1
        %tcl:  create TE line
        fprintf(fid,'\n   l,%d,%d',(1000*(kStation-1)+1),(1000*kStation+1));
        %tcl:  determine the new line number for TE
        fprintf(fid,'\n   *get,z_TE_line(%d),line,,num,max',kStation);
        %tcl:  create flatback spanwise knee lines if neccessary
        thisIsFlatback = strcmp(app.station(kStation).TEtype,'flat');
        nextIsFlatback = strcmp(app.station(kStation+1).TEtype,'flat');
        if thisIsFlatback || nextIsFlatback
            %tcl: create HP_knee line
            fprintf(fid,'\n   l,%d,%d',...
                (1000*(kStation-1)+1+thisIsFlatback),...
                (1000*(kStation  )+1+nextIsFlatback));
            fprintf(fid,'\n   *get,z_HP_knee_line(%d),line,,num,max',kStation);
            %tcl: create LP_knee line, this is complicated by the fact that 
            %  the KP number jumps from the last KP in station to the first
            if thisIsFlatback
                nCoordPairs = numel(KeyPoints(kStation).x);
                lpKneeKPinboard = (1000*(kStation-1)+nCoordPairs);
            else
                lpKneeKPinboard = (1000*(kStation-1)+1);
            end
            if nextIsFlatback
                nCoordPairs = numel(KeyPoints(kStation+1).x);
                lpKneeKPoutboard = (1000*kStation+nCoordPairs);
            else
                lpKneeKPoutboard = (1000*kStation+1);
            end
            fprintf(fid,'\n   l,%d,%d',lpKneeKPinboard,lpKneeKPoutboard);
            fprintf(fid,'\n   *get,z_LP_knee_line(%d),line,,num,max',kStation);
        end
        %tcl:  create LE line
        fprintf(fid,'\n   l,%d,%d',...
            (1000*(kStation-1)+KeyPoints(kStation).LE),...
            (1000*(kStation)+KeyPoints(kStation+1).LE));
        %tcl:  determine the new line number for LE
        fprintf(fid,'\n   *get,z_LE_line(%d),line,,num,max',kStation);
        fprintf(fid,'\n');
    end
    
    
    %tcl:    Generate areas needed for later LAREA commands (spanwise lines that follow blade envelope)
    fprintf(fid,'\n! Generate areas needed for later LAREA commands =======\n');
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 3'); % turn warning dialog off
    fprintf(fid,'\n*endif');
    fprintf(fid,'\ncsys,0');
    fprintf(fid,'\nallsel\n');
    for kStation = 1:TotalStations-1
        %tcl:  create HP area
        thisIsFlatback = strcmp(app.station(kStation).TEtype,'flat');
        nextIsFlatback = strcmp(app.station(kStation+1).TEtype,'flat');
        if thisIsFlatback || nextIsFlatback
            fprintf(fid,'\n   al,z_HP_line(%d),z_LE_line(%d),z_HP_line(%d),z_HP_knee_line(%d)',kStation,kStation,kStation+1,kStation);
        else
            fprintf(fid,'\n   al,z_HP_line(%d),z_LE_line(%d),z_HP_line(%d),z_TE_line(%d)',kStation,kStation,kStation+1,kStation);
        end
        %tcl:  determine the new area number for HP
        fprintf(fid,'\n   *get,z_HP_area(%d),area,,num,max',kStation);
        %tcl:  create LP area
        if thisIsFlatback || nextIsFlatback
            fprintf(fid,'\n   al,z_LP_line(%d),z_LP_knee_line(%d),z_LP_line(%d),z_LE_line(%d)',kStation,kStation,kStation+1,kStation);
        else
            fprintf(fid,'\n   al,z_LP_line(%d),z_TE_line(%d),z_LP_line(%d),z_LE_line(%d)',kStation,kStation,kStation+1,kStation);
        end
        %tcl:  determine the new area number for LP
        fprintf(fid,'\n   *get,z_LP_area(%d),area,,num,max',kStation);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 2'); % turn warning dialog on
    fprintf(fid,'\n*endif\n');
    
    %tcl:    Generate spanwise area-bounding lines with LAREA command and areas created above
    fprintf(fid,'\n! Generate spanwise area-bounding lines with LAREA command\n');
    for kStation = 1:TotalStations-1
        for kArea = 2:numel(SkinAreas(kStation).startIB)
            %tcl:  section loop variable starts at 2 because TE line was created above
            startIB = KeyPoints(kStation  ).DP(SkinAreas(kStation).startIB(kArea));
            startOB = KeyPoints(kStation+1).DP(SkinAreas(kStation).startOB(kArea));
            startA = 1000*(kStation-1) + startIB;
            startB = 1000*(kStation  ) + startOB;
            
%          # create spanwise line which follows areas created above (blade envelope)
%          if {$TEtype($station) == "flat" && $section == 2} {
%             #                 # seemingly unneccesary but ANSYS is choking on spanwise line generation at the knee of the flatback
%             #                 append shell7Contents "
%             #                    l,$startA,$startB
%             #                 "
%          } else {
            if startIB < KeyPoints(kStation).LE && startOB <= KeyPoints(kStation+1).LE
                %tcl:  HP side
                fprintf(fid,'\n  larea,%d,%d,z_HP_area(%d)',startA,startB,kStation);
            elseif startIB > KeyPoints(kStation).LE && startOB >= KeyPoints(kStation+1).LE
                %tcl:  LP side
                fprintf(fid,'\n  larea,%d,%d,z_LP_area(%d)',startA,startB,kStation);
            elseif startIB == KeyPoints(kStation).LE
                if startOB <= KeyPoints(kStation+1).LE
                    fprintf(fid,'\n  larea,%d,%d,z_HP_area(%d)',startA,startB,kStation);
                else
                    fprintf(fid,'\n  larea,%d,%d,z_LP_area(%d)',startA,startB,kStation);
                end
            else
                error('Inter-Station area boundary lines cannot run across the leading edge line.');
            end
        end
        fprintf(fid,'\n');
    end
    
    
    %tcl:    Cleanup entities generated exclusively to facilitate LAREA command
    fprintf(fid,'\n! Cleanup entities generated exclusively to facilitate LAREA command\n');
    for kStation = 1:TotalStations-1
        fprintf(fid,'\n   adel,z_HP_area(%d)',kStation);
        fprintf(fid,'\n   adel,z_LP_area(%d)',kStation);
        fprintf(fid,'\n   ldel,z_HP_line(%d)',kStation);
        fprintf(fid,'\n   ldel,z_LP_line(%d)',kStation);
        fprintf(fid,'\n   ldel,z_LE_line(%d)',kStation);
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n   ldel,z_HP_line(%d)',TotalStations);
    fprintf(fid,'\n   ldel,z_LP_line(%d)',TotalStations);
    fprintf(fid,'\n');
    
    
    %tcl:    COMBINE LINES IN EACH CROSS SECTION INTO LONGER (COMPLEX) LINES.
    fprintf(fid,'\n! COMBINE LINES IN EACH CROSS SECTION INTO LONGER (COMPLEX) LINES\n');
    % create a Component for grouping the tip station lines
    fprintf(fid,'\n  lsel,none');
    fprintf(fid,'\n  cm,tip_station_lines,line');
    prevLine = 0;
    startLine = 1;
    for kStation = 1:TotalStations
        for kDP = KeyPoints(kStation).DP(2:end-1)
            endLine = prevLine + kDP - 1;
            %tcl:  more than a single simple line, combine
            fprintf(fid,'\n   lsel,s,line,,%d,%d',startLine,endLine);
            if endLine > startLine
                fprintf(fid,'\n   lcomb,all');
            end
            startLine = endLine+1;
            if kStation==TotalStations
                % add line to Component tip_station_lines
                fprintf(fid,'\n   cmsel,a,tip_station_lines');
                fprintf(fid,'\n   cm,tip_station_lines,line');
            end
        end
        
        %tcl:  COMBINE LINES IN LAST PORTION OF AIRFOIL IF LAST PORTION CONTAINS MORE THAN ONE LINE
        endLine = prevLine + numel(KeyPoints(kStation).x);
        if endLine > startLine
            fprintf(fid,'\n   lsel,s,line,,%d,%d',startLine,endLine);
            fprintf(fid,'\n   lcomb,all');
            if kStation==TotalStations
              % add line to Component tip_station_lines
              fprintf(fid,'\n   cmsel,a,tip_station_lines');
              fprintf(fid,'\n   cm,tip_station_lines,line');
            end
            fprintf(fid,'\n');
        end
        prevLine = prevLine + numel(KeyPoints(kStation).x);
        startLine = endLine + 1;
    end
    fprintf(fid,'\nallsel');
    fprintf(fid,'\n!   numcmp,kp');
    fprintf(fid,'\n!   numcmp,line\n');
    
    
    %tcl:    GENERATE SKIN AREAS
    fprintf(fid,'\n! GENERATE SKIN AREAS ==================================\n');
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 3'); % turn warning dialog off
    fprintf(fid,'\n*endif\n');
    for kStation = 1:TotalStations-1
        for kArea = 1:numel(SkinAreas(kStation).startIB)
            %jcb:  need to check on how this implementation compares with
            %original tcl code
            try
            startIB = KeyPoints(kStation  ).DP(SkinAreas(kStation).startIB(kArea));
            startOB = KeyPoints(kStation+1).DP(SkinAreas(kStation).startOB(kArea));
            endIB   = KeyPoints(kStation  ).DP(SkinAreas(kStation).endIB(kArea));
            endOB   = KeyPoints(kStation+1).DP(SkinAreas(kStation).endOB(kArea));
            startA = 1000*(kStation-1) + startIB;
            startB = 1000*(kStation  ) + startOB;
            endA   = 1000*(kStation-1) + endIB;
            endB   = 1000*(kStation  ) + endOB;
            catch ME
                kStation %#ok<NOPRT>
                kArea %#ok<NOPRT>
                rethrow(ME);
            end
            
            fprintf(fid,'\n   asel,none');
            fprintf(fid,'\n   lsel,all');
            fprintf(fid,'\n   csys,0');
            
            %tcl: Test for common coordinate values which can cause funny area generation
            % function cartesian_check(app,startIB,endIB,endOB,startOB)
                %tcl:  Base tolerance on chord length (diameter) of first station
                tolerance = 0.005 * app.station(1).Chord;
                %tcl:  convert kp numbers to pair numbers
                %jcb:  no need to convert the unmodified KPs
                sA = startIB;
                eA = endIB;
                eB = endOB;
                sB = startOB;
                %tcl:  determine x values
                station = app.station(kStation);
                keypts = KeyPoints(kStation);
                x = (keypts.x - station.Xoffset) * station.Chord * twistFlag;
                y = (keypts.y              ) * station.Chord;
                twist = twistFlag * station.DegreesTwist * pi/180;
                xt = cos(twist) * x - sin(twist) * y;
                yt = sin(twist) * x + cos(twist) * y;
                sAx = xt(sA);  eAx = xt(eA); 
                sAy = yt(sA);  eAy = yt(eA);
                
                station = app.station(kStation+1);
                keypts = KeyPoints(kStation+1);
                x = (keypts.x - station.Xoffset) * station.Chord * twistFlag;
                y = (keypts.y              ) * station.Chord;
                twist = twistFlag * station.DegreesTwist * pi/180;
                xt = cos(twist) * x - sin(twist) * y;
                yt = sin(twist) * x + cos(twist) * y;
                sBx = xt(sB);  eBx = xt(eB); 
                sBy = yt(sB);  eBy = yt(eB);
                
                xMin = sAx - tolerance;  xMax = sAx + tolerance;
                yMin = sAy - tolerance;  yMax = sAy + tolerance;
                if ((eAx > xMin && eAx < xMax) && (eBx > xMin && eBx < xMax) && (sBx > xMin && sBx < xMax)) ...
                || ((eAy > yMin && eAy < yMax) && (eBy > yMin && eBy < yMax) && (sBy > yMin && sBy < yMax))
                    %tcl:  similar cartesian coords will cause area generation problem, switch to cylindrical
                    fprintf(fid,'\n   csys,1');
                end

            fprintf(fid,'\n   a,%d,%d,%d,%d',startA,endA,endB,startB);
            
            %tcl:  ASSIGN ATTRIBUTES TO AREA JUST CREATED
            % (1000+kStation) is the csys number for this area
            % (ansysRnumber) is the material real constant number
            % (ansysSecNumber) is the material section number
            switch app.ansys.ElementSystem
                case '91'  %jcb: ToDo - SandwichOption?
                    ansysRnumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
%                     n = strcmp(compsInModel(ansysRnumber),app.matlist);
%                     mat = app.matdb(n);
                    if 1   % {$CompositeSandwichOption($SMN($station,$section)) == 0}
                        fprintf(fid,'\n      aatt,1,%d,17,%d',ansysRnumber,(1000+kStation));
                    else   % {$CompositeSandwichOption($SMN($station,$section)) == 1}
                        fprintf(fid,'\n      aatt,1,%d,20,%d',ansysRnumber,(1000+kStation));
                    end
                case {'99','191'}
                    ansysRnumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
                    n = strcmp(compsInModel(ansysRnumber),app.matlist);
                    mat = app.matdb(n);
                    if isequal(mat.thicknessType,'Constant')
                        fprintf(fid,'\n      aatt,1,%d,12,%d',ansysRnumber,(1000+kStation));
                    elseif  isequal(mat.thicknessType,'Tapered')
                        fprintf(fid,'\n      aatt,1,%d,15,%d',ansysRnumber,(1000+kStation));
                    else
                        error('thicknessType not recognized');
                    end
                case '281'
                    ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
                    fprintf(fid,'\n      aatt,,,31,%d,%d',(1000+kStation),ansysSecNumber);
                case '181'
                    ansysSecNumber = find(strcmp(SkinAreas(kStation).Material(kArea),compsInModel)==1);
                    fprintf(fid,'\n      aatt,,,32,%d,%d',(1000+kStation),ansysSecNumber);
                otherwise
                    errordlg(sprintf('Element System %s not yet available',app.ansys.ElementSystem),'write_shell7 error')
                    error('Element System %s not yet available',app.ansys.ElementSystem);
            end
            
        end
        fprintf(fid,'\n');
          
    end
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 2'); % turn warning dialog on
    fprintf(fid,'\n*endif\n');
    fprintf(fid,'\nallsel');
    fprintf(fid,'\n!   numcmp,kp');
    fprintf(fid,'\n!   numcmp,line');
        
    %tcl:    CREATE AREAS FOR SHEAR WEBS
    fprintf(fid,'\n! CREATE AREAS FOR SHEAR WEBS ==========================\n');
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 3'); % turn warning dialog off
    fprintf(fid,'\n*endif\n');
    for kShearweb = 1:TotalShearwebs
      %tcl:  UNSELECT AREAS TO AVOID THE AATT COMMAND FROM ATTACHING SHEAR 
      %tcl:  WEB ATTRIBUTES WITH SURFACE AREAS OR EVEN OTHER SHEAR WEB AREAS
      fprintf(fid,'\n   asel,none');
      fprintf(fid,'\n   lsel,all');
      fprintf(fid,'\n   csys,0\n');
      
      BeginStation = app.shearweb(kShearweb).BeginStation;
      EndStation = app.shearweb(kShearweb).EndStation;
      firstKP = (BeginStation-1)*1000 + KeyPoints(BeginStation).DP(1+app.shearweb(kShearweb).Corner(1));
      secondKP = (BeginStation-1)*1000 + KeyPoints(BeginStation).DP(1+app.shearweb(kShearweb).Corner(2));
      thirdKP = (EndStation-1)*1000 + KeyPoints(EndStation).DP(1+app.shearweb(kShearweb).Corner(3));
      fourthKP = (EndStation-1)*1000 + KeyPoints(EndStation).DP(1+app.shearweb(kShearweb).Corner(4));
      fprintf(fid,'\n   l,%d,%d',firstKP,secondKP);
      fprintf(fid,'\n   l,%d,%d',thirdKP,fourthKP);
      fprintf(fid,'\n   a,%d,%d,%d,%d',firstKP,secondKP,thirdKP,fourthKP);
      
      %tcl:  ASSIGN ATTRIBUTES NEEDED FOR LATER MESHING
      %tcl:    THE MATERIAL NUMBER 1 (first field in aatt) IS MEANINGLES SINCE SHELL99 IS BEING USED
      % (1000+BeginStation) is the csys number for this area
      % (1000+ansysSecNumber) is the material section number
            switch app.ansys.ElementSystem
                case '91'
                    ansysRnumber = find(strcmp(app.shearweb(kShearweb).Material,compsInModel)==1);
                    %jcb - question: above, the aatt type was 17 or 20
                    fprintf(fid,'\n      aatt,1,%d,16,%d\n',ansysRnumber,(1000+BeginStation));
                case {'99','191'}
                    ansysRnumber = find(strcmp(app.shearweb(kShearweb).Material,compsInModel)==1);
                    n = strcmp(compsInModel(ansysRnumber),app.matlist);
                    mat = app.matdb(n);
                    if isequal(mat.thicknessType,'Constant')
                        %jcb - question: above, the aatt type was 12
                        fprintf(fid,'\n      aatt,1,%d,11,%d',ansysRnumber,(1000+BeginStation));
                    elseif  isequal(mat.thicknessType,'Tapered')
                        %jcb - question: above, the aatt type was 15
                        fprintf(fid,'\n      aatt,1,%d,15,%d',ansysRnumber,(1000+BeginStation));
                    else
                        error('thicknessType not recognized');
                    end
                case '281'
                    ansysSecNumber = find(strcmp(app.shearweb(kShearweb).Material,compsInModel)==1);
                    fprintf(fid,'\n      aatt,,,31,%d,%d\n',(1000+BeginStation),(1000+ansysSecNumber));
                case '181'
                    ansysSecNumber = find(strcmp(app.shearweb(kShearweb).Material,compsInModel)==1);
                    fprintf(fid,'\n      aatt,,,32,%d,%d\n',(1000+BeginStation),(1000+ansysSecNumber));
                otherwise
                    errordlg(sprintf('Element System %s not yet available',app.ansys.ElementSystem),'write_shell7 error')
                    error('Element System %s not yet available',app.ansys.ElementSystem);
            end
        
    end
    fprintf(fid,'\n*if,hide_warndlg_areas,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 2'); % turn warning dialog on
    fprintf(fid,'\n*endif\n');
    
    fprintf(fid,'\n! MESH THE AREAS =======================================\n');
    fprintf(fid,'\nallsel');
    
    %tcl: reverse area normals if clockwise blade
    %tcl:    shear web areas are reversed as well - not necessary, just easier
    if strcmp(app.BladeRotation,'cw')
        fprintf(fid,'\n   areverse,all');
    end
   
    %jcb: are these 2 lines necessary now that we have local coordinate
    %  systems to deal with presweep and precurve?
    fprintf(fid,'\n   local,11,CART,0,0,0,90,0,-90');
    fprintf(fid,'\n   esys,11');
    
    switch app.ansys.meshing
        case 'elementsize'
            fprintf(fid,'\n   mshape,0,2d');
            fprintf(fid,'\n   aesize,all,%f',app.ansys.elementsize);
        case 'smartmesh'
            fprintf(fid,'\n   smrtsize,%f',app.ansys.smartmesh);
    end
    fprintf(fid,'\n   amesh,all');
    fprintf(fid,'\n   csys,0\n');

%    if {$ElementSystem == 191} {
%       # brick elements have been selected, invoke macro
%       append shell7Contents "
%          s2s
%       "
%    }    
    
    %LocationZ_lastStation = app.station(TotalStations).LocationZ;
    %fprintf(fid,'\n   local,12,cart,0,0,%f,0,0,0',LocationZ_lastStation);
    fprintf(fid,'\n   csys,12');
    fprintf(fid,'\n   nsel,none');
    fprintf(fid,'\n   n,,0.0,0.0,0.0');
    fprintf(fid,'\n   *get,z_master_node_number,node,,num,max');
    fprintf(fid,'\n   type,21');
    fprintf(fid,'\n   real,999');
    fprintf(fid,'\n   e,z_master_node_number');
    fprintf(fid,'\n   nsel,all');
    fprintf(fid,'\n   csys,0');
    fprintf(fid,'\n   allsel\n');

%jcb: TO BE DELETED
%     fprintf(fid,'\n   nsel,s,loc,z,%.7f,%.7f\n',...
%         (LocationZ_lastStation - 0.0000001),...
%         (LocationZ_lastStation + 0.0000001));

    % select tip station lines and then nodes attached to those lines
    %jcb: I think this can be cleaned up by moving these after the
    %     'e,z_master_node_number' command above and changing 's' to 'a'
    %     below ('nsel' command is then unnecessary because z_master_node 
    %     is already selected)
    fprintf(fid,'\n   cmsel,s,tip_station_lines');
    fprintf(fid,'\n   nsll,s,1');
    fprintf(fid,'\n   nsel,a,node,,z_master_node_number');
    
    
    switch app.ansys.ElementSystem
        case {'91','99','281','181'}
            fprintf(fid,'\n   cerig,z_master_node_number,all,RXYZ\n');
        case '191'
            fprintf(fid,'\n   cerig,z_master_node_number,all,uxyz\n');
    end
    
    if isequal(app.ansys.BoundaryCondition,'cantilevered')
        %jcb: FIXME - nsel could break with swept/bent blades
        fprintf(fid,'\n   nsel,s,loc,z,0');
        fprintf(fid,'\n   d,all,all');
        fprintf(fid,'\n   nsel,all\n');
    end
    
    fprintf(fid,'\nallsel');
    fprintf(fid,'\n!   nummrg,all');
    fprintf(fid,'\n!   numcmp,node');
    fprintf(fid,'\ncsys,0\n');
    
%    if {$BPEstatus == 1} {
%       append shell7Contents "
%          *get,MaxNode,node,,num,max
%          
%          *dim,nodeNum,,MaxNode
%          *dim,xNode,,MaxNode
%          *dim,yNode,,MaxNode
%          *dim,zNode,,MaxNode
%          *dim,nMass,,MaxNode
%          *dim,D1x,,MaxNode
%          *dim,D1y,,MaxNode
%          *dim,D1z,,MaxNode
%          *dim,D2x,,MaxNode
%          *dim,D2y,,MaxNode
%          *dim,D2z,,MaxNode
%          *dim,D3x,,MaxNode
%          *dim,D3y,,MaxNode
%          *dim,D3z,,MaxNode
%          *dim,D4x,,MaxNode
%          *dim,D4y,,MaxNode
%          *dim,D4z,,MaxNode
%          *dim,D5x,,MaxNode
%          *dim,D5y,,MaxNode
%          *dim,D5z,,MaxNode
%          *dim,D6x,,MaxNode
%          *dim,D6y,,MaxNode
%          *dim,D6z,,MaxNode
%       "
%       # populate node list
%       append shell7Contents "
%          *do,j,1,MaxNode
%             nodeNum(j)=j
%             *get,temp,node,j,loc,x
%             xNode(j)=temp
%             *get,temp,node,j,loc,y
%             yNode(j)=temp
%             *get,temp,node,j,loc,z
%             zNode(j)=temp
%          *enddo
%       "
%    }    
    
    % enter POST1 for postprocessing configuration commands
    fprintf(fid,'\nfinish');
    fprintf(fid,'\n/post1\n');   
    
    fprintf(fid,'\nfctyp,dele,all   ! remove all material failure-criteria postprocessing\n');
    for kfc = 1:length(app.ansys.FailureCriteria)
        if app.ansys.FailureCriteria{kfc,2}
            fprintf(fid,'fctyp,add,%s\n',app.ansys.FailureCriteria{kfc,1});
        end
    end
    fprintf(fid,'\nfinish\n')
    %%% Material Properties %%%
    fprintf(fid,'mpwrite,Materials,txt,,\n');
    if ~all(cellfun('isempty',fcvalues))
        fprintf(fid,'/output,Strengths,txt,,\n');
        fprintf(fid,'TBLIST, ,ALL\n');
        fprintf(fid,'/output\n');
    end
    %%% Section Properties %%%
    fprintf(fid,'/output, Sections,txt\n');
    fprintf(fid,'SLIST,,,,FULL\n');
%     fprintf(fid,'SLIST,\n');
    fprintf(fid,'/output\n');

    %%% Element Properties %%%
    fprintf(fid,'/output, Elements,txt\n');
    fprintf(fid,'elist,all,,,0,0 \n');
    fprintf(fid,'/output\n');
    
    %%% Make NLIST file %%%
    fprintf(fid,'!!! BEGIN MAKE_NLIST MACRO TEXT\n');
    fprintf(fid,'ESEL,S,SEC,,1,999   \n');
    % fprintf(fid,'ALLSEL   \n');
    fprintf(fid,'NSLE,S  \n');
    % fprintf(fid,'nsel,u,node,,z_master_node_number\n');
    fprintf(fid,'/output,NLIST,lis\n');
    fprintf(fid,'/page,1e6,,1e6,,\n');
    fprintf(fid,'NLIST,ALL, , ,XYZ,NODE,NODE,NODE\n');
    fprintf(fid,'/output,\n');
    fprintf(fid,'ALLSEL,ALL \n');

    fprintf(fid,'!!! END MAKE_NLIST TEXT\n');
    
    % save database file
    fprintf(fid,'\nfinish');
    fprintf(fid,'\nsave');

%    # unit tip load suite
%    if {$BPEstatus == 1} {
%       append shell7Contents "    
%    ...
%    ...
    
    
catch ME
    % The try..catch..end ensures the file gets closed
    % in case of a programming error.
    fclose(fid);
    rethrow(ME);
end
fclose(fid);
msg = sprintf('shell7.src written to %s',app.settings.job_path);
if app.batchrun
    disp(msg)
else
    msgbox(msg,'Notification');
end

end

