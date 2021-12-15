%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef BladeDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``BladeDef``  A class definition for wind & water turbine blades.
%
% Example:
%
%     ``blade = BladeDef();``
%
% See also ``BladeDef_to_NuMADfile``, ``xlsBlade``, ``AirfoilDef``, 
% ``StationDef``, ``ComponentDef``, ``StackDef``
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = 'public', GetAccess = 'public')
        aerocenter          % Aerodynamic center of airfoil (used only by NuMAD->FAST)
        chord               % Chord distribution [m] 
        chordoffset         % Chordwise offset (in addition to natural offset)
        components          % Blade components such as spar, panels, etc., refer to ``ComponentDef``
        degreestwist        % Twist distribution [degrees]
        ispan               % Spanwise locations of interpolated output
        leband              % Location of keypoint a
        materials           % Material properties, refer to ``MaterialDef``
        mesh = 0.10;        % Approximate element edge size for FE model [m]
        naturaloffset = 1;  % Integar : 1= offset by max thickness location, 0= do not offset to max thickness
        percentthick        % Percent thickness of airfoil [%]
        prebend             % Blade prebend, reference axis location along x2 [m]
        rotorspin = 1;      % Integar: Rotor Spin, 1= CW rotation looking downwind, -1= CCW rotation
        span                % Spanwise location of distributed properties [m]
        sparcapoffset = 0;  % 1 × 2 array. (Does Nothing)
        sparcapwidth        % Locations of keypoints b & c, defines distance between keypoints b & c [mm]. First entry is the HP spar cap. Second entry is the LP spar cap
        stations            % Blade Stations, define the camber and thickness along the blade, refer to `StationDef``
        sweep               % Blade Sweep, Reference axis location along x1 [m] 
        swtwisted = 0;      % Integar : Shear Web, 0 = planar shear webs, 1= shear webs twisted by blade twist
        teband              % Location of keypoint d        
    end
    properties (SetAccess=private)
        idegreestwist       % (read-only) interpolated twist
        ichord              % (read-only) interpolated chord
        ipercentthick       % (read-only) interpolated thickness
        ichordoffset        % (read-only) interpolated offset
        iaerocenter         % (read-only) interpolated aerocenter
        isweep              % (read-only) interpolated sweep
        iprebend            % (read-only) interpolated prebend
        xoffset             % (read-only) natural offset
        profiles            % (read-only) normalized airfoil profiles
        geometry            % (read-only) actual x,y,z geometry
        arclength           % (read-only) surface distance from L.E.
        cpos                % (read-only) chordwise position
        LEindex
        HParcx0
        LParcx0
        keylabels
        keypoints
        keyarcs
        keycpos
        keyareas
        LEbond
        TEbond
        webindices
        webpoints
        webarcs
        webcpos
        webareas
        webwidth
        webbonds
        bom
        bomIndices
        stacks
        swstacks        
    end
    properties (Hidden, SetAccess=private)
        matdb                  % Composite definition for each region at each station
        TEtype                 % trailing edge type; assigned in updateKeypoints
        shearweb = ...         % shearweb structure used in NuMAD v1
            struct('Material','','BeginStation',[],'EndStation',[],'Corner',[]);              
    end
    properties (Hidden)
        bomPlot = struct('kLayer',1,'hgLinesHP',[],'hgLinesLP',[],...
                         'hgPatchHP',[],'hgPatchLP',[],...
                         'uisliderHP',[],'uisliderLP',[],...
                         'hTitleHP',[],'hTitleLP',[]);
        hgGeometry
        hgKeypoints
        job_name = 'numad.nmd'; % (hidden) NuMAD file name -- needed??
        paths = struct('job','','numad','','precomp','','bmodes','','ansys','','batch_run',0)
%         job_path              % (hidden) working NuMAD directory
%         numad_path            % (hidden) NuMAD code directory location
%         ansys_path            % (hidden) ANSYS file path location
%         precomp_path          % (hidden) PreComp .exe file path location
%         bmodes_path           % (hidden) BModes .exe file path location

        % generate ANSYS settings
        ansys = struct('BoundaryCondition','','ElementSystem','','MultipleLayerBehavior','',...
            'meshing','','smartsize',[],'elementsize',[],'shell7gen',[],...
            'dbgen',[],'FailureCriteria',[])
    end
    
    methods
        function obj = BladeDef()
            if nargin > 0
                
            end
            obj.checkNaturalOffset
            obj.checkRotorSpin
            obj.checkSwtwisted
        end
        
        function checkNaturalOffset(obj)
            % This method checks ``naturaloffset`` values
            % 
            if ~(isequal(obj.naturaloffset,0) || isequal(obj.naturaloffset,1))
                error('naturaloffset must be 0 or 1');
            end
        end
        
        function checkRotorSpin(obj)
            % This method checks ``rotorspin`` values
            % 
            if ~(isequal(obj.rotorspin,1) || isequal(obj.rotorspin,-1))
                error('rotorspin must be 1 (cw) or -1 (ccw)');
            end
        end
        
        function checkSwtwisted(obj)
            % This method checks ``swtwisted`` values
            % 
            if ~(isequal(obj.swtwisted,0) || isequal(obj.swtwisted,1))
                error('swtwisted must be 0 or 1');
            end
        end
        
        function addStation(obj,af,spanlocation)
            % This method adds a station
            % 
            % Example:
            %            
            %   ``blade.addStation(af,spanlocation)`` where  ``af`` = airfoil filename or ``AirfoilDef`` object
            %             
            N = numel(obj.stations);
            k = N + 1;
            if k>1
                obj.stations(k) = StationDef(af);
            else
                obj.stations    = StationDef(af);
            end
            obj.stations(k).spanlocation = spanlocation;
            obj.stations(k).parent = obj;
        end
        
        function addComponent(obj,comp)
            % This method adds a Component
            % 
            % Example:
            %            
            %   ``blade.addComponent(comp_struct)`` where ``comp_struct`` = input structure used by ``ComponentDef``
            %               
            N = numel(obj.components);
            k = N + 1;
            if k>1
                obj.components(k) = ComponentDef(comp);
            else
                obj.components    = ComponentDef(comp);
            end
        end
        
        function addMaterial(obj,mat)
            % This method adds Material
            % 
            % Example:
            %            
            %   ``blade.addMaterial(mat_struct)`` where ``mat_struct`` = input structure used by ``MaterialDef``
            %   
            N = numel(obj.materials);
            k = N + 1;
            if k>1
                obj.materials(k) = MaterialDef(mat);
            else
                obj.materials    = MaterialDef(mat);
            end
        end
        
        function updateBlade(obj)
            % This method updates the BladeDef model
            obj.updateGeometry;
            obj.updateKeypoints;
            obj.updateBOM;
        end
            
        function updateGeometry(obj)
            % This method updates the interpolated blade parameters
            
            % update the interpolated station profiles
            nStations = numel(obj.stations);
            if nStations > 0
                nPoints = length(obj.stations(1).airfoil.c);
            else
                error('BladeDef must have at least one station before updating geometry.');
            end
            
            % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % add some error checking -- first station must be at blade
            % root to prevent extrapolation
% %             assert(obj.stations(1).spanlocation==0,'first station must be at the blade root')
            % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
            % Collect parameter tables from the stations.
            spanlocation = [obj.stations.spanlocation];
            c = zeros(nPoints,nStations);
            camber = zeros(nPoints,nStations);
            thickness = zeros(nPoints,nStations);
            tetype = cell(1,nStations);
            for k=1:nStations
                ck = obj.stations(k).airfoil.c;
                if length(ck)~=nPoints
                    error('Station airfoils must have same number of samples.')
                end
                c(:,k) = ck;
                camber(:,k) = obj.stations(k).airfoil.camber;
                thickness(:,k) = obj.stations(k).airfoil.thickness;
                tetype{k} = obj.stations(k).airfoil.TEtype;
            end
            
            % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % fix numerical issue due to precision on camber calculation
            % camber should start and end at y-values of zero
            camber(1,:)   = zeros(1,nStations);
            camber(end,:) = zeros(1,nStations);
            % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
            % Interpolate the station parameter tables.
            % Each column corresponds to an interpolated station.
            ic = transpose(interp1(spanlocation,c',obj.ispan,'pchip'));
            icamber = transpose(interp1(spanlocation,camber',obj.ispan,'pchip'));
            ithickness = transpose(interp1(spanlocation,thickness',obj.ispan,'pchip'));
            obj.cpos = [ -ic(end,:); -flipud(ic);   ic(2:end,:);  ic(end,:)];  % chordwise position of interpolated station points
%             figure(101); surf(repmat(spanlocation,nPoints,1),c,camber,'MeshStyle','column');
%             figure(102); surf(repmat(spanlocation,nPoints,1),c,thickness,'MeshStyle','column');
%             figure(103); surf(repmat(obj.ispan,nPoints,1),ic,icamber,'MeshStyle','column');
%             figure(104); surf(repmat(obj.ispan,nPoints,1),ic,ithickness,'MeshStyle','column');
                                    
            % Adjust the thickness profiles based on TEtype of stations. 
            % This is mainly for transitions to flatbacks were the
            % interpolated airfoil needs to look like a round.
            for k=1:numel(obj.ispan)
                ind = find(obj.ispan(k)<spanlocation,1);
                if isempty(ind) || ind==1
                    continue;
                end
                if isequal(tetype{ind},'flat') && ...
                    isequal(tetype{ind-1},'round')
                        ithickness(end,k) = 0;
                end
            end

            % Interpolate the blade parameter curves.
            % Results are row vectors.
            obj.idegreestwist  = interp1(obj.span,obj.degreestwist,obj.ispan,'pchip');
            obj.ichord         = interp1(obj.span,obj.chord       ,obj.ispan,'pchip');
                absolutethick  = obj.percentthick .* obj.chord  / 100;
                iabsolutethick = interp1(obj.span,absolutethick  ,obj.ispan,'pchip');
            obj.ipercentthick  =   iabsolutethick ./ obj.ichord * 100;
                % ensure that the interpolation doesn't reduce the percent
                % thickness beneath the thinnest airfoil
                obj.ipercentthick(obj.ipercentthick < min(obj.percentthick)) = min(obj.percentthick);
            obj.ichordoffset   = interp1(obj.span,obj.chordoffset ,obj.ispan,'pchip');
            obj.iaerocenter    = interp1(obj.span,obj.aerocenter  ,obj.ispan,'pchip');
            if isempty(obj.sweep),   obj.sweep = zeros(size(obj.span)); end
            if isempty(obj.prebend), obj.prebend = zeros(size(obj.span)); end
            obj.isweep         = interp1(obj.span,obj.sweep       ,obj.ispan,'pchip');
            obj.iprebend       = interp1(obj.span,obj.prebend     ,obj.ispan,'pchip');
            
            % Generate the blade surface geometry.
            N = numel(obj.ispan);
            M = nPoints*2 + 1;
            obj.profiles = zeros(M,2,N);  % 2D normalized profiles
            obj.geometry = zeros(M,3,N);  % actual 3D geometry (scaled, rotated, etc.)
            obj.xoffset = zeros(1,N);
            obj.LEindex = nPoints + 1;  % index of leading edge
            for k=1:N
                %jcb: note that I'm using max thickness about camber
                %instead of overall thickness of airfoil. We may need to
                %change this definition.
                [maxthick,mtindex] = max(ithickness(:,k));
                tratio = obj.ipercentthick(k) / (maxthick * 100);
                thick = ithickness(:,k) * tratio;
                hp = icamber(:,k) - 0.5*thick;
                lp = icamber(:,k) + 0.5*thick;
                c = ic(:,k);
                x = [ c(end); flipud(c) ;  c(2:end);  c(end)];
                y = [ 0     ; flipud(hp); lp(2:end);  0     ];
                obj.profiles(:,1,k) = x;
                obj.profiles(:,2,k) = y;                
                obj.xoffset(1,k) = c(mtindex);
                if obj.naturaloffset
                    x = x - c(mtindex); 
                end
                x = x - obj.ichordoffset(k);   % apply chordwise offset
                x = x * obj.ichord(k) * -1*obj.rotorspin;  % scale by chord
                y = y * obj.ichord(k);                     % scale by chord
                twist = -1*obj.rotorspin * obj.idegreestwist(k);
                % prepare for hgtransform rotate & translate
                coords(:,1) = cosd(twist) * x - sind(twist) * y;
                coords(:,2) = sind(twist) * x + cosd(twist) * y;
                coords(:,3) = zeros(size(x));
                coords(:,4) = ones(size(x));
                
                % use the generating line to translate and rotate the coordinates
                [sweep_rot, prebend_rot] = deal(0);
% jcb: This code, copied from NuMAD 2.0, causes each section to rotate out
% of plane so that its normal follows the generating line direction. Need
% to replace 'twistFlag' with '-1*obj.rotorspin' and calculate the slopes
% based on the available data. For now, default to parallel sections.
%                 if isequal(blade.PresweepRef.method,'normal')
%                     sweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
%                     sweep_rot = atan(sweep_slope*twistFlag);
%                 end
%                 if isequal(blade.PrecurveRef.method,'normal')
%                     prebend_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
%                     prebend_rot = atan(-prebend_slope);
%                 end
                transX = -1*obj.rotorspin*obj.isweep(k);
                transY = obj.iprebend(k);
                transZ = obj.ispan(k);
                R = makehgtform('yrotate',sweep_rot,'xrotate',prebend_rot);
                T = makehgtform('translate',transX,transY,transZ);
                coords = coords * R' * T';
                % save the transformed coordinates
                obj.geometry(:,1,k) = coords(:,1);
                obj.geometry(:,2,k) = coords(:,2);
                obj.geometry(:,3,k) = coords(:,3);
            end
            
            % Calculate the arc length of each curve
            obj.arclength = zeros(M,N);
            obj.HParcx0 = zeros(1,N);
            obj.LParcx0 = zeros(1,N);
%             ptind = 1:M;  % point indices
%             ptindover = 1:0.2:M;  % over-sample point indices
%             LE = find(ptindover==obj.LEindex);
            LE = obj.LEindex;
            for k=1:N
%                 xx = spline(ptind,obj.geometry(:,1,k),ptindover);
%                 yy = spline(ptind,obj.geometry(:,2,k),ptindover);
%                 zz = spline(ptind,obj.geometry(:,3,k),ptindover);
                xx = obj.geometry(:,1,k);
                yy = obj.geometry(:,2,k);
                zz = obj.geometry(:,3,k);

                arclen = sqrt(diff(xx).^2 + diff(yy).^2 + diff(zz).^2);
                arclen = [0; cumsum(arclen(:))];
%                 obj.arclength(:,k) = interp1(ptindover,arclenover,ptind);
                obj.arclength(:,k) = arclen;
                LEarcsum = obj.arclength(obj.LEindex,k);
                obj.arclength(:,k) = obj.arclength(:,k) - LEarcsum;
                % find where x=0 intersects the surface
                obj.HParcx0(1,k) = interp1(xx(2:LE),arclen(2:LE),0) - LEarcsum;
                obj.LParcx0(1,k) = interp1(xx(end-1:-1:LE),arclen(end-1:-1:LE),0) - LEarcsum;
            end
            
        end
        
        function updateKeypoints(obj)
            % This method updates the keypoints (a,b,c,...) which define the blade
            % regions.
            % 
            % Example:
            %            
            %   ``blade.updateKeypoints``            
            %               
                        
            % find the curves which bound each blade region
            N = numel(obj.ispan);   % number of interpolated span stations
            M = 12;      % number of areas around airfoil profile; must be even (see calc of web areas)
            obj.keypoints = zeros(M-2,3,N);  % keypoints in xyz geometry
            obj.keyarcs   = zeros(M+1,N);    % surface arclength distance of keypoints from LE 
            obj.keycpos   = zeros(M+1,N);    % chordwise position of keypoints
            obj.keyareas  = zeros(M,N-1);    % surface area of regions created by keypoints
            obj.keylabels = {'te'; 'e'; 'd'; 'c'; 'b'; 'a'; 'le'; ...
                              'a'; 'b'; 'c'; 'd'; 'e'; 'te'};
            obj.LEbond = zeros(1,N-1);
            obj.TEbond = zeros(1,N-1);
            mm_to_m = 1e-3;
            ns=2; nf=size(obj.geometry,1)-1;  % start and finish indices in geometry/arcs
            for k=1:N
                % allow for separate definitions of HP and LP spar cap
                % width and offset [HP LP]
                if length(obj.sparcapwidth) > 2 || length(obj.sparcapoffset) > 2
                    error('too many entries for spar cap definition')
                end
                scwidth_hp = mm_to_m*obj.sparcapwidth(1);
                scwidth_lp = mm_to_m*obj.sparcapwidth(end);
                scoffset_hp = mm_to_m*obj.sparcapoffset(1);
                scoffset_lp = mm_to_m*obj.sparcapoffset(end);

                n1 = mm_to_m*obj.leband;  % no foam width
                n2 = mm_to_m*obj.teband;  % no foam width
                
                tempTE = obj.getprofileTEtype(k);
                obj.TEtype{k} = tempTE{1};
                if obj.swtwisted
                    % get angle of each xy pair w.r.t. pitch axis (0,0)
                    xyangle = zeros(size(obj.geometry,1),1);
                    for j=1:length(xyangle)
                        xy=obj.geometry(j,1:2,k);
                        xyangle(j) = atan2(obj.rotorspin*xy(2),xy(1));
                    end
                    % unwrap and center around 0
                    xyangle = unwrap(xyangle);
                    xyangle = xyangle - pi * round(xyangle(obj.LEindex)/pi);
                end
                % ****************************************************
                % ==================== HP surface ====================
                if obj.swtwisted
                    % find arclength where xyangle equals normal to chord
                    twistnorm = pi/180*(-obj.idegreestwist(k) - 90);  % angle normal to chord line
                    z = interp1(xyangle(ns:nf),obj.arclength(ns:nf,k),twistnorm);
                else
                    z = obj.HParcx0(1,k);
                end
                z0 = z;            % ble: location where airfoil surface crosses Xglobal=0
                z = z - scoffset_hp;  % positive scoffset moves z toward t.e.
                a = max( (0 - n1)                 , 0.10*obj.arclength(ns,k));
                a = min(a, 0.01*obj.arclength(ns,k));
                b = min( (z + 0.5*scwidth_hp)        , 0.15*obj.arclength(ns,k));
                c = max( (z - 0.5*scwidth_hp)        , 0.80*obj.arclength(ns,k));
                d = min( (obj.arclength(1,k) + n2), 0.85*obj.arclength(ns,k));
                d = max(d, 0.98*obj.arclength(ns,k));
                if strcmp(obj.TEtype(k),'flat')
                    e = obj.arclength(ns,k);
                    obj.keypoints( 1,:,k) = obj.geometry(ns,:,k);
                    obj.keycpos( 2,k) = -1;
                else
%                     e = 0.5 * (d + obj.arclength(ns,k));
                    e = 0.99*obj.arclength(ns,k);
                    obj.keypoints( 1,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),e);
                    obj.keycpos( 2,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),e);
                end
                %              1     -> e
                obj.keypoints( 2,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),d);
                obj.keypoints( 3,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),c);
%                 obj.keypoints(  ,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),z);
                obj.keypoints( 4,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),b);
                obj.keypoints( 5,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),a);
                obj.keyarcs( 1,k)   = obj.arclength(ns,k);
                obj.keyarcs( 2,k)   = e;
                obj.keyarcs( 3,k)   = d;
                obj.keyarcs( 4,k)   = c;
%                 obj.keyarcs(  ,k)   = z;
                obj.keyarcs( 5,k)   = b;
                obj.keyarcs( 6,k)   = a;
                obj.keyarcs( 7,k)   = 0; % le
                obj.keycpos( 1,k) = obj.cpos(ns,k);  % te, hp surface
                %            2   -> e
                obj.keycpos( 3,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),d);
                obj.keycpos( 4,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),c);
%                 obj.keycpos(  ,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),z);
                obj.keycpos( 5,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),b);
                obj.keycpos( 6,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),a);
                obj.keycpos( 7,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),0);  % le
                % ****************************************************
                % ==================== LP surface ====================
                if obj.swtwisted
                    twistnorm = pi/180*(-obj.idegreestwist(k) + 90);  % angle normal to chord line
                    z = interp1(xyangle(ns:nf),obj.arclength(ns:nf,k),twistnorm);
                else
                    z = obj.LParcx0(1,k);
                end
                z0 = z;            % ble: location where airfoil surface crosses Xglobal=0
                z = z + scoffset_lp;  % positive scoffset moves z toward t.e.
                a = min( (0 + n1)                   , 0.10*obj.arclength(nf,k));
                a = max(a, 0.01*obj.arclength(nf,k));
                b = max( (z - 0.5*scwidth_lp)          , 0.15*obj.arclength(nf,k));
                c = min( (z + 0.5*scwidth_lp)          , 0.80*obj.arclength(nf,k));
                d = max( (obj.arclength(end,k) - n2), 0.85*obj.arclength(nf,k));
                d = min(d, 0.96*obj.arclength(nf,k));
                if strcmp(obj.TEtype(k),'flat')
                    e = obj.arclength(nf,k);
                    obj.keypoints(10,:,k) = obj.geometry(nf,:,k);
                    obj.keycpos(12,k) = 1;
                else
%                     e = 0.5 * (d + obj.arclength(nf,k));
                    e = 0.98 * obj.arclength(nf,k);
                    obj.keypoints(10,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),e);
                    obj.keycpos(12,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),e);
                end
                obj.keypoints( 6,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),a);
                obj.keypoints( 7,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),b);
%                 obj.keypoints(  ,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),z);
                obj.keypoints( 8,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),c);
                obj.keypoints( 9,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),d);
                %             10   -> e
                obj.keyarcs( 8,k)   = a;
                obj.keyarcs( 9,k)   = b;
%                 obj.keyarcs( ,k)   = z;
                obj.keyarcs(10,k)   = c;
                obj.keyarcs(11,k)   = d;
                obj.keyarcs(12,k)   = e;
                obj.keyarcs(13,k)   = obj.arclength(nf,k);
                obj.keycpos( 8,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),a);
                obj.keycpos( 9,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),b);
%                 obj.keycpos(  ,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),z);
                obj.keycpos(10,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),c);
                obj.keycpos(11,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),d); 
                %           12   -> e
                obj.keycpos(13,k) = obj.cpos(nf,k);  % te, lp surface
            end
            
            % find the points used by each shear web
            obj.webindices = cell(0);
            obj.webarcs    = cell(0);
            obj.webcpos    = cell(0);
            obj.webpoints  = cell(0);
            obj.webareas   = cell(0);
            obj.webwidth   = cell(0);
            cmpt_groups = [obj.components.group];
            uniq_groups = unique(cmpt_groups);
            uniq_groups(uniq_groups==0) = [];   % group "0" is the blade skins
            for ksw = uniq_groups  % for each shear web
                ksw_cmpts = find(ksw == cmpt_groups);  % find the components that are part of the shear web
                hpextents = unique([obj.components(ksw_cmpts).hpextents]); % get the hp extents
                lpextents = unique([obj.components(ksw_cmpts).lpextents]); % get the lp extents
                assert(numel(hpextents)==1,'HP Extents for components in group %d must be identical and contain no spaces or commas',ksw);
                assert(numel(lpextents)==1,'LP Extents for components in group %d must be identical and contain no spaces or commas',ksw);
                % match extents that have form of either '0.5b-c' or
                % 'b+/-100' or 'b' or 'z+/-100'
                pat = '(?<fraction>\d*[\.]?\d*)(?<pt1>[a-zA-Z]+)-(?<pt2>[a-zA-Z]+)|(?<pt3>[a-zA-Z]+)(?<mm_offset>[+-]\d+)|(?<pt>[a-zA-Z])'; 
                hp = regexp(hpextents{1},pat,'names');
                lp = regexp(lpextents{1},pat,'names');
                le = find(1==strcmpi('le',obj.keylabels));
                % get shear web placement on HP side
                if ~isempty(hp.pt)
                    n = find(1==strcmpi(hp.pt,obj.keylabels(1:le))) + 0;
                    assert(~isempty(n),'HP extent label "%s" not defined.',hp.pt);
                    obj.webindices{ksw}(1,1) = n;
                    obj.webarcs{ksw}(1,:) = obj.keyarcs(n,:);
                    obj.webcpos{ksw}(1,:) = obj.keycpos(n,:);
					n = find(1==strcmpi(hp.pt,obj.keylabels(1:le))) - 1; %% EMA
                    obj.webpoints{ksw}(1,:,:) = obj.keypoints(n,:,:);
                elseif ~isempty(hp.pt1)
                    f = str2double(hp.fraction);
                    if f<=0 || f>=1
                        error('Component group %d: HP extent fraction=%g, which is outside range (0..1)',ksw,f);
                    end
                    n1 = find(1==strcmpi(hp.pt1,obj.keylabels(1:le))) + 0;
                    n2 = find(1==strcmpi(hp.pt2,obj.keylabels(1:le))) + 0;
                    assert(~isempty(n1),'HP extent label "%s" not defined.',hp.pt1);
                    assert(~isempty(n2),'HP extent label "%s" not defined.',hp.pt2);
                    obj.webindices{ksw}(1,1) = nan;
                    p1 = obj.keyarcs(n1,:);
                    p2 = obj.keyarcs(n2,:);
                    p  = (1-f)*p1 + f*p2;
                    obj.webarcs{ksw}(1,:) = p;
                    for k=1:N
                        obj.webcpos{ksw}(1,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),p(k));
                        obj.webpoints{ksw}(1,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),p(k));
                    end
                elseif ~isempty(hp.pt3)
                    n3 = find(1==strcmpi(hp.pt3,obj.keylabels(1:le)))+0;
                    assert(~isempty(n3),'HP extent label "%s" not defined.',hp.pt3);                    
                    obj.webindices{ksw}(1,1) = nan; % web location is not at an existing keypoint
                    p3 = obj.keycpos(n3,:);   % chordwise position of keypoint hp.pt3
                    p  = p3 - str2double(hp.mm_offset)/1000;  % chordwise position of the shear web
                    iMax = find(1==strcmpi('d',obj.keylabels(1:le)))+0;  % index location of keypoint d on HP side
                    pMax = obj.keycpos(iMax,:).*obj.ichord'; % the maximum chord location
                    p(abs(p)>abs(pMax)) = pMax(abs(p)>abs(pMax));
                    iMin = find(1==strcmpi('a',obj.keylabels(1:le)))+0;  % index location of keypoint a on HP side
                    pMin = obj.keycpos(iMin,:).*obj.ichord'; % the minimum chord location
                    p(abs(p)<abs(pMin)) = pMin(abs(p)<abs(pMin));                    
                    obj.webcpos{ksw}(1,:) = p;
                    for k=1:N
                        obj.webarcs{ksw}(1,k) = interp1(obj.cpos(ns:nf,k),obj.arclength(ns:nf,k),p(k));
                        obj.webpoints{ksw}(1,:,k) = interp1(obj.cpos(ns:nf,k),obj.geometry(ns:nf,:,k),p(k));
                    end                    
                    % ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    p3Save{ksw}(1,:) = -1.*p3;
                    pSave{ksw}(1,:) = -1.*p;
                    % ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                else
                    error('Shear web geometry HP extents not defined correctly (e.g., 0.5b-c, b, b+200)')
                end
                % get shear web placement on LP side
                if ~isempty(lp.pt)    
                    n = find(1==strcmpi(lp.pt,obj.keylabels(le:end))) + le - 1;
                    assert(~isempty(n),'LP extent label "%s" not defined.',lp.pt);
                    obj.webindices{ksw}(2,1) = n;
                    obj.webarcs{ksw}(2,:) = obj.keyarcs(n,:);
                    obj.webcpos{ksw}(2,:) = obj.keycpos(n,:);
                    obj.webpoints{ksw}(2,:,:) = obj.keypoints(n,:,:);
                elseif ~isempty(lp.pt1)    
                    f = str2double(lp.fraction);
                    if f<0 || f>1
                        error('Component group %d: LP extent fraction=%g, which is outside range [0..1]',ksw,f);
                    end
                    n1 = find(1==strcmpi(lp.pt1,obj.keylabels(le:end))) + le - 1;
                    n2 = find(1==strcmpi(lp.pt2,obj.keylabels(le:end))) + le - 1;
                    assert(~isempty(n1),'LP extent label "%s" not defined.',lp.pt1);
                    assert(~isempty(n2),'LP extent label "%s" not defined.',lp.pt2);
                    obj.webindices{ksw}(2,1) = nan;
                    p1 = obj.keyarcs(n1,:);
                    p2 = obj.keyarcs(n2,:);
                    p  = (1-f)*p1 + f*p2;
                    obj.webarcs{ksw}(2,:) = p;
                    for k=1:N
                        obj.webcpos{ksw}(2,k) = interp1(obj.arclength(ns:nf,k),obj.cpos(ns:nf,k),p(k));
                        obj.webpoints{ksw}(2,:,k) = interp1(obj.arclength(ns:nf,k),obj.geometry(ns:nf,:,k),p(k));
                    end
                elseif ~isempty(lp.pt3)
                    n3 = find(1==strcmpi(lp.pt3,obj.keylabels(le:end))) + le - 1;
                    assert(~isempty(n3),'LP extent label "%s" not defined.',lp.pt3);                    
                    obj.webindices{ksw}(2,1) = nan; % web location is not at an existing keypoint
                    p3 = obj.keycpos(n3,:);   % chordwise position of keypoint lp.pt3
                    p  = p3 + str2double(lp.mm_offset)/1000;  % chordwise position of the shear web
                    iMax = find(1==strcmpi('d',obj.keylabels(le:end))) + le - 1;  % index location of keypoint d on LP side
                    pMax = obj.keycpos(iMax,:).*obj.ichord'; % the maximum chord location
                    p(abs(p)>abs(pMax)) = pMax(abs(p)>abs(pMax));                    
                    iMin = find(1==strcmpi('a',obj.keylabels(le:end))) + le - 1;  % index location of keypoint d on LP side
                    pMin = obj.keycpos(iMin,:).*obj.ichord'; % the maximum chord location
                    p(abs(p)<abs(pMin)) = pMin(abs(p)<abs(pMin));  
                    obj.webcpos{ksw}(2,:) = p;
                    for k=1:N
                        obj.webarcs{ksw}(2,k) = interp1(obj.cpos(ns:nf,k),obj.arclength(ns:nf,k),p(k));
                        obj.webpoints{ksw}(2,:,k) = interp1(obj.cpos(ns:nf,k),obj.geometry(ns:nf,:,k),p(k));
                    end   
                    % ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                    p3Save{ksw}(2,:) = p3;
                    pSave{ksw}(2,:) = p;
                    % ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                else
                    error('Shear web geometry LP extents not defined correctly (e.g., 0.5b-c, b, b+200)')
                end
            end
            % calculate shell areas
            for kc = 1:N-1
                for kr = 1:M
                    % choose number of points to use in area calculation
                    % jcb: I decided to base this on the number of points
                    % in the interpolated station profile found within the
                    % region of interest.
                    npts = sum( obj.arclength(:,kc) >= obj.keyarcs(kr  ,kc) & ...
                                obj.arclength(:,kc) <= obj.keyarcs(kr+1,kc));
                    npts = max(npts,2); % need at least two points
                    ibarc = linspace(obj.keyarcs(kr  ,kc),...
                                     obj.keyarcs(kr+1,kc),...
                                     npts); % inboard curve arclengths
                    obarc = linspace(obj.keyarcs(kr  ,kc+1),...
                                     obj.keyarcs(kr+1,kc+1),...
                                     npts); % outboard curve arclengths
                    ib = interp1(obj.arclength(ns:nf,kc  ),obj.geometry(ns:nf,:,kc  ),ibarc);  % inboard xyz
                    ob = interp1(obj.arclength(ns:nf,kc+1),obj.geometry(ns:nf,:,kc+1),obarc);  % outboard xyz
                    dspan = sqrt(sum((ob-ib).^2,2));  % "ds" in the span direction
                    obj.keyareas(kr,kc) = ... % treat each "rectangular" area as two triangles
                        0.5*sqrt(sum(diff(ib).^2,2))' * dspan(1:end-1) + ...
                        0.5*sqrt(sum(diff(ob).^2,2))' * dspan(2:end);  
                    if isequal(1,kr)
                        obj.TEbond(kc) = dspan(1);
                    end
                    if isequal(M/2+1,kr)
                        obj.LEbond(kc) = dspan(1);
                    end
                end
            end
            
            % calculate areas used by shear webs
            % jcb: note that these areas come purely from the geometry and
            % do not take into account the thickness of the shell or 
            % sparcap layup.
            for ksw = 1:numel(obj.webpoints)
                for kc = 1:N-1
                    ib = obj.webpoints{ksw}(:,:,kc  );
                    ob = obj.webpoints{ksw}(:,:,kc+1);
                    % treat each "rectangular" area as two triangles
                    b1 = diff(ib);   % vector between inboard points
                    b2 = diff(ob);   % vector between outboard points
                    base1 = sqrt(sum(b1.^2,2));   % length of triangle base
                    base2 = sqrt(sum(b2.^2,2));   % length of triangle base
                    b1 = b1 / base1;  % unit vector
                    b2 = b2 / base2;  % unit vector
                    h1 = abs((ob(1,:)-ib(1,:)) * (1 - b1'));   % height of triangle (perpendicular distance to base)
                    h2 = abs((ib(2,:)-ob(2,:)) * (1 - b2'));   % height of triangle (perpendicular distance to base)
                    obj.webareas{ksw}(1,kc) = 0.5 * (base1*h1 + base2*h2);
                    obj.webwidth{ksw}(1,kc) = base1;
                    % calculate edge (bond-line) lengths
                    obj.webbonds{ksw}(1:2,kc) = sqrt(sum((ob-ib).^2,2));
                end
                obj.webwidth{ksw}(1,N) = base2;
            end
        end
        
        function updateBOM(obj)
            % This method updates the Bill-of-Materials
            % Cell array columns of blade.bom: 
            % 1. Layer #
            % 2. Material ID
            % 3. Component or region name
            % 4. Begin station   (m)
            % 5. End station     (m)
            % 6. Max width       (m)
            % 7. Average width   (m)
            % 8. 3D area         (m^2)
            % 9. Layer thickness (mm)
            % 10. Computed dry layer weight (g)
            %
            
            % initialize structures
            obj.bom = struct('hp',cell(0),'lp',cell(0),'sw',cell(0),...
                             'lebond',[],'tebond',[],'swbonds',[],...
                             'dryweight',[]);
            obj.bomIndices = struct('hp',[],'lp',[],'sw',cell(0));
            % calculate non-dimensional span
            ndspan = (obj.ispan - obj.ispan(1)) ./ (obj.ispan(end) - obj.ispan(1));
            hprow = 1;
            lprow = 1;
            swnum = 0;
            swrow = 1;
            g_to_kg = 1e-3;
            m_to_mm = 1e3;
            mm_to_m = 1e-3;
            nComponents = numel(obj.components);
            for kc = 1:nComponents
                comp = obj.components(kc);
                mat = obj.materials(comp.materialid);
                [hpRegion, lpRegion] = obj.findRegionExtents(obj.keylabels,comp);
                nlay = comp.getNumLayers(ndspan);
                nlay = round(nlay);
%                 if 1
                    layermult = 1;        % multiplier=1 for individual layers
%                 else
%                     layermult = nlay(1);  % multiplies layer thickness
%                     assert(all(nlay==nlay(1)),'updateBOM: component %d, ''T'' flag requires uniform layer thickness',kc);
%                     nlay = ones(size(ndspan));
%                 end
                for klay=1:max(nlay)
                    [beginSta, endSta] = obj.findLayerExtents(nlay,klay);
                    %% EMA original:
%                     for ks = 1:length(beginSta)
                    %% changed to :
                    ksMax = min([length(beginSta),length(endSta)]);
                    for ks = 1:ksMax
                    %% END
                        if comp.group==0 && ~isempty(hpRegion)
                        areas =  obj.keyareas(hpRegion(1):hpRegion(2)-1, beginSta(ks):endSta(ks)-1);
                        regionarea = sum(areas(:));
                        arcs = obj.keyarcs(hpRegion(2), beginSta(ks):endSta(ks)) ...
                             - obj.keyarcs(hpRegion(1), beginSta(ks):endSta(ks));
                        obj.bom(1).hp{hprow, 1} = hprow;
                        obj.bom(1).hp{hprow, 2} = comp.materialid;
                        obj.bom(1).hp{hprow, 3} = comp.name;
                        obj.bom(1).hp{hprow, 4} = obj.ispan(beginSta(ks));
                        obj.bom(1).hp{hprow, 5} = obj.ispan(endSta(ks));
                        obj.bom(1).hp{hprow, 6} = max(arcs);
                        obj.bom(1).hp{hprow, 7} = mean(arcs);
                        obj.bom(1).hp{hprow, 8} = regionarea;
                        obj.bom(1).hp{hprow, 9} = layermult * mat.layerthickness;
                        obj.bom(1).hp{hprow,10} = layermult * mat.drydensity * regionarea;
                        obj.bomIndices(1).hp(hprow,1:4) = [beginSta(ks), endSta(ks), hpRegion];
                        hprow = hprow + 1;
                        end
                        if comp.group==0 && ~isempty(lpRegion)
                        areas =  obj.keyareas(lpRegion(1):lpRegion(2)-1, beginSta(ks):endSta(ks)-1);
                        regionarea = sum(areas(:));
                        arcs = obj.keyarcs(lpRegion(2), beginSta(ks):endSta(ks)) ...
                             - obj.keyarcs(lpRegion(1), beginSta(ks):endSta(ks));
                        obj.bom(1).lp{lprow, 1} = lprow;
                        obj.bom(1).lp{lprow, 2} = comp.materialid;
                        obj.bom(1).lp{lprow, 3} = comp.name;
                        obj.bom(1).lp{lprow, 4} = obj.ispan(beginSta(ks));
                        obj.bom(1).lp{lprow, 5} = obj.ispan(endSta(ks));
                        obj.bom(1).lp{lprow, 6} = max(arcs);
                        obj.bom(1).lp{lprow, 7} = mean(arcs);
                        obj.bom(1).lp{lprow, 8} = regionarea;
                        obj.bom(1).lp{lprow, 9} = layermult * mat.layerthickness;
                        obj.bom(1).lp{lprow,10} = layermult * mat.drydensity * regionarea;
                        obj.bomIndices(1).lp(lprow,1:4) = [beginSta(ks), endSta(ks), lpRegion];
                        lprow = lprow + 1;
                        end
                        if comp.group>0
                            if swnum ~= comp.group
                                swnum = comp.group;
                                swrow = 1;
                                swBeginSta(swnum) = beginSta;
                                swEndSta(swnum) = endSta;
                            end
                            % EMA original:
%                             swBeginSta(swnum) = min(beginSta,swBeginSta(swnum));
%                             swEndSta(swnum) = max(endSta,swEndSta(swnum));
                            % changed to:
                            swBeginSta(swnum) = min([beginSta,swBeginSta(swnum)]);
                            swEndSta(swnum) = max([endSta,swEndSta(swnum)]);
                            % END
                            areas = obj.webareas{swnum}(beginSta(ks):endSta(ks)-1);
                            regionarea = sum(areas(:));
                            obj.bom(1).sw{swnum}{swrow, 1} = swrow;
                            obj.bom(1).sw{swnum}{swrow, 2} = comp.materialid;
                            obj.bom(1).sw{swnum}{swrow, 3} = comp.name;
                            obj.bom(1).sw{swnum}{swrow, 4} = obj.ispan(beginSta(ks));
                            obj.bom(1).sw{swnum}{swrow, 5} = obj.ispan(endSta(ks));
                            obj.bom(1).sw{swnum}{swrow, 6} = max(obj.webwidth{swnum});
                            obj.bom(1).sw{swnum}{swrow, 7} = mean(obj.webwidth{swnum});
                            obj.bom(1).sw{swnum}{swrow, 8} = regionarea;
                            obj.bom(1).sw{swnum}{swrow, 9} = layermult * mat.layerthickness;
                            obj.bom(1).sw{swnum}{swrow,10} = layermult * mat.drydensity * regionarea;
                            obj.bomIndices(1).sw{swnum}(swrow,1:2) = [beginSta(ks), endSta(ks)];
                            swrow = swrow + 1;
                        end
                    end
                end
            end
            obj.bom.lebond = sum(obj.LEbond)*m_to_mm;
            obj.bom.tebond = sum(obj.TEbond)*m_to_mm;
            obj.bom.dryweight = g_to_kg*(sum(cell2mat(obj.bom.hp(:,10))) ...
                                       + sum(cell2mat(obj.bom.lp(:,10))));
            for k=1:numel(obj.bom.sw)
                obj.bom.dryweight = obj.bom.dryweight + sum(cell2mat(obj.bom.sw{k}(:,10)));
                
                C = obj.webbonds{k}(:,swBeginSta(k):swEndSta(k)-1);
                obj.bom.swbonds{k} = m_to_mm*sum(C,2);
            end
            
            % build the material stack for each area
            nSegments = size(obj.keyareas,1);                   
            nStations = size(obj.keyareas,2);
            nWebs = numel(obj.bomIndices(1).sw);
            segmentLabels = {'HP_TE_FLAT','HP_TE_REINF','HP_TE_PANEL','HP_SPAR','HP_LE_PANEL','HP_LE',...
                             'LP_LE','LP_LE_PANEL','LP_SPAR','LP_TE_PANEL','LP_TE_REINF','LP_TE_FLAT'};
            obj.stacks = StackDef(); % clear old values
            obj.stacks(nSegments,nStations) = StackDef();  % allocate array
            for kr = 1:nSegments
                for kc = 1:nStations
                    % name the stacks <mm_span_location>_<segmentLabel>
                    obj.stacks(kr,kc).name = sprintf('%06d_%s',fix(m_to_mm*obj.ispan(kc)),segmentLabels{kr});
                    obj.stacks(kr,kc).indices = [kc, kc+1, kr, kr+1];
                end
            end
            for k = 1:size(obj.bom.hp,1)
                % for each row in the BOM, get the ply definition ...
                ply.component  = obj.bom.hp{k,3};   % parent component of ply
                ply.materialid = obj.bom.hp{k,2};   % materialid of ply
                ply.thickness  = obj.bom.hp{k,9};   % thickness [mm] of single ply
                ply.angle      = 0; % TO-DO, set to 0 for now, obj.bom.hp(k, );
                ply.nPlies     = 1; % default to 1, modified in addply() if necessary
                % ... and add the ply to every area that is part of the region
                ind = obj.bomIndices.hp(k,:);
                for kr = ind(3):ind(4)-1    
                    for kc = ind(1):ind(2)-1
                        obj.stacks(kr,kc) = obj.stacks(kr,kc).addply(ply);
                    end
                end
            end
            for k = 1:size(obj.bom.lp,1)
                % for each row in the BOM, get the ply definition ...
                ply.component  = obj.bom.lp{k,3};   % parent component of ply
                ply.materialid = obj.bom.lp{k,2};   % materialid of ply
                ply.thickness  = obj.bom.lp{k,9};   % thickness [mm] of single ply
                ply.angle      = 0; % TO-DO, set to 0 for now, obj.bom.hp(k, );
                ply.nPlies     = 1; % default to 1, modified in addply() if necessary
                % ... and add the ply to every area that is part of the region
                ind = obj.bomIndices.lp(k,:);
                for kr = ind(3):ind(4)-1    
                    for kc = ind(1):ind(2)-1
                        obj.stacks(kr,kc) = obj.stacks(kr,kc).addply(ply);
                    end
                end
            end
            
            obj.swstacks = cell(0); % clear old values
            for kw = 1:nWebs
                obj.swstacks{kw}(1,nStations) = StackDef(); % allocate array
                for kc = 1:nStations
                    % name the stacks <mm_span_location>_SW#
                    obj.swstacks{kw}(1,kc).name = sprintf('%06d_SW%d',fix(m_to_mm*obj.ispan(kc)),kw);
                    ind = obj.webindices{kw};  % currently, the shearweb indices do not change down the span
                    obj.swstacks{kw}(1,kc).indices = [kc,kc+1,ind(1),ind(2)];
                end
                for k = 1:size(obj.bom.sw{kw},1)
                    % for each row in the BOM, get the ply definition ...
                    ply.component  = obj.bom.sw{kw}{k,3};   % parent component of ply
                    ply.materialid = obj.bom.sw{kw}{k,2};   % materialid of ply
                    ply.thickness  = obj.bom.sw{kw}{k,9};   % thickness [mm] of single ply
                    ply.angle      = 0; % TO-DO, set to 0 for now, obj.bom.sw{kw}{k, };
                    ply.nPlies     = 1; % default to 1, modified in addply() if necessary
                    % ... and add the ply to every area that is part of the region
                    ind = obj.bomIndices.sw{kw}(k,:);
                    for kc = ind(1):ind(2)-1
                        obj.swstacks{kw}(1,kc) = obj.swstacks{kw}(1,kc).addply(ply);
                    end
                end
            end
            
            % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % need to add the 'MatDB' information which stores composite stack
            % information in each region at each station
            % prepare material database ==========================================
            obj.matdb = struct('type'     ,'','name'     ,'','reference','',...
                'dens'     ,[],'nuxy'     ,[],'ex'       ,[],'ey'       ,[],...
                'ez'       ,[],'gxy'      ,[],'gyz'      ,[],'gxz'      ,[],...
                'prxy'     ,[],'pryz'     ,[],'prxz'     ,[],'xten'     ,[],...
                'xcmp'     ,[],'yten'     ,[],'ycmp'     ,[],'zten'     ,[],...
                'zcmp'     ,[],'xy'       ,[],'yz'       ,[],'xz'       ,[],...
                'xycp'     ,[],'yzcp'     ,[],'xzcp'     ,[],'xzit'     ,[],...
                'xzic'     ,[],'yzit'     ,[],'yzic'     ,[],'g1g2'     ,[],...
                'etal'     ,[],'etat'     ,[],'alp0'     ,[],...
                'thicknessType',[],'uniqueLayers' ,[],'symmetryType' ,[],...
                'layer'        ,[]);
                        
            for k=1:length(obj.materials)
                obj.matdb(k).name = obj.materials(k).name;
                obj.matdb(k).type = obj.materials(k).type;
                obj.matdb(k).ex   = obj.materials(k).ex;
                obj.matdb(k).ey   = obj.materials(k).ey;
                obj.matdb(k).ez   = obj.materials(k).ez;
                obj.matdb(k).gxy  = obj.materials(k).gxy;
                obj.matdb(k).gyz  = obj.materials(k).gyz;
                obj.matdb(k).gxz  = obj.materials(k).gxz;
                if strcmpi(obj.matdb(k).type,'isotropic')
                    obj.matdb(k).nuxy = obj.materials(k).prxy;
                else
                    obj.matdb(k).prxy = obj.materials(k).prxy;
                    obj.matdb(k).pryz = obj.materials(k).pryz;
                    obj.matdb(k).prxz = obj.materials(k).prxz;
                end
                obj.matdb(k).dens = obj.materials(k).density;
                obj.matdb(k).reference = obj.materials(k).reference;
            end
            N = numel(obj.matdb);
            for k=1:numel(obj.stacks)
                obj.matdb(N+k).name          = obj.stacks(k).name;
                obj.matdb(N+k).type          = 'composite';
                obj.matdb(N+k).reference     = 'Reference text';
                obj.matdb(N+k).thicknessType = 'Constant';
                obj.matdb(N+k).uniqueLayers  = numel(obj.stacks(k).plygroups);
                obj.matdb(N+k).symmetryType  = 'none';
                for j=1:obj.matdb(N+k).uniqueLayers
                    matid = obj.stacks(k).plygroups(j).materialid;
                    obj.matdb(N+k).layer(j).layerName  = obj.matdb(matid).name;
                    obj.matdb(N+k).layer(j).thicknessA = mm_to_m * obj.stacks(k).plygroups(j).thickness;
                    obj.matdb(N+k).layer(j).thicknessB = obj.matdb(N+k).layer(j).thicknessA;
                    obj.matdb(N+k).layer(j).quantity   = obj.stacks(k).plygroups(j).nPlies;
                    obj.matdb(N+k).layer(j).theta      = obj.stacks(k).plygroups(j).angle;
                end
            end
            for kw=1:numel(obj.swstacks)
                N = numel(obj.matdb);
                for k=1:numel(obj.swstacks{kw})
                    obj.matdb(N+k).name          = obj.swstacks{kw}(k).name;
                    obj.matdb(N+k).type          = 'composite';
                    obj.matdb(N+k).reference     = 'Reference text';
                    obj.matdb(N+k).thicknessType = 'Constant';
                    obj.matdb(N+k).uniqueLayers  = numel(obj.swstacks{kw}(k).plygroups);
                    obj.matdb(N+k).symmetryType  = 'none';
                    for j=1:obj.matdb(N+k).uniqueLayers
                        matid = obj.swstacks{kw}(k).plygroups(j).materialid;
                        obj.matdb(N+k).layer(j).layerName  = obj.matdb(matid).name;
                        obj.matdb(N+k).layer(j).thicknessA = mm_to_m * obj.swstacks{kw}(k).plygroups(j).thickness;
                        obj.matdb(N+k).layer(j).thicknessB = obj.matdb(N+k).layer(j).thicknessA;
                        obj.matdb(N+k).layer(j).quantity   = obj.swstacks{kw}(k).plygroups(j).nPlies;
                        obj.matdb(N+k).layer(j).theta      = obj.swstacks{kw}(k).plygroups(j).angle;
                    end
                end
            end            
            % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     
            
            
            % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % shearweb information from NuMAD v1 is formatted in a specific
            % way, recreating that here
            % recreating data.shearweb ====================================
            ctr=1;
            for kw=1:numel(obj.swstacks)
                ind = obj.webindices{kw};
                for k=1:numel(obj.swstacks{kw})
                    if ~isempty(obj.swstacks{kw}(k).plygroups)
                        obj.shearweb(ctr).Material     = obj.swstacks{kw}(k).name;
                        obj.shearweb(ctr).BeginStation = obj.swstacks{kw}(k).indices(1);%=k
                        obj.shearweb(ctr).EndStation   = obj.swstacks{kw}(k).indices(2);%=k+1
                        obj.shearweb(ctr).Corner       = ind([2,1,1,2])-1; % dp number is offset by 1 in NuMAD v1
                        ctr=ctr+1;
                    end
                end
            end            
            % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>     
        end
        
        function readYAML(obj,file)
           % This method calls a separate matlab function, sends
           % the YAML file ('file') and blade object ('obj'), and returns the
           % blade object generated using the YAML file converter           
           obj = YAML_to_BladeDef(obj, file); 
        end
        
        function writeYAML(obj,file)
           % This method calls a separate matlab function to write the 
           % the Blade Object ('obj') to a YAML file ('file')           
           BladeDef_to_YAML(obj,file); 
        end
        
                
        function bmodesFrequencies = generateBeamModel(obj)
            % This method generates blade sectional properties used for
            % aeroelastic analyses
            
            global precompPath
            global bmodesPath
            
            parID = gcp('nocreate'); 
            if ~isempty(parID) % running in parallel, but not in a worker
                batchRun = true;
            else
                batchRun = false;
            end
                        
            % NOTES: ******************************************************
            % 1. FIXED -- needs to read MatDBsi.txt file, store this
            % internally (blade.matdb)
            % 2. FIX -- the code creates files prepmat -- why? is this
            % needed elsewhere?
            % 3. using blade.profiles instead of data.station(ii).coords
            % 4. FIX -- web data not saved in same format as NuMAD
            % (data.shearweb structure in NuMAD)
            % *************************************************************
            
            % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
            % code from BladeDef_to_NuMADfile and NuMAD2PreComp files
            % need to add material properties for end station
            if size(obj.stacks,2) < length(obj.ispan)
                for ii = 1:size(obj.stacks,1)
                    obj.stacks(ii,length(obj.ispan)).name = '**UNSPECIFIED**';
                end
            end
            % original code from NuMAD files ended ------------------------
            % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            
            % create precomp input files
            precomp = Blade2PreComp(obj,obj.matdb);
            % run precomp
            PreComp_SectionData = runPreCompAnalysis(precomp,precompPath,batchRun);
            
            % Use BModes to calculate mode shapes based on section properties from PreComp
            bmodesFrequencies = BModes4Blade2PreComp2FASTBlade(obj,bmodesPath,PreComp_SectionData.data);
            % fit polynomials to the modes
            if batchRun
                modeShapes=polyfitmodes([1 3 2]); % should change this to automatically assign mode order
            else
                modeShapes=polyfitmodes;
            end
                        
            % Generate FAST Blade file from this analysis
            warning('need PresweepRef and PrecurveRef variables from NuMAD')
            PreComp2FASTBlade_BladeDef(obj,PreComp_SectionData.data,modeShapes,batchRun);
            
            disp('FAST Blade file has been written: FASTBlade_precomp.dat')
            
            % Input file cleanup
            if batchRun
                qu='Yes'; % always cleanup during batch runs
            else
                qu=questdlg('Delete all miscellaneous input/output files?','File Cleanup...','Yes','No','No');
            end
            switch qu
                case 'Yes'
                    delete('*.pci')
                    delete('*.inp')
                    delete PrepMat.txt
                    delete bmodes*.bmi
                    delete blade_sec_props.dat
                    delete bmodes*.echo
                case 'No'
            end
            
            % Estimate the flexural bending distance
            warning('this code is not correct -- should be based on distance to neutral axis, not purely geometry')
            make_c_array_BladeDef(obj)
        end
        
        function generateFEA(obj) 
            % This method generates FEA    
            
            % 1. check that functionality from original code is not needed
            % 2. FIXED -- shell7 needs to read MatDBsi.txt file, store this
            % internally if the file isn't used elsewhere (blade.matdb)
            % 3. FIX -- in shell7 file, need to define PrecurveRef and
            % PresweepRef (currently set to zero)
            % 4. FIX -- web data not saved in same format as NuMAD
            % (data.shearweb structure in NuMAD)
            %             
            
            % NOTE:can add flags into the call -- e.g., element type, ...
            global ansysPath
            % define ANSYS model settings (can be options in generateFEA)
            config.ansys.BoundaryCondition = 'cantilered';
            config.ansys.ElementSystem = '181';
            config.ansys.MultipleLayerBehavior = 'distinct';
            config.ansys.meshing = 'elementsize';
            config.ansys.smartsize = 5;
            config.ansys.elementsize = 0.2026;
            config.ansys.shell7gen = 1;
            config.ansys.dbgen = 1;            
            fcopts = {'EMAX','SMAX','TWSI','TWSR','HFIB','HMAT','PFIB','PMAT',...
                'L3FB','L3MT','L4FB','L4MT','USR1','USR2','USR3','USR4',...
                'USR5','USR6','USR7','USR8','USR9'};
            config.ansys.FailureCriteria = cell(numel(fcopts),2);
            config.ansys.FailureCriteria(:,1) = fcopts';
            config.ansys.FailureCriteria(:,2) = deal({false});
            
                       
            
            
            % Generate a mesh using shell elements -- can add options here
            shell7_name = 'shell7.src'; ansys_product = 'ANSYS';
            obj.paths.job = pwd;% ble: is this needed? FIX THIS -- should update with parallel simulations??
            filename = fullfile(obj.paths.job,shell7_name);
                        
            develop__write_shell7(obj,filename);
            
            if obj.ansys.dbgen
                if isempty(ansysPath)
                    errordlg('Path to ANSYS not specified. Aborting.','Operation Not Permitted');
                    return;
                end                
                try
                    %tcl: exec "$ANSYS_path" -b -p $AnsysProductVariable -I shell7.src -o output.txt
                    ansys_call = sprintf('"%s" -b -p %s -I %s -o output.txt',...
                        ansysPath,ansys_product,shell7_name);
                    [status,result] = dos(ansys_call);  % the windows system call to run the above ansys command
                    
                    if isequal(status,0)
                        % dos command completed successfully; log written to output.txt
                        if 1%obj.batchrun
                            disp('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.');
                        else
                            helpdlg('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.','ANSYS Call Completed');
                        end
                    end
                    if isequal(status,7)
                        % an error has occured which is stored in output.txt
                        if 1%app.batchrun
                            disp('Could not complete ANSYS call. See "output.txt" for details.');
                        else
                            warndlg('Could not complete ANSYS call. See "output.txt" for details.','Error: ANSYS Call');
                        end
                    end
                catch ME
                    rethrow(ME);
                end
            end
        end
        
        function writeBOMxls(obj,file)
            % This method writes the bill-of-materials out to a spreadsheet.
            %            
            % Example:
            %            
            %   ``bladeDef.writeBOMxls('bom.xlsx')``
            
            m_to_mm = 1e3;            
            if exist('BOM_template.xlsx','file')
                copyfile('BOM_template.xlsx',file)
            end
            header = {'Layer #','Material ID','Component','Begin Station','End Station','Max width','Mean width','3D area','Layer Thickness','Computed layer weight';
                      ''       ,''           ,''         ,'(m)'          ,'(m)'        ,'(m)'      ,'(m)'       ,'(m^2)'  ,'(mm)'           ,'(g)'};
            % LP skin table      
            array = [header; obj.bom.lp];
            xlswrite(file,array,'LP skin');
            
            % HP skin table
            array = [header; obj.bom.hp];
            xlswrite(file,array,'HP skin');
            
            % shear web table
            array = [{'SW #';''}, header];  % insert another column
            for k=1:numel(obj.bom.sw)
                nr = size(obj.bom.sw{k},1);
                array = [array; [repmat({k},nr,1), obj.bom.sw{k}]]; %#ok<AGROW>
            end
            xlswrite(file,array,'shear webs');
            
            % bond line lengths
            array = {''       ,'Length';
                     ''       ,'(mm)';
                     'Root diameter',obj.ichord(1)*m_to_mm;
                     'LE bond',round(obj.bom.lebond);
                     'TE bond',round(obj.bom.tebond)};
            for r=1:2
                for k=1:numel(obj.bom.swbonds)
                    surfs = {'HP','LP'};
                    str = sprintf('%s bond, SW %d',surfs{r},k);
                    cellrow = {str,round(obj.bom.swbonds{k}(r))};
                    array = [array; cellrow]; %#ok<AGROW>
                end
            end
            xlswrite(file,array,'lengths');
            
            
        end
        
        function writePlot3D(obj,file,breakpoints)
            % Write the current blade geometry in Plot3D format.
            % breakpoints is a list of chord fractions at which the
            % surface geometry is divided into blocks
            %
            % Examples:
            %            
            %   ``BladeDef.writePlot3D(filename,[breakpoints])``
            %                        
            %   ``BladeDef.writePlot3D('file.p3d',[-.3, .3]);``
            %
            
            if ~exist('breakpoints','var')
                breakpoints = [];
            end
            indicesOfBreakpoints = zeros(1,numel(breakpoints));
            % get the chordwise spacing of points, assuming identical
            % spacing for all stations
            chordspacing = obj.cpos(:,1);
            for kBreakpoint = 1:numel(breakpoints)
                bp = breakpoints(kBreakpoint);
                [~,ind] = min(sqrt((chordspacing-bp).^2));
                indicesOfBreakpoints(kBreakpoint) = ind;
            end
            [N,M] = size(obj.cpos);
            INCLUDE_REPEATS = false;
            if INCLUDE_REPEATS
                indicesOfBreakpoints = unique([1, indicesOfBreakpoints, N]); %#ok<UNRCH>
            else
                indicesOfBreakpoints = unique([2, indicesOfBreakpoints, N-1]);
            end
            
            fid = fopen(file,'wt'); % open for Writing in Text mode
            if (fid == -1)
                error('Could not open file "%s"',file);
            end
            
            % output the data in Plot3d format
            try
                nBlocks = numel(indicesOfBreakpoints)-1; % number of xyz blocks
                fprintf(fid,'%d\n',nBlocks);
                for kblock = 1:nBlocks
                    a = indicesOfBreakpoints(kblock);
                    b = indicesOfBreakpoints(kblock+1);
                    fprintf(fid,'%d  %d  %d\n',1+b-a,M,1);
                end
                
                columnsPerLine = 5;
                for kblock = 1:nBlocks
                    a = indicesOfBreakpoints(kblock);
                    b = indicesOfBreakpoints(kblock+1);
                    obj.fprintf_matrix(fid,obj.geometry(a:b,1,:),columnsPerLine);
                    obj.fprintf_matrix(fid,obj.geometry(a:b,2,:),columnsPerLine);
                    obj.fprintf_matrix(fid,obj.geometry(a:b,3,:),columnsPerLine);
                end
                
            catch ME
                % The try..catch..end ensures the file gets closed
                % in case of a programming error.
                fclose(fid);
                rethrow(ME);
            end
            fclose(fid);
        end
        
        function tetype = getprofileTEtype(obj,k)
            if any(k<1 | k>size(obj.profiles,3))
                error('One of the requested indices is out of range.')
            end 
            tetype = cell(1,length(k));
            for i = 1:length(k)
                xy = obj.profiles(:,:,k(i));
                tetype{i} = obj.getTEtype(xy);
            end
        end
        
        function coords = downsampleProfile(obj,k,n_points)
            % ble: can this be deleted?? Doesn't precisely control the
            % number of points. replace with resampleAirfoil_ble
            % currently used by BladeDef_to_NuMADfile and develop__write_shell7
            N = size(obj.profiles,1);
            assert(isscalar(k),'Profile index "k" must be scalar');
            assert(n_points<N,'n_points must be less than length of profile');
            dk = round((N-2)/n_points);
            LE = obj.LEindex;
            tetype = obj.getprofileTEtype(k);
            switch tetype{1}
                case 'flat'
                    ind = unique([1, 2:dk:LE, LE, fliplr(N-1:-dk:LE)]);
                case {'sharp','round'}
                    ind = unique([   2:dk:LE, LE, fliplr(N-2:-dk:LE)]);
            end
            coords = obj.profiles(ind,:,k);
        end
        
        function delete(obj)
            try %#ok<TRYNC>
                delete(obj.bomPlot.hgLinesHP);
                delete(obj.bomPlot.hgLinesLP);
                delete(obj.bomPlot.hgPatchHP);
                delete(obj.bomPlot.hgPatchLP);
                delete(obj.bomPlot.uisliderHP);
%                 delete(obj.bomPlot.uisliderLP);
                delete(obj.bomPlot.hTitleHP);
                delete(obj.bomPlot.hTitleLP);
            end
        end
        
        function varargout=surf(obj)
            h=surf(squeeze(obj.geometry(:,3,:)),...
                 squeeze(obj.geometry(:,1,:)),...
                 squeeze(obj.geometry(:,2,:)),...
                 'MeshStyle','column');
            if nargout > 0
                varargout = {h};
            end
        end
        
        function plotregions(obj)
            try %#ok<TRYNC>
                delete(obj.hgKeypoints)
            end
            M = size(obj.keypoints,1);
            for km = 1:M
                obj.hgKeypoints(km) = line(...
                    squeeze(obj.keypoints(km,3,:)), ...
                    squeeze(obj.keypoints(km,1,:)), ...
                    squeeze(obj.keypoints(km,2,:)) );
            end
        end
        
        function plotgeom(obj)
            try %#ok<TRYNC>
                delete(obj.hgGeometry)
            end
            N = size(obj.geometry,3);
            for k = 1:N
                obj.hgGeometry(k) = line(...
                    obj.geometry(:,3,k),...
                    obj.geometry(:,1,k),...
                    obj.geometry(:,2,k));
            end
        end
        
        function plotbom(obj,k)
            
            if ~exist('k','var') || isempty(k)
                k = 1;
            end
            if isequal(k,'cb')
                k = fix(get(obj.bomPlot.uisliderHP,'Value'));
            end
            
            if isempty(obj.bomPlot.hgLinesHP) || ~all(ishandle(obj.bomPlot.hgLinesHP))
            clf;
            obj.bomPlot.axHP = axes('Position',[0.1 0.6 0.8 .3]);
            obj.bomPlot.hgLinesHP = ...
                plot(obj.ispan,obj.keyarcs(1,:)-obj.HParcx0,'k-.',...
                     obj.ispan,obj.keyarcs(7,:)-obj.HParcx0,'k-.');
            obj.bomPlot.hTitleHP = title('','Interpreter','none');
            
            obj.bomPlot.axLP = axes('Position',[0.1 0.2 0.8 .3]);
            obj.bomPlot.hgLinesLP = ...
                plot(obj.ispan,obj.keyarcs(13,:)-obj.LParcx0,'k-.',...
                     obj.ispan,obj.keyarcs( 7,:)-obj.LParcx0,'k-.');
            obj.bomPlot.hTitleLP = title('','Interpreter','none');
            
            n = size(obj.bom.hp,1);
            obj.bomPlot.uisliderHP = uicontrol('Style','slider',...
                'Min',1,'Max',n,'Value',1,'SliderStep',[1/(n-1) 10/(n-1)],...
                'Units','normalized','Position',[.02 .02 .96 .04],...
                'Callback',@(src,evt) obj.plotbom('cb'));
%             obj.bomPlot.uisliderLP = uicontrol('Style','slider',...
%                 'Min',1,'Max',n,'Value',1,'SliderStep',[1/(n-1) 10/(n-1)],...
%                 'Units','normalized','Position',[.02 .02 .96 .04],...
%                 'Callback',@(src,evt) obj.plotbom('cb'));
            end
            
            k = min(k,size(obj.bom.hp,1));
            obj.bomPlot.kLayer = k;
            set(obj.bomPlot.uisliderHP,'Value',k);
            str = sprintf('HP Layer %d: %s',k,obj.bom.hp{k,3});
            set(obj.bomPlot.hTitleHP,'String',str);
            str = sprintf('LP Layer %d: %s',k,obj.bom.lp{k,3});
            set(obj.bomPlot.hTitleLP,'String',str);
            
            hp = obj.bomIndices.hp(k,:);
            x = obj.ispan(hp(1):hp(2));
            y1 = obj.keyarcs(hp(3),hp(1):hp(2))-obj.HParcx0(hp(1):hp(2));
            y2 = obj.keyarcs(hp(4),hp(1):hp(2))-obj.HParcx0(hp(1):hp(2));
            if ishandle(obj.bomPlot.hgPatchHP)
                delete(obj.bomPlot.hgPatchHP)
            end
            axes(obj.bomPlot.axHP);
            obj.bomPlot.hgPatchHP = patch([x,fliplr(x)], [y1,fliplr(y2)], 'b');
            
            lp = obj.bomIndices.lp(k,:);
            x = obj.ispan(lp(1):lp(2));
            y1 = obj.keyarcs(lp(3),lp(1):lp(2))-obj.LParcx0(lp(1):lp(2));
            y2 = obj.keyarcs(lp(4),lp(1):lp(2))-obj.LParcx0(lp(1):lp(2));
            if ishandle(obj.bomPlot.hgPatchLP)
                delete(obj.bomPlot.hgPatchLP)
            end
            axes(obj.bomPlot.axLP);
            obj.bomPlot.hgPatchLP = patch([x,fliplr(x)], [y1,fliplr(y2)], 'b');
        end
        
        function plotprofile(obj,k)
            % This method plots profiles
            %            
            % Examples:
            %            
            %   ``blade.plotprofile(1);``
            %
            %   ``blade.plotprofile(1:N);``
            %                        
            plot(squeeze(obj.profiles(:,1,k)),...
                 squeeze(obj.profiles(:,2,k)),'.-')
        end
            
        function [beginSta,endSta] = findLayerExtents(obj,layerDist,layerN)
            % This method... 

            assert(isscalar(layerN),'second argument ''layerN'' must be a scalar');
            staLogical = layerDist >= layerN;
            prev = 0;
            beginSta = []; endSta = [];
            for k=1:length(staLogical)
                if staLogical(k)==1 && prev==0
                    beginSta(end+1) = k; %#ok<AGROW>
                end
                if staLogical(k)==0 && prev==1
                    endSta(end+1) = k; %#ok<AGROW>

                elseif k==length(staLogical) && prev==1
                    endSta(end+1) = k; %#ok<AGROW>
                end
                prev = staLogical(k);

            end
        end


        function [hpRegion, lpRegion, swRegion] = findRegionExtents(obj,keylabels,comp)
            % This method...

            le = find(1==strcmpi('le',keylabels));
            % "keylabels" is expected to wrap from te on hp side around to te on lp side
            if length(comp.hpextents)==2
                hp1 = find(1==strcmpi(comp.hpextents{1},keylabels(1:le)));
                assert(~isempty(hp1),'HP extent label "%s" not defined.',comp.hpextents{1});
                hp2 = find(1==strcmpi(comp.hpextents{2},keylabels(1:le)));
                assert(~isempty(hp2),'HP extent label "%s" not defined.',comp.hpextents{2});
                hpRegion = sort([hp1 hp2]);
            else
                hpRegion = [];
            end
            if length(comp.lpextents)==2
                lp1 = find(1==strcmpi(comp.lpextents{1},keylabels(le:end))) + le-1;
                assert(~isempty(lp1),'LP extent label "%s" not defined.',comp.lpextents{1});
                lp2 = find(1==strcmpi(comp.lpextents{2},keylabels(le:end))) + le-1;
                assert(~isempty(lp2),'LP extent label "%s" not defined.',comp.lpextents{2});
                lpRegion = sort([lp1 lp2]);
            else
                lpRegion = [];
            end
        %     if length(comp.hpextents)==1 && length(comp.lpextents)==1
        %         sw1 = find(1==strcmpi(comp.hpextents{1},keylabels(1:le)));
        %         assert(~isempty(sw1),'HP extent label "%s" not defined.',comp.hpextents{1});
        %         sw2 = find(1==strcmpi(comp.lpextents{1},keylabels(le:end))) + le-1;
        %         assert(~isempty(sw2),'LP extent label "%s" not defined.',comp.lpextents{1});
        %         swRegion = [sw1 sw2];
        %     else
                swRegion = [];
        %     end
        end


        function tetype = getTEtype(obj,xy)
            % This method...

            if abs(xy(2,2)-xy(end-1,2)) > 1e-5
                % y-diff of second and end-1 points is non-zero for flatback
                tetype = 'flat';
                disp('FLATBACK AIRFOIL')
            else
                % y-diff of first two points is zero otherwise (point
                % is duplicated)
                hp_angle = atan2(xy(2,2)-xy(3,2),xy(2,1)-xy(3,1));
                lp_angle = atan2(xy(end-2,2)-xy(end-1,2),xy(end-1,1)-xy(end-2,1));
                if (hp_angle + lp_angle) > 0.8*pi
                    % if angle is approaching 180deg, then treat as
                    % 'round'
                    % jcb: it may be better to base this decision on
                    % continuity of slope or curvature
                    tetype = 'round';
                else
                    tetype = 'sharp';
                end
            end
        end

        function fprintf_matrix(obj,fid,matrixData,columnsPerLine)
            % This method...

            kColumn = 1;
            for kData=1:numel(matrixData)
                fprintf(fid,'%g ',matrixData(kData));
                kColumn = kColumn + 1;
                if kColumn > columnsPerLine
                    fprintf(fid,'\n');
                    kColumn = 1;
                end
            end
        end
    end
end
