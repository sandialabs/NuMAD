%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef AirfoilDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``AirfoilDef``  A class definition for airfoil profiles.
% 
% Examples: 
% 
%	``af = AirfoilDef();``
% 
%	``af = AirfoilDef(FILENAME);`` where FILENAME is the file containing airfoil profile data in various
%   2-column formats including NuMAD-xml
%
% See also ``AirfoilDef.resample``, ``AirfoilDef.writeAirfoil``,
% ``AirfoilDef.convertAirfoil``, ``AirfoilDef.plot``
% ``xlsBlade``, ``BladeDef``, ``StationDef``
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    properties
        name   = ''     % String: User selected airfoil name
        reference       % Header info in file
        coordinates     % Profile data in two columns
        c               % Computed by NuAMD
        camber          % Camber line as a function of x. Distance in percent chord between LP and HP curves. Computed by NuAMD.
        thickness       % Float: Relative thickness as a function of the x coordinate. Values between 0 and 1, where 1 corresponds to maximum thickness. Computed by NuAMD.
        percentthick    % Float: Maximum airfoil thickness as a percentage of chord length [%]
        maxthick        % Float: Airfoil thickness as a percentage of chord. Values between 0 and 1.
		TEtype  = ''	% String: Options, ``‘round'``, ``'sharp'``, or ``'flat'``
    end
    
    properties (Dependent)
        x               % Horizontal axis of Airfoil shape coordinates Working clockwise starting from the TE to the LE and back to the TE. LE must be at (1,0) and TE at (0,0). Needed only by ``AirfoilDef.plot``
        y               % Vertical axis of Airfoil shape coordinates Working clockwise starting from the TE to the LE and back to the TE. LE must be at (1,0) and TE at (0,0). Needed only by ``AirfoilDef.plot``
    end
    
    methods
        function obj = AirfoilDef(filename)
            if ~exist('filename','var') || isempty(filename)
                obj.name = 'circular';
                obj.reference = '';
				theta = linspace(0,pi,50);
                theta = [theta (pi+theta(2:end-1))];
                xcoord = 0.5*cos(-theta) + 0.5;
                ycoord = 0.5*sin(-theta);
                obj.coordinates = [xcoord(:) ycoord(:)];

            else
                % Open the file and read the entire contents
                fid = fopen(filename);
                assert(fid ~= -1,'FileIO:FileOpenError',...
                    'Could not open file "%s"',filename);
                filecontents = fread(fid,inf,'uint8=>char')';
                fclose(fid);
                
                % Use the base filename as the airfoil name
                [~,fn,~] = fileparts(filename);
                obj.name = fn;
                
                % Read the file contents
                try
                    if ~isempty(regexp(filecontents,'<coords>','once'))
                        % files is numad's xml-style
                        [coords, ref] = readAirfoilXML(filecontents);
                        obj.reference = ref;
                        obj.coordinates = coords;
                    else
                        [coords, ref] = readAirfoilColumns(filecontents);
                        obj.reference = ref;
                        obj.coordinates = coords;
                    end
                catch ME
                    fprintf('Error interpreting airfoil file %s\n',filename);
                    rethrow(ME);
                end
            end
%             obj.resample(400);
        end
        
        function resample(obj,n_samples,spacing)
            % AirfoilDef.resample
            % af.resample(n_samples,spacing)
            %   spacing = 'auto' | 'half-cosine' | 'cosine'
            %
            % af.resample(200,'half-cosine');
            if ~exist('n_samples','var') || isempty(n_samples)
                n_samples = 400;
            end
            if ~exist('spacing','var') || isempty(spacing)
                spacing = 'auto';
            end
            for k=1:numel(obj)
                af_out = resampleAirfoil(obj(k).coordinates,n_samples,spacing);
                xcoord = af_out(:,1);
                ycoord = af_out(:,2);
                % obj(k).percentthick = (max(ycoord) - min(ycoord))*100;
                [obj(k).c, obj(k).camber, obj(k).thickness] = computeCamberAndThickness(xcoord,ycoord);
                [m,i] = max(obj(k).thickness);
                obj(k).percentthick = m*100;
                obj(k).maxthick = obj(k).c(i);
                obj(k).TEtype = getTEtype(af_out);
            end
        end
        
        function adjustTE(obj,tet,tes,onset)
            % AirfoilDef.adjustTE
            % af.adjustTE(TE_thick,[TE_slope],[onset])
            %   TE_thick is the amount of TE thickness to add
            %   TE_slope is the slope of the added thickness profile at TE,
            %            defaults to 5/3 * TE_thick
            %   onset    is the chord fraction where adjustment begins,
            %            defaults to location of max thickness
            %   
            % af.adjustTE(0.02);
            % af.adjustTE(0.02,0);
            % af.adjustTE(0.02,[],0.8);
            if ~exist('tes','var') || isempty(tes)
                tes = 5/3 * tet;       % slope of TE adjustment; 5/3*tet is "natural"
            end
            if ~exist('onset','var') || isempty(onset)
                USEMAXTHICK = true;
            else
                USEMAXTHICK = false;  % use the given 'onset' instead
            end
            % continuous first & second derivatives at 'onset'
            % maintain second & third derivative at mc==1 (TE)
            % adjust slope at mc==1 (TE) by tes
            A = [1,  1,  1,   1;
                 3,  4,  5,   6;
                 6, 12, 20,  30;
                 6, 24, 60, 120];
            d = [tet; tes; 0; 0];
            p = A\d;
            for k=1:numel(obj)
                if USEMAXTHICK
                    onset = obj(k).maxthick;  % start of TE modification, measured from LE
                end
                mc = max((obj(k).c - onset)/(1-onset),0);  % coordinate for TE mod
                temod = [mc.^3, mc.^4, mc.^5, mc.^6] * p;
                obj(k).thickness = obj(k).thickness + temod;
            end
        end
        
        function xcoord = get.x(obj)
            cc = obj.c;
            xcoord = [cc(end); flipud(cc); cc(2:end); cc(end)];
        end
        
        function ycoord = get.y(obj)
            lp = obj.camber+obj.thickness/2;
            hp = obj.camber-obj.thickness/2;
            ycoord = [0; flipud(hp); lp(2:end); 0];
        end
        
        function writeAirfoil(obj, filename, fileformat)
            % AirfoilDef.writeAirfoil
            % af.writeAirfoil(FILENAME,fileformat)
            %  fileformat = 'xfoil' | 'rfoil' | 'selig'
            %             = 'xfoil-reverse' | 'rfoil-reverse'
            %             = 'numad'
            %             = 'lednicer'
            %
            % af.writeAirfoil('af-test.txt','xfoil');
            fid = fopen(filename,'wt');
            assert(fid ~= -1,'FileIO:FileOpenError',...
                'Could not open file "%s"',filename);
            % determine precision of data
            tmp = obj.coordinates(:);
            for k=1:14  % maximum 14 digits of precision
                prec = k;
                if any(rem(tmp,10^(-k))) == 0
                    break
                end
            end
            % create the format spec based on the precision
            % e.g. '%10.7f %10.7f\n'
            formatspec = sprintf('%%%d.%df %%%d.%df\\n',prec+3,prec,prec+3,prec);
            % output the specified file format
            try
                switch lower(fileformat)
                    case {'xfoil','rfoil','selig'}
                        % write header lines on one line
                        ref = regexprep(obj.reference,'[\n\r]*',' | ');
                        fprintf(fid,'%s\n',ref);
                        % write data
                        coords = flipud(obj.coordinates);
                        fprintf(fid,formatspec,coords');
                    case {'xfoil-reverse','rfoil-reverse'}
                        % write header lines on one line
                        ref = regexprep(obj.reference,'[\n\r]*',' | ');
                        fprintf(fid,'%s\n',ref);
                        % write data
                        fprintf(fid,formatspec,obj.coordinates');    
                    case {'numad'}
                        % write header lines on separate lines
                        fprintf(fid,'<reference>\n');
                        fprintf(fid,'%s\n',obj.reference);
                        fprintf(fid,'</reference>\n');
                        fprintf(fid,'<coords>\n');
                        coords = obj.coordinates;
                        if isequal(coords(1,:),coords(end,:))
                            coords(end,:) = [];
                        end
                        fprintf(fid,formatspec,coords');
                        fprintf(fid,'</coords>');
                    case {'lednicer'}
                        coords = obj.coordinates;
                        dc = diff(coords(:,1));
                        dc = dc * sign(dc(1));
                        k = find(dc<0,1);
                        HPside = coords(1:k,:);
                        LPside = coords(k:end,:);
                        % write header lines on one line
                        ref = regexprep(obj.reference,'[\n\r]*',' | ');
                        fprintf(fid,'%s\n',ref);
                        % write table sizes
                        tmp = sprintf('%%%dd %%%dd\\n',prec+3,prec+3);  % e.g. '%10d %10d\n'
                        fprintf(fid,tmp,size(LPside,1),size(HPside,1));
                        fprintf(fid,'\n');
                        % write LP data, LE to TE
                        fprintf(fid,formatspec,LPside');
                        fprintf(fid,'\n');
                        % write HP data, LE to TE
                        fprintf(fid,formatspec,flipud(HPside)');
                    otherwise
                        error('Unknown format: %s',fileformat)
                end
                fclose(fid);
            catch ME
                fclose(fid);
                rethrow(ME);
            end
            
        end
        
        function other = copyobj(obj)
            other      = AirfoilDef;
            prop_array = properties(obj);
            for k=1:numel(prop_array)
                prop = prop_array{k};
                other.(prop) = obj.(prop);
            end
        end
        
        function plot(obj)
            % AirfoilDef.plot
            % plot(af)  Plots an AirfoilDef object (or an array of them)
            N = numel(obj);
            if N==1
                nc=1;
                nr=1;
            else
                nc = fix(sqrt(N));
                nr = ceil(N/nc);
            end
            clf;
            for k=1:N
                subplot(nr,nc,k);
                if ~isempty(obj(k).c)
                    plot(obj(k).x,obj(k).y,'.-');
                    line(obj(k).c,obj(k).camber,'LineStyle','-.','Color','k');
                    mtx = obj(k).maxthick*[1 1];
                    kn = find(obj(k).c>=obj(k).maxthick,1);
                    mty = obj(k).camber(kn) + obj(k).thickness(kn)*[.5 -.5];
                    line(mtx,mty,'LineStyle',':','Color','k');
                else
                    plot(obj(k).coordinates(:,1),obj(k).coordinates(:,2),'.-');
                end
                axis equal;
                title(obj(k).name,'Interpreter','none');
            end
        end
    end
    
    methods (Static)
        function convertAirfoil(infile,outfile,fileformat)
            % AirfoilDef.convertAirfoil(INFILE,OUTFILE,fileformat)
            %  fileformat = 'xfoil' | 'rfoil' | 'selig'
            %             = 'xfoil-reverse' | 'rfoil-reverse'
            %             = 'numad'
            %             = 'lednicer'
            %
            % AirfoilDef.convertAirfoil('airfoils/e205.dat','af-test.txt','xfoil')
            if isequal(infile,outfile)
                error('Aborting: writing over input file not allowed')
            end
            af = AirfoilDef(infile);
            af.writeAirfoil(outfile,fileformat);
        end
    end
    
end % classdef

%==========================================================================
function [coords, reference] = readAirfoilXML(filecontents)
% [coords, reference] = readAirfoilXML(filecontents)

% The following regular expression pattern matches any number of characters
% found between the opening and closing "reference" tags
pattern = '<reference>(.*)</reference>';
t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
try
    % jcb: is there a better way to extract the contents of t?
    reference = cell2mat(t{1});
    reference = regexprep(reference,'[\n\r]*','\n');
    reference = strtrim(reference);
catch ME  %#ok
    reference = '';
end

% The following regular expression pattern matches any number of characters
% found between the opening and closing "coords" tags
pattern = '<coords>(.*)</coords>';
t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
try
    % jcb: is there a better way to extract the contents of t?
    coord_text = cell2mat(t{1});
catch ME
    error('Airfoil <coords>..</coords> not found in file "%s".  Please check the file format.',filename);
end
% Convert coord_text to floating point; coords is an Nx2 matrix
% coords(:,1) are the x coordinates and coords(:,2) are the y coordinates
coords = cell2mat(textscan(coord_text,'%f %f'));
end

function [coords, reference] = readAirfoilColumns(filecontents)
% [coords, reference] = readAirfoilColumns(filecontents)

    % All of these file formats assume that the
    % LE is at (0,0) and the TE is at (1,0)
    raw = regexp(filecontents,'[^\n\r]*','match');  % get lines
    Nraw = numel(raw);
    kh = 1;  % index counter for header lines
    kt = 1;  % index counter for tables
    kr = 1;  % index counter for rows
    header = cell(0);
    table = cell(0);
    for k=1:Nraw
        % try to read pairs of coordinates
        pair = cell2mat(textscan(raw{k},'%f %f'));
        if isempty(pair)  % if empty...
            if kt>1 || kr>1  % ...and we already have some data
                % then move to a new table
                kt = kt + 1;
                kr = 1;
            else
                % otherwise keep reading the header
                header(kh) = raw(k);
                kh=kh+1;
            end
        else
            % place coordinate pair in table
            table{kt}(kr,:) = pair;
            kr = kr + 1;
        end
    end
    if numel(table) == 1
        % assume points wrap around either LE or TE
        coords = table{1};
        if coords(1,1) < 0 || coords(1,1) > 1
            error('First x-coordinate not in range 0..1')
        end
        dc = diff(coords(:,1));
        dc = dc * sign(dc(1));
        k = find(dc<0,1);
        sideA = coords(1:k,2);
        sideB = coords(k:end,2);
        if mean(sideA) > mean(sideB)
            % LP (upper) surface given first, so flipud
%             disp('LP first');
            coords = flipud(coords);
            k = size(coords,1) - k + 1;
        end
        if (1 - coords(1,1)) > 0.5
            % coordinates begin at LE and wrap around TE
%             disp('TE wrap');
            if isequal(coords(1,:),coords(end,:))
                coords(end,:) = [];
            end
            coords = [coords(k:-1:1,:); coords(end:-1:k+1,:)];
        end
    elseif numel(table) == 3
        % assume "Lednicer's" format 
        %(see http://www.ae.illinois.edu/m-selig/ads.html)
        npoints = table{1};   % one line gives number of points in each table
        lp = table{2}; % first table is upper surface LE to TE
        hp = table{3}; % second table is lower surface LE to TE
        if size(npoints,1) ~= 1 
            error('Format similar to "Lednicers", but more than one row found for table sizes')
        end
        if isequal(hp(1,:),lp(1,:))
            lp(1,:) = [];
        end
        coords = [flipud(hp); lp];
    else
        error('File format not recognized')
    end
    if numel(header)>=1
        reference = header{1};
        for k=2:numel(header)
           reference = sprintf('%s\n%s',reference,header{k});
        end
    else
        reference = '';
    end
end

%==========================================================================
function af_out = resampleAirfoil(af_in, n_samples, spacing)
%RESAMPLEAIRFOIL_NMD  Resample airfoil coordinates
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% af_out = resampleAirfoil_nmd(af_in, n_samples, spacing)
%
%       af_in = Nx2 array containing N normalized xy airfoil points
%    n_samples = number of points to be created around surface
%     spacing = spacing routine to be used
%                 'cosine', 'half-cosine', 'constant', 'auto'
%
%      af_out = array containing n_samples+1 airfoil points
%
% Cosine spacing: puts higher density of points at both LE and TE;
%       constant arc length point spacing around a perfect circle.
% Half-cosine spacing: puts higher density of points at LE and lesser
%       density of points at TE
% constant spacing: constant spacing of points along chord line
% auto: choose between Cosine and Half-cosine based on steepness of TE
%
% Assumes coordinates begin at trailing edge (x,y)=(1,0) and trace out the HP
% surface, then the LP surface.
%
% Flatback airfoil inputs are designated by ensuring the following:
%    Point #1 = (1,0)   (center of trailing edge)
%    Point #2 = (1,-y) where y~=0  (HP corner of flatback trailing edge)
%    Point #N = (1,y)  where y~=0  (LP corner of flatback trailing edge)
%
% Notes:
%  * This routine enforces LE point is at (0,0) - Warning, this may have
%    complications when it comes ot CFD analyses!!
%  * spline_type = spline algorithm to be used for oversampling:
%       'linear', 'pchip', 'spline'; right now, the function is hard-coded
%       for 'spline', but others can be used by changing the te_type setting
%       in the code.
%
% Assumes leading edge at (0,0) and trailing edge at (1,0)
%
%

% JP: initial creation
% BRR: modified for release 1/25/2011

% Error checking
% if airfoil coordinates are not in Nx2 array
if size(af_in,2)~=2
    tmpN=size(af_in,1);
    tmpM=size(af_in,2);
    warning(['af_in array was defined in ' num2str(tmpN) 'x' num2str(tmpM) ' array.  Automatically changing it to an ' num2str(tmpM) 'x' num2str(tmpN) ' array.'])
    af_in=af_in';
end
% End error checking routines
xy=af_in;
% ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% add a (1,0) point for the TE if it does not exist; important for flatback
% airfoil calculations
if xy(1,1)==1 && xy(1,2)~=0
    xy = [1 0; xy];
elseif xy(1,1)~=1
    error('airfoil coordinates should include x=1 values for the TE.')
end
% ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
xyTE=xy(1,:); % need to store trailing edge point for flatbacks
% determine trailing edge type
% flatback: x-coord of 1,2,last point are equal, y-coord not equal
if (xy(2,1)==xy(1,1)) && (xy(2,2)~=xy(1,2)) ...
        && (xy(end,1)==xy(1,1)) && (xy(end,2)~=xy(1,2))
    te_type = 'flatback';
    xy(1,:) = [];  % do not include first point in spline
else
    te_type = 'normal';
    % make sure the airfoil is closed (so spline extends to TE)
    if (xy(end,1)~=xy(1,1))
        xy(end+1,:) = xy(1,:);
    end
end

%Calculate arc length of xy points clockwise from trailing edge
n_points=size(xy,1);
t=zeros(1,n_points);
for i=2:n_points
    %formula: t(i) = hypot( x(i)-x(i-1) , y(i)-y(i-1) ) + t(i-1);
    t(i) = hypot(xy(i,1)-xy(i-1,1),xy(i,2)-xy(i-1,2)) + t(i-1);
end

%Get total arc length
arc_length=t(end);

%Spline airfoil with many points
oversample=10e3;
delta=arc_length/(oversample-1);
%  The manypoints range from 0 to total arc_length, adding a bit on each 
%  side so that flatbacks extend past x=1 after rotation corrections.
manypoints = -delta:delta:arc_length+delta;
spline_type='spline';
switch spline_type
    case {'linear','pchip','spline'}
        xxyy = interp1(t,xy,manypoints,spline_type);
    otherwise
        fprintf('Airfoil oversampling algorithm specified is not an available option. Defaulting to "spline".');
        xxyy = interp1(t,xy,manypoints,'spline');
end

% Normalize the airfoil: 
%   correct rotation so that LE is at (0,0) and TE is at (1,0).
% jcb: Technially, the LE is at the point of max curvature, but that
% definition can produce situations that break the interpolation step.
% Instead, we define the LE as the point that is the furthest distance from
% the TE.
xxyy = xxyy - repmat(xyTE,size(xxyy,1),1);
rays = hypot( xxyy(:,1), xxyy(:,2) ); % distance of each point from the TE
[max_ray, max_point] = max(rays);
ray_angle = atan2( xxyy(max_point,2), -xxyy(max_point,1) );
xxyy = rotate2d(xxyy,ray_angle);
xxyy = xxyy/max_ray + repmat([1,0],size(xxyy,1),1);

%Separate into high and low pressure surfaces
HP=xxyy(1:max_point,:);  % HP points progress from TE (x=1) to LE (x=0)
LP=xxyy(max_point:end,:);  % LP points progress from LE (x=0) to TE (x=1)

% if 'auto', determine which spacing algorithm to use
if strcmp('auto',spacing)
    dx = xxyy(2,1)-xxyy(3,1);
    % If x-spacing of the oversampled data at the TE is below a threshold, 
    % assume that cosine spacing would be best choice for the profile.
    if dx < 1/(10*oversample)  
        spacing = 'cosine';
    else
        spacing = 'half-cosine';
    end
end
%Calculate x points based on spacing algorithm specified
n_panels = fix(n_samples/2) - 1;
switch spacing
    case 'cosine'
        beta=linspace(0,pi,n_panels+1);
        x_fwd=0.5*(1-cos(beta));
    case 'half-cosine'
        beta=linspace(0,pi/2,n_panels+1);
        x_fwd=(1-cos(beta));
    case 'constant'
        x_fwd=linspace(0,1,n_panels+1);
    otherwise
        error('Resampling algorithm specified is not an available option');
end
x_fwd = x_fwd(:);       % make x_fwd a column vector
x_rev = flipud(x_fwd);  % x_rev values are nominally 1 to 0

%Calculate interpolated airfoil points. For sharp trailing edge airfoils,
%the trailing edge point is not repeated
LP_new=[x_fwd  interp1(LP(:,1),LP(:,2),x_fwd)];
HP_new=[x_rev  interp1(HP(:,1),HP(:,2),x_rev)];
% Make sure that LE point is at (0,0)
HP_new(end,:)=[0 0];
xyTE = [1 0];

%Assemble the two curves into a continuous line
af_out=[xyTE; HP_new(:,:); LP_new(2:end,:); xyTE];

%Get index of leading edge point in resampled coordinates
% LE = 1 + (n_panels + 1);  % or 1 + length(HP_new)

if 0  % debugging plots
    figure(701)
    subplot(2,1,1)
    plot(t,xy(:,1),'k-x',t,xy(:,2),'r-x',manypoints,xxyy(:,1),'b-',manypoints,xxyy(:,2),'g-')
    xlabel('airfoil surface distance from first point')
    ylabel('oversampled coordinate value')
    legend('af\_in x','af\_in y','oversampled x','oversampled y','Location','North')
    subplot(2,1,2)
    plot(xy(:,1),xy(:,2),'r-x',xxyy(:,1),xxyy(:,2),'b',af_out(:,1),af_out(:,2),'k-o')
    axis equal
    legend('af\_in points',['af\_in points oversampled by factor of ' num2str(oversample)],'af\_out','Location','SouthEast')
end

end

function xyout = rotate2d(xyin,angle)
    xyout(:,1) = cos(angle) * xyin(:,1) - sin(angle) * xyin(:,2);
    xyout(:,2) = sin(angle) * xyin(:,1) + cos(angle) * xyin(:,2);
end

function [c, camber, thickness] = computeCamberAndThickness(x,y)
    n_samples = length(x);
    LE = fix(n_samples/2) + 1;
    xhp = x(LE:-1:2);
    xlp = x(LE:end-1);
    assert(length(xhp)==length(xlp),'Error computing camber and thickness.');
    assert(sum(xhp-xlp)<eps,'Upper and lower surface x-coordinates must align.');
    yhp = y(LE:-1:2);
    ylp = y(LE:end-1);
    
    c = xlp;
    camber = (yhp + ylp)/2;
    thickness = abs(ylp - yhp);
end

function tetype = getTEtype(xy)
    if abs(xy(1,2)-xy(2,2)) > eps(1)
        % y-diff of first two points is non-zero for flatback
        tetype = 'flat';
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