function afdb = resampleAirfoilDB(afdb, n_panels, spacing)
%RESAMPLEAIRFOILDB  Resample airfoil coordinates of entire database
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% afdb = resampleAirfoilDB(afdb, n_panels, spacing)
%
%        afdb = airfoil database structure
%    n_panels = number of panels to be created per surface
%     spacing = spacing routine to be used
%                 'cosine', 'half-cosine', 'constant'
%
% See also resampleAirfoil_nmd

n_airfoils = numel(afdb);  % number of airfoils to loop through

% determine number of points per airfoil
max_points=-inf;
for kaf=1:n_airfoils
    max_points = max(size(afdb(kaf).coords,1),max_points);
end
if 2*(n_panels+1) < max_points
    warndlg(['The number of panels requested for resampling is less than ' ...
             'the current number of panels for some airfoils in database'], ...
            'resample airfoil warning');
end

for kaf=1:n_airfoils
    xy=afdb(kaf).coords;
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
            xy(end+1,:) = xy(1,:); %#ok<AGROW>
        end
    end
    
    %Calculate arc length of xy points clockwise from trailing edge
    n_points=size(xy,1);
    t=zeros(n_points,1);
    for i=2:n_points
        %formula: t(i) = hypot( x(i)-x(i-1) , y(i)-y(i-1) ) + t(i-1);
        t(i) = hypot(xy(i,1)-xy(i-1,1),xy(i,2)-xy(i-1,2)) + t(i-1);
    end

    %Get total arc length
    arc_length=t(end);
    
    %Spline airfoil with many points
    oversample=100;
    manypoints=linspace(0,arc_length,oversample*n_points);
    spline_type='spline';
    switch spline_type
        case {'linear','pchip','spline'}
            try
            xxyy = interp1(t,xy,manypoints,spline_type);
            catch
                keyboard%ble
            end
        otherwise
            fprintf('Airfoil oversampling algorithm specified is not an available option. Defaulting to "spline".');
            xxyy = interp1(t,xy,manypoints,'spline');
    end
    
    %Separate into high and low pressure surfaces
    [min_xx min_point]=min(xxyy(:,1));
    HP=xxyy(1:min_point,:);  % HP points progress from TE (x=1) to LE (x=0)
    LP=xxyy(min_point:end,:);  % LP points progress from LE (x=0) to TE (x=1)
    %jcb: technially, the LE is at the point of max curvature, but that
    %definition can produce situations that break the interpolation step
    
    %Calculate x points based on spacing algorithm specified
    switch spacing
        case 'cosine'
            beta=linspace(0,pi,n_panels+1);
            x=0.5*(1-cos(beta));
        case 'half-cosine'
            beta=linspace(0,pi/2,n_panels+1);
            x=(1-cos(beta));
        case 'constant'
            x=linspace(0,1,n_panels+1);
        otherwise
            error('Resampling algorithm specified is not an available option');
            return;  %#ok
    end
    %Compensate for TE not at x=1 and LE not at x=0
    max_xx=xxyy(1,1);
    x = x(:); % "x(:)" ensures column vector
    x_fwd=x*(max_xx-min_xx)+min_xx;  % x_fwd values are nominally 0 to 1; 
    x_rev=flipud(x_fwd);  % x_rev values are nominally 1 to 0 

    %Calculate interpolated airfoil points. For sharp trailing edge airfoils,
    %the trailing edge point is not repeated
    LP_new=[x_fwd  interp1(LP(:,1),LP(:,2),x_fwd)];
    HP_new=[x_rev  interp1(HP(:,1),HP(:,2),x_rev)];
    % Make sure that LE point is at (0,0) - Warning, this may have
    % complications when it comes to CFD analyses!!
    % HP_new(end,:)=[0 0];
    %JCB: 2012-8-24, forcing LE at (0,0) has caused spike at LE in some 
    %models e.g. CX-100 airfoil NPS_1200_100f
    
    %Assemble the two curves into a continuous line
    afdb(kaf).xy=[xyTE; HP_new(:,:); LP_new(2:end,:); xyTE];
    x_chord     =[1; flipud(x); x(2:end); 1];
    switch te_type
        case 'flatback'
            TEtype = 'flat';
        case 'normal'
            hpTEslope = diff(xy(1:2,2)) / diff(xy(1:2,1));
            lpTEslope = diff(xy(end-1:end,2)) / diff(xy(end-1:end,1));
            % jcb: attempt to guess 'round' or 'sharp' based on difference
            %      of slopes;  2.00 is arbitrary threshold
            if hpTEslope-lpTEslope > 2.00   
                TEtype = 'round';
            else
                TEtype = 'sharp';
            end
    end
    if isfield(afdb(kaf),'AirfoilName')
        % TEtype should already be defined for stations
        afname = afdb(kaf).AirfoilName; %#ok<NASGU>
    else
        % define TEtype for airfoil database entries
        afdb(kaf).TEtype = TEtype;
        afname = afdb(kaf).name; %#ok<NASGU>
    end
    
    %Get index of leading edge point in resampled coordinates
    LE = 1 + (n_panels + 1);  % or 1 + length(HP_new)
    afdb(kaf).LE = LE;
    
    %Calculate arc length of xy points clockwise from trailing edge
    n_points=size(afdb(kaf).xy,1);
    t=zeros(n_points,1);
    for i=2:n_points
        %formula: t(i) = hypot( x(i)-x(i-1) , y(i)-y(i-1) ) + t(i-1);
        t(i) = hypot(afdb(kaf).xy(i,1)-afdb(kaf).xy(i-1,1),afdb(kaf).xy(i,2)-afdb(kaf).xy(i-1,2)) + t(i-1);
    end
    t = t - t(LE);  % offset so that leading edge is datum
    afdb(kaf).s = t;  % signed arc length (HP-,LP+) for interpolation
    
    % signed chord (HP-,LP+) for interpolation
    afdb(kaf).c = [-1*x_chord(1:LE-1); x_chord(LE:end)];
    
    % signed chord of original airfoil coordinates
    % jcb: these calculations are not verified and are not currently used
    % elsewhere in NuMAD. I started writing these lines when dealing with
    % leading edge shell7 problems on 2012-08-28
%     oc = afdb(kaf).coords(:,1)*(max_xx-min_xx)+min_xx;
%     k = find((oc < 0.1) & (afdb(kaf).coords(:,2) > afdb(kaf).xy(LE,2)),1);
%     afdb(kaf).oc = [-1*oc(1:k-1); oc(k:end)];
    
    
    if 0  % debugging plots
        figure(701); %#ok<UNRCH>
%         subplot(2,1,1);
%         plot(t,xy(:,1),'k-x',t,xy(:,2),'r-x',manypoints,xxyy(:,1),'b-',manypoints,xxyy(:,2),'g-')
%         xlabel('airfoil surface distance from first point')
%         ylabel('oversampled coordinate value')
%         legend('af\_in x','af\_in y','oversampled x','oversampled y','Location','North')
%         subplot(2,1,2)
        plot(xy(:,1),xy(:,2),'r-x',...
            xxyy(:,1),xxyy(:,2),'b',...
            afdb(kaf).xy(:,1),afdb(kaf).xy(:,2),'k-o',...
            afdb(kaf).xy(LE,1),afdb(kaf).xy(LE,2),'c*')
%plot oc:  afdb(kaf).coords(:,1),interp1(afdb(kaf).c(2:end-1),afdb(kaf).xy(2:end-1,2),afdb(kaf).oc),'bs',...
        axis equal
        legend('af\_in points',['af\_in points oversampled by factor of ' num2str(oversample)],'af\_out','Location','SouthEast')
        title(afname,'interpreter','none');
        pause
    end

end