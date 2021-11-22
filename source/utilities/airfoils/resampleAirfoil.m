function af_out = resampleAirfoil(af_in, n_panels, spacing)
% RESAMPLEAIRFOIL  Resample airfoil coordinates
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% af_out = resampleAirfoil(af_in, n_panels, spacing)
%
%       af_in = Nx2 array containing N normalized xy airfoil points
%    n_panels = number of panels to be created per surface
%     spacing = spacing routine to be used
%                 'cosine', 'half-cosine', 'constant', %ble: 'arc'
%
%      af_out = array containing n_panels+1 airfoil points
%
% Cosine spacing: puts higher density of points at both LE and TE;
%       constant arc length point spacing around a perfect circle.
% Half-cosine spacing: puts higher density of points at LE and lesser
%       density of points at TE
% constant spacing: constant spacing of points along chord line
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
%       for 'pchip', but others can be used by changing the te_type setting
%       in the code.
%
% Assumes leading edge at (0,0) and trailing edge at (1,0)

% Error checking
% if airfoil coordinates are not in Nx2 array
if size(af_in,2)~=2
    tmpN=size(af_in,1);
    tmpM=size(af_in,2);
    warning(['af_in array was defined in ' num2str(tmpN) 'x' num2str(tmpM) ' array.  Automatically changing it to an ' num2str(tmpM) 'x' num2str(tmpN) ' array.'])
    af_in=af_in';
end
% End error checking routines

% Determine te_type
if af_in(2,1)==1 && af_in(end,1)==1  && (af_in(2,2)~=0 || af_in(end,2)~=0)
    % it is flatback
    te_type='flatback';
    xy=af_in(2:end,:);
else
    % not flatback.  Either sharp or round
    xy=[af_in; af_in(1,:)];
    te_type='normal';
end

n_points=size(xy,1);
t=zeros(1,n_points);

%Calculate arc length of xy points clockwise from trailing edge
for i=2:n_points
    t(i)=((xy(i,1)-xy(i-1,1))^2+(xy(i,2)-xy(i-1,2))^2)^.5+t(i-1);
end

%Calculate total arc length
length_xy=t(i);

%Spline airfoil with many points
oversample=100;
manypoints=linspace(0,length_xy,oversample*n_points);
spline_type='pchip';
switch spline_type
    case 'linear'
        pp = interp1(t,xy,'linear','pp');
        xxyy = ppval(pp, manypoints);
    case 'pchip'
        pp = pchip(t,xy');
        xxyy = ppval(pp, manypoints)';
    case 'spline'
        pp = spline(t,xy');
        xxyy = ppval(pp, manypoints)';
    otherwise
        error('Airfoil oversampling algorithm specified is not an available option');
        return;
end

%Separate into high and low pressure surfaces
min_xxyy=min(xxyy(:,1));
pointer=find(xxyy(:,1)==min_xxyy);
HP=xxyy(1:pointer,:);
LP=xxyy(pointer:end,:);

%Calculate x points based on spacing algorithm specified
switch spacing
    case 'cosine'
        beta=linspace(min_xxyy,pi,n_panels+1)';
        x=0.5*(1-cos(beta));
        x_rev=x(end:-1:1);
    case 'half-cosine'
        beta=linspace(min_xxyy,pi/2,n_panels+1)';
        x=1-cos(beta);
        x_rev=x(end:-1:1);
    case 'constant'
        x=linspace(min_xxyy,1,n_panels+1)';
        x_rev=x(end:-1:1);
    case 'arc'% ble
        % track the arc length to evenly sample to find x-values
        % constant for HP and LP surfaces
        xxyy_hp = xxyy(1:pointer,:);
        points_arc=size(xxyy_hp,1);        
        s=zeros(points_arc,1);
        for ii = 2:points_arc
            s(ii,1)=sqrt((xxyy(ii,1)-xxyy(ii-1,1))^2+(xxyy(ii,2)-xxyy(ii-1,2))^2)+s(ii-1);
        end
        beta=linspace(0,s(end),n_panels+1)';
        x_rev = interp1(s,xxyy_hp(:,1),beta);        
        x=x_rev(end:-1:1);
    otherwise
        error('Resampling algorithm specified is not an available option');
end

%Calculate interpolated airfoil points. For sharp trailing edge airfoils,
%the trailing edge point is not repeated
HP_new=[x_rev  interp1(HP(:,1),HP(:,2),x_rev)];
% Make sure that LE point is at (0,0) - Warning, this may have
% complications when it comes ot CFD analyses!!
HP_new(end,:)=[0 0];
LP_new=[x  interp1(LP(:,1),LP(:,2),x)];
% ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% % switch te_type
% %     case 'flatback'
% %         af_out=[1 0; HP_new(:,:); LP_new(2:end,:)];
% %     case 'normal'
% %         af_out=[1 0; HP_new(:,:); LP_new(2:end-1,:); 1 0];
% % end
switch te_type
    case 'flatback'
        af_out=[1 0; HP_new(:,:); LP_new(2:end,:)];
    case 'normal'
        af_out=[HP_new(:,:); LP_new(2:end-1,:)];
end
% ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

end