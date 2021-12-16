    function  [newaf,newchord] = interpAirfoil( coords, chord, spanloc, newspanloc)
% INTERPAIRFOIL Interpolate existing shapes to form new airfoil shape.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% output = interpAirfoil( coords1, coords2, spanlocs, method)
%   Generate airfoil coordinates representing an interpolation
%  between two different airfoils.  Interpolation algorithm is determined
%  by the number of airfoils provided.  If information for two airfoils is
%  provided, then a linear interpolation is performed. If four
%  airfoils are provided, then a combination of pchip and splines are used
%  to create the new shape.
%
% Interpolation is performed as follows:
%  1. Determine camber lines and thickness distribution for supplied airfoils
%  2. Interpolate to find camber line, chord, percent thickness distribution,
%     and absolute thickness for new airfoil.
%  3. Scale percent thickness distribution of new airfoil in order to
%     ensure that it is at the desired thickness/chord to fit in with 
%     blade geometry.
%  4. Reconstruct new airfoil based on new camber line and new thickness
%     distribution
%
%  coords  = x,y pairs of normalized coordinates for airfoils; format as follows:
%                [ x1 y1 x2 y2]   - method 1
%                     or,
%                [ x1 y1 x2 y2 x3 y3 x4 y4]  -  method 2
%              where each column contains N coordinate data points; note
%              that this format assumes the airfoils have been resampled
%              into points that map between all airfoils.
%            Method 1--Linear interpolation between two adjacent airfoils
%            Method 2--Interpolation using additional pair of adjacent
%                      airfoils
%  spanloc = 1x2 or 1x4 array of span locations for existing airfoils
%            described by coords
%  newspanloc = Span location for new airfoil; must fall in middle of all
%               given airfoils
%
%  newaf    = x,y pairs of coordinates for interpolated airfoil at location
%             determined by newspanloc
%  newchord = chord length for new airfoil
%

%===== CREDITS & CHANGELOG ================================================
%  1/25/2011 BRR: Created

% Error Checking
if size(coords,2)/2~=size(spanloc,2)
    error('Number of supplied airfoils must equal number of supplied span locations.')
    return
else
    Npt=size(coords,1);
    Naf=size(coords,2)/2;
end

switch Naf
    case 2
        if (newspanloc>spanloc(2) || newspanloc<spanloc(1))
            error('Quit 1')
            return
        end
    case 4
        if (newspanloc>spanloc(3) || newspanloc<spanloc(2))
            error('Quit 2')
            return
        end
end
% End Error Checking

% find camber line and thickness of airfoils
xcamb=coords(1:Npt/2+1,1);
camb=zeros(Npt/2+1,Naf);
perthick=zeros(Npt/2+1,Naf);
for i=2:Npt/2
    switch Naf
        case 2
            camb(i,:) = ( coords(i,[2 4]) + coords(end-i+2,[2 4]) ) /2;
            perthick(i,:) = (camb(i,:)-coords(i,[2 4]) )*2;
        case 4
            camb(i,:) = ( coords(i,[2 4 6 8]) + coords(end-i+2,[2 4 6 8]) ) /2;
            perthick(i,:) = (camb(i,:)-coords(i,[2 4 6 8]) )*2;
    end
end

% Interpolate to find camber line, chord, percent thickness distribution,
%  and absolute thickness for new airfoil
ycamb=zeros(1,Npt/2+1);
yperthick=zeros(1,Npt/2+1);
for i=1:Npt/2+1
    switch Naf
        case 2
            yabsthick=interp1(spanloc,max(perthick).*chord,newspanloc,'linear');
            ychord=interp1(spanloc,chord,newspanloc,'linear');
            ycamb(i)=interp1(spanloc,camb(i,:),newspanloc,'linear');
            yperthick(i)=interp1(spanloc,perthick(i,:),newspanloc,'linear');
        case 4
            yabsthick=interp1(spanloc,max(perthick).*chord,newspanloc,'pchip');
            ychord=interp1(spanloc,chord,newspanloc,'pchip');
            ycamb(i)=interp1(spanloc,camb(i,:),newspanloc,'spline');
            yperthick(i)=interp1(spanloc,perthick(i,:),newspanloc,'pchip');
    end
end

% scale the thickness to match the desired t/c
toverc=yabsthick./ychord;
yperthick=toverc/max(yperthick)*yperthick;

% reconstruct the new airfoil based on new camber line and new thickness distribution
hp=ycamb-yperthick/2;
lp=ycamb+yperthick/2;
y=[hp' ; lp(end-1:-1:2)'];
x=coords(:,1);
af3=[x  y];

if 0
    figure(2003)
    plot(coords(:,1),coords(:,4),'b--')
    hold on
    plot(coords(:,1),coords(:,6),'g--')
    plot(xcamb,ycamb,'k')
    plot(x,y,'r')
    hold off
end

newaf=af3;
newchord=ychord;

end
