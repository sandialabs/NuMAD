function  newaf = scaleAirfoilThickness(coords,newthickness)
% SCALEAIRFOILTHICKNESS Scale thicknes of airfoil shape.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% newaf = scaleAirfoilThickness( coords, newthickness )
%   Determine existing airfoil thickness and then scale to a new thickness.
%  Shapes are scaled about the airfoil camber line.
%
%  coords  = x,y pairs of normalized coordinates for airfoils where each
%            column contains N coordinate data points; note that this
%            format assumes the airfoils have been resampled into chord
%            coordinates that are common between HP and LP surfaces
%  newthickness = desired new airfoil thickness in percent; eg a 23%
%            airfoil would take 23 as newthickness, not 0.23
%  newaf    = x,y pairs of coordinates for scaled airfoil 


% find camber line and thickness of airfoils
Npt=length(coords);
xcamb=coords(1:Npt/2+1,1);
camb=zeros(Npt/2+1,1);
perthick=zeros(Npt/2+1,1);

for i=2:Npt/2
    camb(i) = ( coords(i,2) + coords(end-i+2,2) ) /2;
    perthick(i) = (camb(i)-coords(i,2) )*2;
end

t=max(perthick)*100;
% scale the thickness to match the desired t/c
scalefactor=newthickness/t;

% reconstruct the new airfoil based on new camber line and new thickness distribution
hp=camb-scalefactor*perthick/2;
hp=hp';
lp=camb+scalefactor*perthick/2;
lp=lp';
y=[hp' ; lp(end-1:-1:2)'];
x=coords(:,1);
af3=[x  y];

if 0
    figure(2003)
    plot(coords(:,1),coords(:,2),xcamb,camb,'--',xcamb,perthick,x,y,'k-o',xcamb,perthick*scalefactor,'k')
    legend('Original airfoil','Camber line','Original % thickness distribution','Scaled airfoil','New % thickness distribution')
end

newaf=af3;

end
