function [ output_args ] = getLCS( adfn )
%GETLCS Compute lift curve slope of AeroDyn airfoil performance data
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   LCS = getLCS(adfn)
%   Created for use with AeroDyn v13.00.00
% 
%   Lift curves for all airfoil performance data files listed in the 
%   AeroDyn input file are plotted. User must click the cursor at minimum 
%   and maximum AOA for straight line fit of each lift curve slope.
% 
%   adfn = Name of AeroDyn input file containing pointers to N individual 
%       airfoil performance files 
%   LCS = vector of length N with LCS estimates for each airfoil in the
%       AeroDyn input file
% 

ad=readFastAD(adfn);% get LCS values

for i=1:length(ad.FoilNm)
    foil=ad.FoilNm(i);
    foil=char(strrep(foil,'"',''));
    af=readAirfoilData(foil);
    x=af.AoA;
    y=af.CL;
    plot(x,y,'k-x')
    title('Choose two points to define range for linear fit')
    d=ginput(2);
    d=[d(1,1) d(2,1)];
    d=sort(d);
    pointer1=min(find(x>=d(1)));
    pointer2=max(find(x<=d(2)));
    xs=x(pointer1:pointer2);
    ys=y(pointer1:pointer2);
    p=polyfit(xs,ys,1);
    LCSvals(i)=p(1)*180/pi;
    xn=[xs(1) xs(end)];
    yn=polyval(p,xn);
    plot(x,y,'k-x',xn,yn,'r')
    pause(1)
end

% for i=1:length(ad.NFoil)
%     tmp(i)=LCSvals(ad.NFoil(i));
% end

LCS=LCSvals;

save LCSvals LCS

end

