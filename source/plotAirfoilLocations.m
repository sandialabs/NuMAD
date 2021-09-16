function output = plotAirfoilLocations(ADFilename,FASTFilename)
%PLOTAIRFOILLOCATIONS Create illustration of airfoil placement based on
%   AeroDyn input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   output = plotAirfoilLocations(ADFilename)
%   Supply the full name of an AeroDyn input file.  A plot is produced
%   which shows the distribution of airfoils along the span of the blade.
%
%      ADFilename = full AeroDyn input filename; eg 'Micon_AD.ipt'
%      FASTFilename = Optional; full FAST input filename; eg 'Micon.fst';
%               used to get hub radius for aspect ratio calculation
%

ad=readFastAD(ADFilename);

edges=ad.RNodes-ad.DRNodes/2;
tip=ad.RNodes(end)+ad.DRNodes(end)/2;
edges=[edges ; tip];
colors={'r','b','g','c','m','k'};
for i=1:ad.BldNodes
    patch([edges(i) edges(i+1) edges(i+1) edges(i)],[0 0 ad.Chord(i) ad.Chord(i)],colors{ad.NFoil(i)})
    area(i)=ad.Chord(i)*(edges(i+1)-edges(i));
    hold on
end
axis equal
hold off
names={};
for i=1:length(ad.FoilNm)
    tmp=regexp(ad.FoilNm{i},'\\(.*)"','tokens');
    tmp=tmp{1};
    tmp=strrep(tmp,'_','\_');
    names=[names tmp{1}];
end

text(2,-2,names)
xlabel('meters')
ylabel('meters')

if exist('FASTFilename','var')
    fast=readFastMain(FASTFilename);
    L=tip-fast.TurbConf.HubRad;
    A=sum(area);
    disp(['Blade Length = ' num2str(L)])
    disp(['Area ~= ' num2str(A)])
    disp(['Aspect ratio = ' num2str(L^2/A)])
end