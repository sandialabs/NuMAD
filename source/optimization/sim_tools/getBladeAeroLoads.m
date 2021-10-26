function output = getBladeAeroLoads(adfile,outfile,elmfile,result_type)
%GETBLADEAEROLOADS Get magnitude and location of blade aero loads (AeroDyn)
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   output = getBladeAeroLoads(adfile,outfile,elmfile,result_type)
%
%   Computes flapwise and edgewise loads on each aerodynamic element.
%   adfile = file name of AeroDyn input
%   outfile = FAST *.out file name
%   elmfile = FAST/AeroDyn *.elm file name (created when AeroDyn "PrnElm"
%         is set to PRINT)
%   result_type = 
%           'mean' = compute and report mean aero load for entire time
%               history
%           'maxRootM' = compute and report aero load from entire time
%               history that yields the highest blade root moment.
% 
%   *.out must contain a column for 'BldPitch'
%   Writes 'forces.forces' file which contains the following columns:
%     Z, Fx, Fy, M, Alpha, x_off, y_off
%   **Note: Alpha, x_off, and y_off inactive as of 9/12/2012
%

ad = readFastAD(adfile);
out = loadFASTOutData(outfile);
elm = loadElmData(elmfile);

colPitch = find(1==strncmp('BldPitch',out.list,8),1);

aero.R = ad.RNodes;
aero.Z = ad.RNodes - 0.5*ad.DRNodes;
aero.Z(end+1) = ad.RNodes(end) + 0.5*ad.DRNodes(end);
aero.R = aero.R - aero.Z(1);  % subtract hub radius
aero.Z = aero.Z - aero.Z(1);  % subtract hub radius
aero.Z(1) = [];     % first entry not needed

FT = elm.data(:,elm.ForcT);
FN = elm.data(:,elm.ForcN);
[Fflap Fedge] = deal(zeros(size(FT)));
for k=1:size(FT,1)
    adtime = elm.data(k,1);
    theta = pi/180*interp1(out.data(:,1),out.data(:,colPitch),adtime,'linear','extrap');
    Fflap(k,:) =   FN(k,:)*cos(theta) + FT(k,:)*sin(theta);
    Fedge(k,:) = (-FN(k,:)*sin(theta) + FT(k,:)*cos(theta))*-1;
end

switch result_type
    case 'mean'
        aero.Fy = mean(Fflap);
        aero.Fx = mean(Fedge);
        [aero.M aero.alpha aero.xoff aero.yoff] = deal(zeros(size(aero.Fx)));
        aero.RootM = aero.Fy*aero.R;
    case 'maxRootM'
        RootM = Fflap*aero.R;
        [y,i] = max(RootM);
        aero.RootM = y;
        aero.i = i;
        aero.Fy = Fflap(i,:);
        aero.Fx = Fedge(i,:);
        [aero.M aero.alpha aero.xoff aero.yoff] = deal(zeros(size(aero.Fx)));
end

fid=fopen('forces.forces','wt');
fprintf(fid,'Z, Fx, Fy, M, Alpha, x_off, y_off\n');

for i=1:length(aero.Z)
    fprintf(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\n',aero.R(i),aero.Fx(i),aero.Fy(i),aero.M(i),0,0,0);
end

fclose all;

end