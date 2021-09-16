function output = BPEPost2FASTBlade(varargin)
%FUNCTIONNAME  One-line function description
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   output = functionName(input)
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%      input  = description of inputs
%  
%      output = description of outputs
%

%===== CREDITS & CHANGELOG ================================================
% 2011.06.09  brr: initial release
% yyyy.mm.dd  initials: description

load('BPE_SectionData')
load('polyfitdata')
N=length(zbeamnode_mid);

% Calculate AeroCent at desired locations
if ischar(varargin{1})
    nmdfn = varargin{1};
    a=readNuMADinput(nmdfn);
else
    a=varargin{1};  % assume the 'stations' structure from NuMAD is provided
end
for i=1:length(a)
    tmp(i,:)=[a(i).LocationZ a(i).AeroCenter a(i).Xoffset];
end
tmp(:,4)=0.25-tmp(:,3)+tmp(:,2); % FAST AeroCenter = 0.25-PitchAx + Real AeroCenter
AeroCent=interp1(tmp(:,1),tmp(:,4),zbeamnode_mid,'linear'); 
disp('X-offset locations:')
interp1(tmp(:,1),tmp(:,3),zbeamnode_mid,'linear')

% fill in elements of blade data structure
blade.title{1}=sprintf('Properties generated using BPE 1.31 and ANSYS on %s',datestr(now));

blade.NBlInpSt=N;
blade.CalcBMode='F';

blade.BldFlDmp(1)=1.5;
blade.BldFlDmp(2)=1.5;
blade.BldEdDmp(1)=1.5;

blade.FlStTunr(1)=1;
blade.FlStTunr(2)=1;
blade.AdjBlMs=1;
blade.AdjFlSt=1;
blade.AdjEdSt=1;

blade.prop.BlFract=zbeamnode_mid/zbeamnode_mid(end);
blade.prop.AeroCent=AeroCent;  % need to get this right
%blade.prop.StrcTwst=rot_i*57.3*rotdir;  % need to verify the signs in this one
blade.prop.StrcTwst=rot_aero*57.3*-1*rotdir;
blade.prop.BMassDen=massperlen;
blade.prop.FlpStff=EI_flap;
blade.prop.EdgStff=EI_edge;
blade.prop.GJStff=GJ;
blade.prop.EAStff=EA;
blade.prop.Alpha=Alphaa;
blade.prop.FlpIner=Iyyperlen_cg;
blade.prop.EdgIner=Ixxperlen_cg;
blade.prop.PrecrvRef=0*ones(N,1);
blade.prop.PreswpRef=0*ones(N,1);
blade.prop.FlpcgOf=xcg;  % switch of x-y coordinate system
blade.prop.EdgcgOf=ycg*rotdir;  % switch of x-y coordinate system
blade.prop.FlpEAOf=x_el_off_sect;  % switch of x-y coordinate system
blade.prop.EdgEAOf=y_el_off_sect*rotdir;  % switch of x-y coordinate system
% blade.prop.FlpEAOf=x_sh_off_sect;  % switch of x-y coordinate system
% blade.prop.EdgEAOf=y_sh_off_sect*rotdir;  % switch of x-y coordinate system

% Use blade data structure in FAST blade file
writeFastBlade(blade,'FASTBlade_bpe.dat')

% plot data for review
names = fieldnames(blade.prop);
logplotnames='BMassDen FlpStff EdgStff GJStff EAStff';
N=length(names);
m=floor(sqrt(N));
n=floor(N/m);
figure('Name','FAST/ADAMS Blade File Parameter Inspection')
for j=2:length(names)
    subplot(m,n,j-1)
    x=getfield(blade.prop,names{1});
    y=getfield(blade.prop,names{j});
    if strfind(logplotnames,names{j})
        semilogy(x,y,'b-s')
    else
        plot(x,y,'b-s')
    end
    xlabel(names{1})
    ylabel(names{j})
end

end