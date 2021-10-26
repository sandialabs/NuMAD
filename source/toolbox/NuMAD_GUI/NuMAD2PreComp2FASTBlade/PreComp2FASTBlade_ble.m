function PreComp2FASTBlade_ble(varargin)

load PreComp_SectionData
load polyfitdata
N=size(PreComp_SectionData,1);

data = varargin{1};  % expects "data" from NuMAD
a=data.station;
ref=data.blade;
ref=calcGenLinePP(ref);

% Calculate AeroCent, sweep, and curve at desired locations
L=a(end).LocationZ;
for i=1:length(a)
    tmp(i,:)=[a(i).LocationZ a(i).AeroCenter a(i).Xoffset];
end
tmp(:,1)=tmp(:,1)./L; % convert to BlFract from meters
tmp(:,4)=0.25-tmp(:,3)+tmp(:,2); % FAST AeroCenter = 0.25-PitchAx + Real AeroCenter
AeroCent=interp1(tmp(:,1),tmp(:,4),PreComp_SectionData(:,1),'linear');
sweep = ppval(ref.PresweepRef.pp,PreComp_SectionData(:,1)*L);
curve = ppval(ref.PrecurveRef.pp,PreComp_SectionData(:,1)*L);

% fill in elements of blade data structure
blade.title{1}=['Properties generated using NuMAD2PreComp2FASTBlade on ' date];

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

blade.prop.BlFract=PreComp_SectionData(:,1);
blade.prop.AeroCent=AeroCent;
blade.prop.StrcTwst=PreComp_SectionData(:,3);  % degrees
blade.prop.BMassDen=PreComp_SectionData(:,18);
blade.prop.FlpStff=PreComp_SectionData(:,4);
blade.prop.EdgStff=PreComp_SectionData(:,5);
blade.prop.GJStff=PreComp_SectionData(:,6);
blade.prop.EAStff=PreComp_SectionData(:,7);
blade.prop.Alpha=PreComp_SectionData(:,11);  % need to check this out
blade.prop.FlpIner=PreComp_SectionData(:,19);
blade.prop.EdgIner=PreComp_SectionData(:,20);
blade.prop.PrecrvRef=curve;
blade.prop.PreswpRef=sweep;
blade.prop.FlpcgOf=PreComp_SectionData(:,22);
blade.prop.EdgcgOf=PreComp_SectionData(:,23);
if 1  % using tension center as elastic center
    blade.prop.FlpEAOf=PreComp_SectionData(:,16);
    blade.prop.EdgEAOf=PreComp_SectionData(:,17);  % switch of x-y coordinate system
else  % using shear center as elastic center
    blade.prop.FlpEAOf=PreComp_SectionData(:,14);
    blade.prop.EdgEAOf=PreComp_SectionData(:,15);  % switch of x-y coordinate system
end

% Use blade data structure in FAST blade file
writeFastBlade(blade,'FASTBlade_precomp.dat')

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