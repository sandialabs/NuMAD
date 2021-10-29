function PreComp2FASTBlade_BladeDef(blade,PreComp_SectionData,fastBlade,BATCH_RUN)

N=size(PreComp_SectionData,1);
% Calculate AeroCent, sweep, and curve at desired locations
L = blade.ispan(end);
for i=1:length(blade.ispan)
    tmp(i,:)=[blade.ispan(i) blade.iaerocenter(i) blade.xoffset(i)];
end
tmp(:,1)=tmp(:,1)./L; % convert to BlFract from meters
tmp(:,4)=0.25-tmp(:,3)+tmp(:,2); % FAST AeroCenter = 0.25-PitchAx + Real AeroCenter
AeroCent=interp1(tmp(:,1),tmp(:,4),PreComp_SectionData(:,1),'linear');

% ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% change code from NuMAD to use bladeDef variables
% ref=data.blade;
% ref=calcGenLinePP(ref);
% sweep = ppval(ref.PresweepRef.pp,PreComp_SectionData(:,1)*L);
% curve = ppval(ref.PrecurveRef.pp,PreComp_SectionData(:,1)*L);
sweep = interp1(blade.ispan,blade.isweep,PreComp_SectionData(:,1),'pchip');
curve = interp1(blade.ispan,blade.iprebend,PreComp_SectionData(:,1),'pchip');
% ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% fill in elements of blade data structure
fastBlade.title{1}=['Properties generated using NuMAD2PreComp2FASTBlade on ' date];

fastBlade.NBlInpSt=N;
fastBlade.CalcBMode='F';

fastBlade.BldFlDmp(1)=1.5;
fastBlade.BldFlDmp(2)=1.5;
fastBlade.BldEdDmp(1)=1.5;

fastBlade.FlStTunr(1)=1;
fastBlade.FlStTunr(2)=1;
fastBlade.AdjBlMs=1;
fastBlade.AdjFlSt=1;
fastBlade.AdjEdSt=1;

fastBlade.prop.BlFract=PreComp_SectionData(:,1);
fastBlade.prop.AeroCent=AeroCent;
fastBlade.prop.StrcTwst=PreComp_SectionData(:,3);  % degrees
fastBlade.prop.BMassDen=PreComp_SectionData(:,18);
fastBlade.prop.FlpStff=PreComp_SectionData(:,4);
fastBlade.prop.EdgStff=PreComp_SectionData(:,5);
fastBlade.prop.GJStff=PreComp_SectionData(:,6);
fastBlade.prop.EAStff=PreComp_SectionData(:,7);
fastBlade.prop.Alpha=PreComp_SectionData(:,11);  % need to check this out
fastBlade.prop.FlpIner=PreComp_SectionData(:,19);
fastBlade.prop.EdgIner=PreComp_SectionData(:,20);
fastBlade.prop.PrecrvRef=curve;
fastBlade.prop.PreswpRef=sweep;
fastBlade.prop.FlpcgOf=PreComp_SectionData(:,22);
fastBlade.prop.EdgcgOf=PreComp_SectionData(:,23);
if 1  % using tension center as elastic center
    fastBlade.prop.FlpEAOf=PreComp_SectionData(:,16);
    fastBlade.prop.EdgEAOf=PreComp_SectionData(:,17);  % switch of x-y coordinate system
else  % using shear center as elastic center
    fastBlade.prop.FlpEAOf=PreComp_SectionData(:,14);
    fastBlade.prop.EdgEAOf=PreComp_SectionData(:,15);  % switch of x-y coordinate system
end

% Use blade data structure in FAST blade file
writeFastBlade(fastBlade,'FASTBlade_precomp.dat')

% plot data for review
if ~BATCH_RUN
    names = fieldnames(fastBlade.prop);
    logplotnames='BMassDen FlpStff EdgStff GJStff EAStff';
    N=length(names);
    m=floor(sqrt(N));
    n=floor(N/m);
    figure('Name','FAST/ADAMS Blade File Parameter Inspection')
    for j=2:length(names)
        subplot(m,n,j-1)
        x=getfield(fastBlade.prop,names{1});
        y=getfield(fastBlade.prop,names{j});
        if strfind(logplotnames,names{j})
            semilogy(x,y,'b-s')
        else
            plot(x,y,'b-s')
        end
        xlabel(names{1})
        ylabel(names{j})
    end
end

end