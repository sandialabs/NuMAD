function [designvar] = layupDesignFASTmass(blade)

disp('Creating FAST Blade file using NuMAD and PreComp...')
numad('numad.nmd','precomp',[1 3 2])
bld=readFastBlade('FASTBlade_precomp.dat');

tmp=(bld.prop.BMassDen(1:end-1)+bld.prop.BMassDen(2:end))/2;
mass_per_meter=sum(tmp.*diff(bld.prop.BlFract));
designvar=mass_per_meter;


end