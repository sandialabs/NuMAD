function mass_per_meter = layupDesign_FASTmass_BLE(bladeFilename)

% disp('Creating FAST Blade file using NuMAD and PreComp...')
% numad('numad.nmd','precomp',[1 3 2])

hm=pwd;
cd ..

disp('Calculating blade mass from FAST blade input file.')
% read in blade input file and calculate blade mass
blade = readFastBlade(bladeFilename);
mass_per_meter = trapz(blade.prop.BlFract, blade.prop.BMassDen.*blade.AdjBlMs);

cd(hm)


end