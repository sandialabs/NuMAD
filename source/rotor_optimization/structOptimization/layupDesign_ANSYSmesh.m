function layupDesign_ANSYSmesh(blade,config)
% Generate the ANSYS mesh and nodal listing

disp(' '); disp('Creating ANSYS model...')
fprintf('Mesh size setting = %0.4f\n',blade.mesh)
while 1
    numad('numad.nmd','ansys',blade.mesh) 
    %numad_multiLayer('numad.nmd','ansys',blade.mesh) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%TEMP
    if exist('master.db','file')
        disp('ANSYS model created')
        break;
    end
    fprintf('%s: layupDesign_ANSYSmesh: Waiting for NUMAD to create ANSYS model...\n',datestr(now))
    pause(3);
end


