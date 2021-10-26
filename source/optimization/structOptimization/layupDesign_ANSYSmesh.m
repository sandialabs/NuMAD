function layupDesign_ANSYSmesh(blade,varargin)
% Generate the ANSYS mesh and nodal listing

defaultfn = 'numad.nmd';


% hard-code filename if not specified
if nargin==1
    fileName = defaultfn;
else
    fileName = varargin{1};
end

disp(' '); disp('Creating ANSYS model...')
fprintf('Mesh size setting = %0.4f\n',blade.mesh)
while 1
    numad(fileName,'ansys',blade.mesh) 
    if exist('master.db','file')
        disp('ANSYS model created')
        break;
    end
    fprintf('%s: layupDesign_ANSYSmesh: Waiting for NUMAD to create ANSYS model...\n',datestr(now))
    pause(3);
end


