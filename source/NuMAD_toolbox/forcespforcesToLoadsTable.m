function loads = forcespforcesToLoadsTable(forcesfile)
%Reads the old forces.forces file and transforms to 
%the new loads_table format. The code assumes that
%forces.forces is in the ANSYS coordinate system
%and the loads_table is in FAST blade coordinates.


% open file selector if 'forcesfile' not given
if ~exist('forcesfile','var') || isempty(forcesfile)
    [fn,pn] = uigetfile( ...
        {'*.forces','AD line load (*.forces)'; ...
        '*.*','All files (*.*)'},...
        'Select the forces file');
    if isequal(fn,0) || isequal(pn,0)
        disp('Operation canceled by user.')
        return;
    end
    forcesfile = fullfile(pn,fn);
end
% read the forces file
forces = read_forces(forcesfile);

loads.rBlade = forces{1}';
loads.Fxb = forces{2}';
loads.Fyb = forces{3}';
loads.Fzb = zeros(1,length(forces{1}));
loads.Mxb = zeros(1,length(forces{1}));
loads.Myb = zeros(1,length(forces{1}));
loads.Mzb = forces{4}';
loads.Alpha = forces{5}';
loads.presweep = forces{6}';
loads.prebend = forces{7}';

%%%%Transform from ANSYS coordinates to FAST

beta=[0 -1 0; %Direction cosine matrix for a 90 degree clockwise rotation
      1 0 0;  %about the z axis. 
      0 0 1]; 
for i=1:length(loads.rBlade)
    F=beta'*[loads.Fxb(i); loads.Fyb(i); loads.Fzb(i)];
    M=beta'*[loads.Mxb(i); loads.Myb(i); loads.Mzb(i)];


    %Overwrite loads in the FAST CSYS with the ANSYS CSYS
    loads.Fxb(i)=F(1);
    loads.Fyb(i)=F(2);
    loads.Fzb(i)=F(3);

    loads.Mxb(i)=M(1);
    loads.Myb(i)=M(2);
    loads.Mzb(i)=M(3);
end




%Back calculate input span and moment
z0=0;
zout=loads.rBlade;
rBladeForce_x=loads.Fxb;
rBladeForce_y=loads.Fyb;

rBladeMoment=zeros(numel(rBladeForce_x),1);
for i=1:numel(rBladeMoment)
   rBladeMoment(i)=mean([zout(i),z0]);
   z0=zout(i);
end


bladeLength=100;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recheck Moment Dist
rBlade = [rBladeMoment; bladeLength];

A=zeros(numel(zout));
for i=1:numel(zout)
    A(i,i:end)=(zout(i:end)-rBlade(i))';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
loads.input.rGage=rBladeMoment'; 
loads.input.Mxb=A*rBladeForce_y';
loads.input.Myb=A*rBladeForce_x';
%<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

function forces = read_forces(filename)

% Open the file
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
header = fgetl(fid);  %#ok (header not used)
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

% 'Z (m)	Fx (N)	Fy (N)	M (N-m)	Alpha	x_off	y_off'
forces = textscan(filecontents,repmat('%f',1,7));
end

end
