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

% error('calculate blade mass somewhere in this script')
bladeMass = 0;

%% ************************************************************************
% ================= CREATE NODAL LISTING FROM ANSYS MESH =================
% run make_nlist.mac

ansysFilename = config.ansys.analysisFileName;
ansys_path = blade.paths.ansys;
ansys_product = 'ANSYS';

script_name='makeNLIST.mac';
script_out='output1.txt';

fid=fopen(script_name,'w+');
fprintf(fid,'resume,master,db\n');
fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);   %From master, change the jobname

fprintf(fid,'!!! BEGIN MAKE_NLIST MACRO TEXT\n');
fprintf(fid,'ESEL,S,SEC,,1,999   \n');
% fprintf(fid,'ALLSEL   \n');
fprintf(fid,'NSLE,S  \n');
% fprintf(fid,'nsel,u,node,,z_master_node_number\n');
fprintf(fid,'/output,NLIST,lis\n');
fprintf(fid,'/page,1e6,,1e6,,\n');
fprintf(fid,'NLIST,ALL, , ,XYZ,NODE,NODE,NODE\n');
fprintf(fid,'/output,\n');
fprintf(fid,'ALLSEL,ALL \n');

% fprintf(fid,'allsel\n'); %Select all nodes
% fprintf(fid,'nsel, all, node\n'); %Select all nodes
% fprintf(fid,'esel, all, elem\n'); %Select all elements
% fprintf(fid,'*get, nnode, node,0,count\n'); %Find the number (quantity) of nodes selected
% fprintf(fid,'*get, nel, ELEM,0,count\n'); %Find the number (quantity) of nodes selected

% fprintf(fid,'/output,nel,txt\n');
% fprintf(fid,'*status, nel\n');
% fprintf(fid,'/output\n');
% 
% fprintf(fid,'/output,nnode,txt\n');
% fprintf(fid,'*status, nnode\n');
% fprintf(fid,'/output\n');

fprintf(fid,'!!! END MAKE_NLIST TEXT\n');
fclose(fid);

ansys_call = sprintf('"%s" -b -p %s -I %s -o %s',ansys_path,ansys_product,script_name,script_out);
disp(ansys_call)%ble
pse=3;
while 1
    [status,~] = dos(ansys_call);  % the windows system call to run the above ansys command
    if status==0
        disp(' ')
        disp('Node listing created')
        %delete(script_name);
        break;
    end
    fprintf('%s: layupDesign_ANSYSmesh: Waiting for ANSYS to create nlist... ...\n',datestr(now))
    pause(pse);
end

