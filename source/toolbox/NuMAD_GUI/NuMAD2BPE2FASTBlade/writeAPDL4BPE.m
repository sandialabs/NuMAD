function writeAPDL4BPE(loads)
% writeAPDL4BPE Write APDL code to use at end of shell7.src file for BPE analysis
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   writeAPDL4BPE(loads)
%
%   Write APDL code to use at end of shell7.src file for BPE analysis.
%
%    loads, Nx6 matrix of forces to apply in each of the N static analysis
%    runs.  Columns correspond to Fx Fy Fz Mx My Mz.  Default input is a
%    6x6 diagonal matrix of ones.
%

%===== CREDITS & CHANGELOG ================================================
% 2011.12.21  brr: initial creation of this script
% yyyy.mm.dd  initials: description

% hard-code the loads for now:
loads=diag(ones(6,1))*1e3;

% read in shell7 file from NuMAD (without BPE commands)
file1 = 'shell7.src';
file2 = 'shell7bpe.src';
a=dir(file1);
if isempty(a)
    errordlg(sprintf('Could not open file "%s"',file1),'Error');
    return
end
disp(sprintf('Using shell7 file created on %s',a.date))
disp(sprintf('Appending BPE APDL commands to "%s"',file1))
copyfile(file1,file2);
fid = fopen(file2,'at');

numcases=size(loads,1);

fprintf(fid,'\n\n!=========== Set up and perform static analyses to enable BPE analysis ===========\n\n');

fprintf(fid,'*get,MaxNode,node,,num,max\n');

fprintf(fid,'*dim,nodeNum,,MaxNode\n');
fprintf(fid,'*dim,xNode,,MaxNode\n');
fprintf(fid,'*dim,yNode,,MaxNode\n');
fprintf(fid,'*dim,zNode,,MaxNode\n');
fprintf(fid,'*dim,nMass,,MaxNode\n');

for j=1:numcases
    fprintf(fid,'*dim,D%ix,,MaxNode\n',j);
    fprintf(fid,'*dim,D%iy,,MaxNode\n',j);
    fprintf(fid,'*dim,D%iz,,MaxNode\n',j);
end
fprintf(fid,'\n');

fprintf(fid,'*do,j,1,MaxNode\n');
fprintf(fid,'  nodeNum(j)=j\n');
fprintf(fid,'  *get,temp,node,j,loc,x\n');
fprintf(fid,'  xNode(j)=temp\n');
fprintf(fid,'  *get,temp,node,j,loc,y\n');
fprintf(fid,'  yNode(j)=temp\n');
fprintf(fid,'  *get,temp,node,j,loc,z\n');
fprintf(fid,'  zNode(j)=temp\n');
fprintf(fid,'*enddo\n\n');

fprintf(fid,'/filname,utlsuite\n');
fprintf(fid,'/solution\n');
fprintf(fid,'lumpm,on\n\n');

fprintf(fid,'nsel,s,loc,z,0\n');
fprintf(fid,'d,all,all\n');
fprintf(fid,'nsel,all\n\n');

fprintf(fid,'sfedele,all,all,all\n\n');
fprintf(fid,'nlgeom,off\n\n');

for j=1:numcases
    fprintf(fid,'f,z_master_node_number,fx,%3.1f\n',loads(j,1));
    fprintf(fid,'f,z_master_node_number,fy,%3.1f\n',loads(j,2));
    fprintf(fid,'f,z_master_node_number,fz,%3.1f\n',loads(j,3));
    fprintf(fid,'f,z_master_node_number,mx,%3.1f\n',loads(j,4));
    fprintf(fid,'f,z_master_node_number,my,%3.1f\n',loads(j,5));
    fprintf(fid,'f,z_master_node_number,mz,%3.1f\n',loads(j,6));
    fprintf(fid,'allsel\n');
    fprintf(fid,'solve\n\n');
end

fprintf(fid,'lsclear,all\n');
fprintf(fid,'d,all,all\n');
fprintf(fid,'acel,1\n');
fprintf(fid,'solve\n\n');

fprintf(fid,'finish\n\n');
fprintf(fid,'/post1\n\n');

fprintf(fid,'!=========== Load results into output files ===========\n\n');

for j=1:numcases
    fprintf(fid,'set,%i\n\n',j);
    fprintf(fid,'*do,j,1,MaxNode\n');
    fprintf(fid,'  *get,temp,node,j,u,x\n');
    fprintf(fid,'  D%ix(j)=temp\n',j);
    fprintf(fid,'  *get,temp,node,j,u,y\n');
    fprintf(fid,'  D%iy(j)=temp\n',j);
    fprintf(fid,'  *get,temp,node,j,u,z\n');
    fprintf(fid,'  D%iz(j)=temp\n',j);
    fprintf(fid,'*enddo\n\n');
end

fprintf(fid,'set,%i\n',numcases+1);
fprintf(fid,'*vget,nMass(1),node,1,rf,fx\n\n');

fprintf(fid,'*cfopen,nodeloc,txt\n');
fprintf(fid,'*vwrite,nodeNum(1),xNode(1),yNode(1),zNode(1)\n');
fprintf(fid,'(F8.0,1X,E16.8,1X,E16.8,1X,E16.8)\n');
fprintf(fid,'*cfclos\n\n');

fprintf(fid,'*cfopen,fdisp,txt\n');
fprintf(fid,'*vwrite,D1x(1),D1y(1),D1z(1),D2x(1),D2y(1),D2z(1),D3x(1),D3y(1),D3z(1)\n');
fprintf(fid,'(E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16)\n');
fprintf(fid,'*cfclos\n\n');

fprintf(fid,'*cfopen,mdisp,txt\n');
fprintf(fid,'*vwrite,D4x(1),D4y(1),D4z(1),D5x(1),D5y(1),D5z(1),D6x(1),D6y(1),D6z(1)\n');
fprintf(fid,'(E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16,1X,E24.16)\n');
fprintf(fid,'*cfclos\n\n');

fprintf(fid,'*cfopen,nMasses,txt\n');
fprintf(fid,'*vwrite,nMass(1)\n');
fprintf(fid,'(F12.8)\n');
fprintf(fid,'*cfclos\n\n');

fprintf(fid,'allsel\n\n');
fprintf(fid,'finish\n');
fprintf(fid,'save\n');

fclose(fid);
end