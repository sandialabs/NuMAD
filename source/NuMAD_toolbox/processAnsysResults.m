% NOT A FUNCTION
% This script template can be used to interface ANSYS results to 
%  computations performed in Matlab.  Then back into ANSYS for display of 
%  results on the original model.
%
% Basic steps:
%  1) Run analysis in ANSYS
%  2) Export results and model information from ANSYS to files that can be
%     read by functions readANSYS*.m from the SNL Wind Turbine Toolbox
%  3) Modify the USER INPUT SECTION of this script to perform the
%  calculations required for your analysis
%  4) Run this script.  Two files are saved:
%     - results.txt, containing numerical values for entry into an ANSYS
%     Element Table
%     - readem.mac, an ANSYS macro that will read results.txt and place the
%     data in the ANSYS Element Table named 'z_results'
%  5) In ANSYS, load the original model and run the macro readem.mac
%
% B.Resor 8/9/2011

%%%%%%%% USER INPUT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in ANSYS output data lists (Use functions readANSYS*.m)
s1=readANSYSStrains('Strains1.txt');
s2=readANSYSStrains('Strains2.txt');

% Define names of computed data columns
names={'dEPELX','dEPELY','dEPELZ','dEPELXY','dEPELYZ','dEPELXZ'};

% Perform computations
% 'd' is an array of computed results for reloading back into ANSYS
% Column 1 of 'd' is the list of element numbers
% Columns 2 thru 'outs'+1 are computed data
% 'outs' is calculated later automatically based on the size of d that has been computed
d(:,1)=s1(:,1);
d(:,2:7)=s2(:,2:7)-s1(:,2:7);

%%%%%%%% END USER INPUT SECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Error checking
for j=1:length(names)
   if length(names{j})>8
       error(sprintf('Names are not set up correctly: %s has more than 8 characters',names{j}))
   end
end

if size(d,2)-1~=length(names)
    error('Problem is not set up correctly: ''size(d,2)-1~=length(names)''')
end

% Save calculated results in txt file
n=size(d,1);
outs=size(d,2)-1;
str='%7.0f';
for j=1:outs
    str=[str ' %+10.4e'];
end
str=[str '\n'];

fid=fopen('results.txt','wt');
for j=1:n
    fprintf(fid,str,d(j,:));
end
fclose(fid);

% Generate macro (readem.mac) to read in the results file into ANSYS element table
fid=fopen('readem.mac','wt');
fprintf(fid,'esel,all\n');
fprintf(fid,'/POST1\n');
fprintf(fid,'*get,z_num_elem,elem,0,count\n');
fprintf(fid,'*dim,z_results,array,z_num_elem,%i\n',outs+1);
fprintf(fid,'*VREAD,z_results(1,1),results,txt,,JIK,%i,z_num_elem\n',outs+1);
fprintf(fid,'(F8.0,%iE13.4)\n',outs);
for j=1:outs
    fprintf(fid,'ETABLE,%s\n',names{j});
end
fprintf(fid,'*do,i,1,z_num_elem\n');
for j=1:outs
    fprintf(fid,'  detab,z_results(i,1),%s,z_results(i,%i)\n',names{j},j+1);
end
fprintf(fid,'*enddo\n');

fclose(fid);

disp('          processAnsysResults.m Success')

