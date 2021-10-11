function [feaData] = eb131readfea(namefea)
% matlab file eb131readfea
%  reads all output from Numad/Ansys run for use in Equivalent beam code, Ebeam
%  data = node id, initial coordinates, six sets of global displacements.
%  10  Oct  2002 (correction to transformation of displ due to Fz)
%  20 nov 2003, v27
%  14 april 2004, v30, outer profile added
%  22 july 2005, v31, reference to profile removed (solid elements)
%  02 sept 2009, v131, brr, cleaned up for NuMAD usage
%  28 june 2012, bco, descriptive variable names and more comments

temp = loadAscii(namefea); % opens file with name of basic data

[numFeaNode,numColsFeaNode] =size(temp);

feaNodeId = temp(:,1);              %node ID #
coord = temp(:,2:4);                %nodal coordinate
displ = temp(:,5:numColsFeaNode-1); %nodal displacements
feaMass = temp(:,numColsFeaNode);   %nodal mass
temp=[];

% transform to IEC coordinates
temp       = coord(:,1); %  initial coords of fe model %stores NUMAD 'X' in temp
coord(:,1) = coord(:,2); % sets 'X' to NUMAD 'Y'
coord(:,2) = temp;       % sets 'Y' to NUMAD 'X'
temp=[];
    %note, this transformation does not seem to be  consistent with IEC
    %coordinates it would seem that it should becoord(:,2) = -temp

%displ due to Fx and Fy
temp       = displ(:,1:3);
displ(:,1) = displ(:,5);
displ(:,2) = displ(:,4);
displ(:,3) = displ(:,6);
displ(:,4) = temp(:,2);
displ(:,5) = temp(:,1);
displ(:,6) = temp(:,3);
temp=[];

% due to Fz
temp       = displ(:,7);
displ(:,7) = displ(:,8);
displ(:,8) = temp;
%displ(:,9)=displ(:,9);
temp=[];

% due to Mx and My
temp        = displ(:,10:12);
displ(:,10) = -displ(:,14);
displ(:,11) = -displ(:,13);
displ(:,12) = -displ(:,15);
displ(:,13) = -temp(:,2);
displ(:,14) = -temp(:,1);
displ(:,15) = -temp(:,3);
temp=[];

% due to Mz
temp        = displ(:,16:18);
displ(:,16) = -temp(:,2);
displ(:,17) = -temp(:,1);
displ(:,18) = -temp(:,3);
temp=[];

[nrdispl,ncdispl] = size(displ);

disp(['Size of nodal displacement matrix: ' num2str(nrdispl) 'x' num2str(ncdispl) ' (rxc)']);

%store data in feaData struct   bco 06/20/12
feaData.nfeanode = numFeaNode;
feaData.feanodeid = feaNodeId;
feaData.coord = coord;
feaData.displ = displ;
feaData.feamass = feaMass;


end



