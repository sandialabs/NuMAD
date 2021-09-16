function eb131readfeaNew(namefea)
% matlab file eb131readfea
%  reads all output from Numad/Ansys run for use in Equivalent beam code, Ebeam
%  data = node id, initial coordinates, six sets of global displacements.
%  10  Oct  2002 (correction to transformation of displ due to Fz)
%  20 nov 2003, v27
%  14 april 2004, v30, outer profile added
%  22 july 2005, v31, reference to profile removed (solid elements)
%  02 sept 2009, v131, brr, cleaned up for NuMAD usage
%
global nbeamnode nfeanode zbeamnode feanodeid coord displ precurve feamass feaflag givenkflag wtflag expflag sizedispl

% temp = loadAscii(namefea); % opens file with name of basic data
temp = dlmread(namefea);

[nfeanode,ncfeanode] =size(temp);

feanodeid = temp(:,1);
coord = temp(:,2:4);

% I am going to reprocess the X Y and Z displacements so that they are X Y
% and Z for each run.
tDispl = temp(:,5:ncfeanode-1);
displ = zeros(size(tDispl));
nRuns = size(tDispl,2)/3;
for iR = 1:nRuns
    displ(:,(1:3)+(iR-1)*3) = ...
        [tDispl(:,iR) tDispl(:,nRuns+iR) tDispl(:,2*nRuns+iR)];
end

feamass = temp(:,ncfeanode);
temp=[];

% transform to IEC coordinates
tmp = coord(:,1);
coord(:,1) = coord(:,2);
coord(:,2) = -tmp;

tmp = displ(:,1:3:end);
displ(:,1:3:end) = displ(:,2:3:end);
displ(:,2:3:end) = -tmp;

[nrdispl,ncdispl] = size(displ);
disp(['Size of nodal displacement matrix: ' num2str(nrdispl) 'x' num2str(ncdispl) ' (rxc)']);

return



