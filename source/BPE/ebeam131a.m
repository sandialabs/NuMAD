%function Ebeam131a
% BPE MAIN ROUTINE
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% Matlab code to extract equivalent beam stiffnesses from ansys output.
% Original code and concept by David Malcolm, DNV.
% Adapted here for use with a more modern version of NuMAD.
% 
% For full background see,
%       Malcolm, D. J. & Laird, D. L. "Extraction of Equivalent Beam 
%         Properties from Blade Models." Wind Energy, 2007, 10, 135-137.
% 
%  This routine reads in basic data and calls other routines;
%  Includes section stiffness routine and inertial properties routine
%  10 oct 2002: version 20: optional over-ride by specified element 
%       stiffness matrices
%  14 march 03, v 23, added input and output lines for precurve offset
%  8 may 2003, v 24, mass weighting added in solution for section 
%       displacements
%  24 june 2003, v25, additional output file (section k with inertia/
%       length), flag added for no /blank input file
%  latest: 11july 2003, mass mts of inertia moved to be about c of g
%  14 july 2003, v26,  mass-weighting option, exponential interpolation 
%       option added
%  16 july 2003, dislsize (no of terms in  section displ solution) variable 
%       added
%  19 aug 2003, sort.m altered to detect near-circular shape
%  20 nov 2003, v 27, specialized to be linked with Legacy version of Numad
%              all givenk features, exp interpolation, and # of z-displ 
%              terms removed or hardwired.
%              Basic input assumed in "BPE.txt".  FEA displacement assumed 
%              in "displacement.txt"
%  16 jan 2004, dm, v28, includes "singular" routine to expand singular 
%       elements until non-singular
%  3 feb 2004, dm, v29, includes "convert" routine to transform matrices to 
%       precurved configuration a lot of renaming of  arrays
%  13 april 2004, dm, v30, includes outer profile input and output
%  27 april 2004, dm, singular.m: new K put back into storestiffsym
%  30 april 2004, dm, ebeam30a.m, modified for all "old" inputs(title, 
%       feaname etc)
%  22 july 2005, dm, ebeam31a.m, solid element version, modifications to 
%       remove reference to profile
%  02 sept 2009, brr, ebeam131a.m, cleaned up and arranged for better 
%       output to NuMAD for FAST Blade input
%  07 dec 2010, brr, ebeam131a.m, changed structure of input files
%                       -eliminate need for Ebeam131.inp
%  25 june 2012, bco, clean up and eliminate globals
%

% list of some variables
%       nbeamnode       =  no of nodes in beam (incl start & end)
%       nfeanode        =  no of nodes in fea model
%       zbeamnode       =  spanwise coords of input beam nodes
%       feanodeid       =  id no of fea node given by ansys
%       coord                    =  x y z coordinates of all fea nodes
%       displ                    =  x y z displacements of all fea nodes (6 load cases)
%       sectionnodecount=  no of fea nodes at each beam node
%       sectionnodecoord=  x y z coords of each fea node at each beam node
%       sectionnodedispl=  x y z displacments of each fea node at each beam node, each of 6 loads
%       sectionmaxchord =  maximum chord (distance between any 2 points) at a section
%       sectionchordrot =  positive rotation of the chord at a section
%       sectioncofpx    =  local x coordinate of center of pressure at section
%       sectioncofpy    =  local y coordinate of center of pressure at section
%       sectionLE       =  distance of leading edge from local origin
%       sectiondispl    =  6 overall displacements at each beam node for each of 6 loads
%       beamstiff       =  all stiffness matrices
%       storestiff      =  stores all beam stiffness matrices
%       storesectstiff  =  stores all section stiffness matrices
%       storestiffsym   =  stores symmetric version of all element stiffness matrices

fid2 = fopen('bpe.txt','r'); % opens basic input file
if fid2==-1, error('BPE input file "bpe.txt" not found'), end
namefea='displacement.txt';
feaflag=1;
nbeamnode = str2num(scandm2(fid2)); %#ok<ST2NM> % no of beam nodes desired
string = scandm2(fid2); % spanwise coords of beam nodes

[zbeamnode nbeamnode] = sscanf(string,'%f,',nbeamnode); % puts node values into vector
nodeData.zbeamnode = zbeamnode;
nodeData.nbeamnode = nbeamnode;

string = scandm2(fid2); % precurve offsets
[precurve number] = sscanf(string,'%f,',nbeamnode); % puts precurve offsets into vector

wtflag = str2num(scandm2(fid2)); %read flag for mass weighting ("0" = No, "1" = Yes)

% expflag = str2num(scandm2(fid2)) %read flag for exponential interpolation ("0" = No, "1" = Yes)
expflag  = 0;

sizedispl = str2num(scandm2(fid2)); % read size of terms for displ solution.  Min = 6, max=13

convertflag = str2num(scandm2(fid2)); %read flag for conversion to local coordinates ("0" = No, "1" = Yes)

givenkflag = 0;

flagData.expflag = expflag;
flagData.convertFlag = convertflag;
flagData.wtflag = wtflag;
flagData.feaflag = feaflag;
flagData.givenkflag = givenkflag;

kData.givenk = [];
kData.npairk = [];
kData.zpairk = [];

xyData.xlocate = [];
xyData.ylocate = [];
xyData.xrotate = [];
xyData.yrotate = [];

status=fclose(fid2);

[feaData] = eb131readfea(namefea); % read ansys output file

[sectionData] = eb131sort(feaData,nodeData); % sorts ansys data

[sectionData] = eb131sectdispl(nodeData,sectionData,flagData,sizedispl); % calcs the overall displacement at each beam node

[storeData,singularData] = eb131stiff(nodeData,sectionData,flagData,kData); % calculates all beam stiffnesses

% These three lines need to be commented out if method 2 or 3 is used in the matrix inverse, see eb131stiff.m for calls for method 2 and 3
% if(singularData.nsingular)>0 % if any elements are singular, call singular
%     [storeData] = eb131singular(feaData,nodeData,sectionData,flagData,singularData,storeData,sizedispl);
% end

[storeData] = eb131sectstiff(nodeData,storeData); % calculates all section stiffnesses

[massData] = eb131mass(feaData,nodeData,sectionData,flagData,kData); % calculates all inertial section propeties

if convertflag==1
    [xyData] = eb131convert; % transforms to local coordinates
end

eb131out(nodeData,sectionData,flagData,storeData,massData,xyData,precurve); % outputs data for wtprep and xls and givenk wrt global or local coords

disp('**** BPE Analysis Finished ****')

status=fclose('all');