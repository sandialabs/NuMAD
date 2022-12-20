%% Add paths to NuMAD source
cd('..\..');
addNumadPaths;

%% Return to the example directory where the needed files are located
cd('examples\2022bladeWorkshop');

%% Create a blade object as an instance of BladeDef and read in data from the .yaml file

blade = BladeDef;
fileName = 'myBlade.yaml';
blade.readYAML(fileName);

%% Set the global element size for the shell mesh

blade.mesh = 0.2;

%% Generate the shell mesh

adhes = 1;
[nodes,elements,OSSets,OSNodes,SWSets,adNds,adEls] = blade.getShellMesh(adhes);

%% Print first 10 nodes coordinates

nodes(1:10,:)

%% Print first 10 element connectivities

elements(1:10,:)

%% Print first 10 elements in section 4, spanwise segment 10 (pressure side spar cap)

OSSets(4,10).elementList(1:10)

%% Print first 10 elements in spanwise segment 10 of first shear web

SWSets{1}(10).elementList(1:10)


