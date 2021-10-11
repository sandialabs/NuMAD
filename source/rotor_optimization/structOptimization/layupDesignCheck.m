function layupDesignCheck(xlsFile)

mesh=0.10;
addpath('C:\Users\brresor\Documents\SNL_NuMAD_ToolBox_PreNuMAD')

disp('Reading XLS input file...')
blade = xlsBlade(xlsFile);
disp('Creating bladeDef and writing NMD files...')
blade.updateGeometry
blade.updateKeypoints
blade.updateBOM
blade.mesh=mesh;
BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')

disp('Checking OoP deflection using FAST...')
designvar1=layupDesignFAST(blade);
disp('Checking buckling criteria and blade mass using ANSYS...')
[designvar2,mass]=layupDesignANSYS(blade);
copyfile('file.rst','file.rst.buckle');
% % % % disp('Computing blade mass using ANSYS...')
% % % % mass=layupDesignANSYSmass;
disp('Computing blade frequencies using PreComp and BModes...')
[~,freqs] = layupDesignFASTBlade(blade)

disp(sprintf('RESULTS: Buckle %2.4f, OoPDefl %2.4f, Mass %6.0fkg',designvar2(1),designvar1,mass))

fid=fopen('layupDesignCheck.txt','w+');
fprintf(fid,'RESULTS: Buckle %2.4f, OoPDefl %2.4f, Mass %6.0fkg\n\n',designvar2(1),designvar1,mass);
fprintf(fid,'Frequencies:\n');
fprintf(fid,'  %f\n',freqs);
fclose(fid);

end