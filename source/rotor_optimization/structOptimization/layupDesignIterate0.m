function layupDesignIterate0(xlsFile)

mesh=0.10;
addpath('C:\Users\brresor\Documents\SNL_NuMAD_ToolBox_PreNuMAD')

disp('Reading XLS input file...')
blade = xlsBlade(xlsFile);
blade.mesh=mesh;

% refine spar cap 
blade = layupDesignIterate1(blade);

% refine te reinf
% blade = layupDesignIterate2(blade);

% refine te panels
% blade = layupDesignIterate3(blade);

% save blade_refined blade

end