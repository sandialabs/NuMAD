function bladeDesignIterate(xlsFile)

mesh=0.10;
addpath('C:\Users\brresor\Documents\SNL_NuMAD_ToolBox_PreNuMAD')

disp('Reading XLS input file...')
blade = xlsBlade(xlsFile);
disp('Creating bladeDef and writing NMD files...')

for param=49:15:10000
    %update layups
        cp=6;  % panels
        blade.components(cp).cp(3,2)=param;
        blade.components(cp).cp(4,2)=param;
        blade.components(cp).cp(5,2)=param*27/49;
        blade.components(cp).cp(6,2)=param*27/49;
        blade.components(cp).name
        blade.components(cp).cp
%     
%     cp=5;  % spar cap
%     blade.components(cp).cp(2,2)=param;
%     blade.components(cp).cp(3,2)=122/108*param;
%     blade.components(cp).name
%     blade.components(cp).cp
%     
%         cp=7;   % te-reinf
%         blade.components(cp).cp(2,2)=param;
%         blade.components(cp).cp(3,2)=param;
%         blade.components(cp).name
%         blade.components(cp).cp
    
    % update bladeDef
    blade.updateGeometry
    blade.updateKeypoints
    blade.updateBOM
    blade.mesh=mesh;
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')
    
%         disp('Checking frequencies using NuMAD, PreComp, and BModes...')
%         [~,freqs]=layupDesignFASTBlade(blade);
%         designvar=freqs(2)/freqs(1);
        
%     disp('Checking OoP deflection using FAST...')
%     designvar=layupDesignFAST(blade)

            disp('Checking buckling criteria using ANSYS...')
            designvar=layupDesignANSYS(blade);

%         disp('Computing blade mass using ANSYS...')
%         mass=layupDesignANSYSmass;

%         disp('Checking fatigue values using FAST...')
%     tmp=layupDesignFASTFatigue(blade);
%     designvar=max(tmp(11:20,2));
    
    figure(1)
    plot(param,designvar,'o')
    grid on
    hold on
    
    if designvar(1) > 1.62
        disp('Design criteria has been met.')
        break;
    end
    
end