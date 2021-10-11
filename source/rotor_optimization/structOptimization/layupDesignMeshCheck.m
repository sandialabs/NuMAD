function bladeMeshCheck(xlsFile,aeSizes,plotOnly)

if ~plotOnly
    numelem=aeSizes.^-2;
   
    addpath('C:\Users\brresor\Documents\SNL_NuMAD_ToolBox_PreNuMAD')
    
    disp('Reading XLS input file...')
    blade = xlsBlade(xlsFile);
    disp('Creating bladeDef and writing NMD files...')
    blade.updateGeometry
    blade.updateKeypoints
    blade.updateBOM
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')
else
    load meshConvergenceResults
end

for i=1:length(mesh)
    if runANSYS
        blade.mesh=mesh(i);
        disp('Checking buckling criteria using ANSYS...')
        result=bladeDesignANSYS(blade);
        buckle(i,:)=result';
    end
    
    
    figure(1)
    set(1,'Name','Mesh convergence results')
    subplot(2,2,1)
    semilogx(mesh(1:i),buckle(1:i,:),'-o')
    ylabel('Buckling load factor')
    grid on
    subplot(2,2,2)
    semilogx(numelem(1:i),buckle(1:i,:),'-o')
    grid on
    legend('Mode 1','Mode 2','Mode 3','Mode 4','Mode 5')
    if i>1
        % percent change for a doubling of number of elements
        tmp1=diff(buckle)./buckle(1:end-1,:);
        tmp2=log(numelem(2:end)./numelem(1:end-1)) ./ log(2);
        tmp2=[ tmp2 tmp2 tmp2 tmp2 tmp2 ];
        perChange=tmp1./tmp2*100;
        subplot(2,2,3)
        semilogx(mesh(1:i-1),perChange(1:i-1,:),'-o')
        ylabel('% change due to doubling the elements #')
        xlabel('Global mesh size setting, m')
        grid on
        subplot(2,2,4)
        semilogx(numelem(1:i-1),perChange(1:i-1,:),'-o')
        ylabel('% change due to doubling the elements #')
        xlabel('Equivalent elements # per m^2')
        grid on
        
    end
    
    
end

save meshConvergenceResults numelem mesh buckle

end