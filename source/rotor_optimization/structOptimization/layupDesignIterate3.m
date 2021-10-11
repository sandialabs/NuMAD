function blade = bladeDesignIterate3(blade)

cp=6;  % te panels
criteria = -1.62;

disp('Working on spar cap...')
disp('Creating bladeDef and writing NMD files...')

tmp1=blade.components(cp).cp(3,2);
tmp2=blade.components(cp).cp(5,2);

for param=tmp1:-1:20
    blade.components(cp).cp(3,2)=param;
    blade.components(cp).cp(4,2)=param;
    blade.components(cp).cp(5,2)=param*tmp2/tmp1;
    blade.components(cp).cp(6,2)=param*tmp2/tmp1;
    blade.components(cp).name
    blade.components(cp).cp
    
    % update bladeDef
    blade.updateGeometry
    blade.updateKeypoints
    blade.updateBOM
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')
    
    disp('Checking buckling criteria using ANSYS...')
    [designvar,mass]=-1*layupDesignANSYS(blade);
    
    figure(1)
    subplot(3,1,3)
    plot(param,designvar,'o')
    grid on
    hold on
    title( blade.components(cp).name )
    
    designvar=designvar(1);
    
    if designvar > criteria
        disp('Design criteria has been met.')
        int_param=interp1([prev_designvar designvar],[prev_param param],criteria);
        int_param=ceil(int_param)
        blade.components(cp).cp(3,2)=int_param;
        blade.components(cp).cp(4,2)=int_param;
        blade.components(cp).cp(5,2)=int_param*tmp2/tmp1;
        blade.components(cp).cp(6,2)=int_param*tmp2/tmp1;
        blade.components(cp).name
        blade.components(cp).cp
        blade.updateGeometry
        blade.updateKeypoints
        blade.updateBOM
        return;
    end
    
    prev_designvar=designvar;
    prev_param=param;
    
end