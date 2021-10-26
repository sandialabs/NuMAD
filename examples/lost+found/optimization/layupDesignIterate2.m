function blade = bladeDesignIterate2(blade)

cp=7;   % te-reinf
criteria=-1.3;

disp('Working on te reinf...')
disp('Creating bladeDef and writing NMD files...')

tmp1=blade.components(cp).cp(2,2);

for param=tmp1:-5:0
    blade.components(cp).cp(2,2)=param;
    blade.components(cp).cp(3,2)=param;
    blade.components(cp).name
    blade.components(cp).cp
    
    % update bladeDef
    blade.updateGeometry
    blade.updateKeypoints
    blade.updateBOM
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')
    
    disp('Checking frequencies using NuMAD, PreComp, and BModes...')
    [~,freqs]=layupDesignFASTBlade(blade);
    designvar=-1*freqs(2)/freqs(1);
    
    %         disp('Checking fatigue values using FAST...')
    %     tmp=layupDesignFASTFatigue(blade);
    %     designvar=max(tmp(1:10,1));   % <-- double check this to make    
    
    figure(1)
    subplot(3,1,2)
    plot(param,designvar,'o')
    grid on
    hold on
    title( blade.components(cp).name )
    
    if designvar > criteria
        disp('Design criteria has been met.')
        int_param=interp1([prev_designvar designvar],[prev_param param],criteria);
        int_param=ceil(int_param)
        blade.components(cp).cp(2,2)=int_param;
        blade.components(cp).cp(3,2)=int_param;
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