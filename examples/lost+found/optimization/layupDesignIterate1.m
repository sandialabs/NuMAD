function blade = bladeDesignIterate1(blade)

cp=5;  % spar cap
mode='OoPDefl';

switch mode
    case 'freq'
        criteria = -0.25*3*1.1;
    case 'fatigue'
        criteria = 1.0;
    case 'OoPDefl'
        criteria = 7.07*1.10;
end

disp(sprintf('Working on spar cap driven by %s...',mode))
disp('Creating bladeDef and writing NMD files...')

tmp1=blade.components(cp).cp(2,2);
tmp2=blade.components(cp).cp(3,2);

for param=tmp1:-1:0
    blade.components(cp).cp(2,2)=param;
    blade.components(cp).cp(3,2)=tmp2/tmp1*param;
    blade.components(cp).name
    blade.components(cp).cp
    
    % update bladeDef
    blade.updateGeometry
    blade.updateKeypoints
    blade.updateBOM
    BladeDef_to_NuMADfile(blade,'numad.nmd','MatDBsi.txt')
    
    switch mode
        case 'freq'
            disp('Checking frequencies using NuMAD, PreComp, and BModes...')
            [~,freqs]=layupDesignFASTBlade(blade);
            designvar=-1*freqs(1);
        case 'fatigue'
            disp('Checking fatigue values using FAST...')
            tmp=layupDesignFASTFatigue(blade);
            designvar=max(tmp(11:20,2));
        case 'OoPDefl'
            disp('Checking OoP deflection using FAST...')
            designvar=layupDesignFAST(blade)
    end
    
    figure(1)
    subplot(3,1,1)
    plot(param,designvar,'o')
    grid on
    hold on
    title( blade.components(cp).name )
    
    if designvar > criteria
        disp('Design criteria has been met.')
        int_param=interp1([prev_designvar designvar],[prev_param param],criteria);
        int_param=ceil(int_param)
        blade.components(cp).cp(2,2)=int_param;
        blade.components(cp).cp(3,2)=tmp2/tmp1*int_param;
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