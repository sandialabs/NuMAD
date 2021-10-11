function [designvar,freqs] = layupDesign_FASTmodal(bmodesFrequencies)

% NOTE: this script could be used to generate rotor modes from the blade
% modes calculated with BModes


% disp('Creating FAST Blade file using NuMAD and PreComp...')
% delete bmodesFrequencies.mat
% numad('numad.nmd','precomp',[1 3 2])
designvar=bmodesFrequencies(2)/bmodesFrequencies(1);
freqs=bmodesFrequencies;

end

