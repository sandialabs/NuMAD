addpath(genpath('/Users/clkell/Documents/git_design_codes/'))
% Create blade in prenumad
blade = BladeDef;
blade.readYAML('BAR00.yaml')
surf(blade)
view(-90,90)
return
% Change blade object dimensions
% Example make both shear web cores 10% thicker
% blade.components(13).cp(:,2) = blade.components(13).cp(:,2).*1.1;
% blade.components(16).cp(:,2) = blade.components(16).cp(:,2).*1.1;
% blade.updateGeometry;
% blade.updateKeypoints;
% blade.updateBOM;

% Write modified blade object to yaml file
blade.writeYAML('BAR00_mod.yaml')