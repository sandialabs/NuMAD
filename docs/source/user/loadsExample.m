%Get extremum loads
fast_gage=get_fast_gage(params.momentMaxRotation);
output=layupDesign_FASTanalysis(blade,'1.3',0,0);

% Loads to be applied to blade
halfdz=2.5; %for best results halfdz should be a multiple of L and should be <= L/2

bladeLength = blade.ispan(end);
r=2.5:5:bladeLength; %Location where output forces and moments will be calculated

loadsTable = FastLoads4ansys(output,fast_gage,r);