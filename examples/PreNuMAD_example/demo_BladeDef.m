clear blade
% blade = xlsBlade('cylinder.xlsx');   % surface_area = pi*0.92*23.934 % = 69.1756
% blade = xlsBlade('revolvesurf.xlsx');  % surface_area = 277.8651
% blade = xlsBlade('metablade.xlsx');
blade = xlsBlade('NREL5MW.xlsx');
% blade = xlsBlade('100m.xlsx');

%%
% blade.ispan = linspace(blade.span(1),blade.span(end),21);  % define output stations
% blade.ispan = blade.span;
blade.updateGeometry
blade.updateKeypoints
% fprintf('Surface area: %g\n',sum(blade.keyareas(:)))
%%
% figure(1);
% surf(blade); axis tight equal
% view(0,90);
% blade.plotregions;

%%
blade.updateBOM
blade.writeBOMxls('bom.xlsx')
% BladeDef_to_NuMADfile(blade,'BladeDef_test.nmd','MatDBsi.txt');
BladeDef_to_NuMADfile(blade,'NuMAD/BladeDef_test.nmd','NuMAD/MatDBsi.txt','NuMAD/airfoils/');

return

%%
figure(2); blade.plotbom;

%%
figure(3); clf; blade.components(3).plotcp
xlabel('Blade Fraction'); ylabel('# of layers');