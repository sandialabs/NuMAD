clear blade
blade = xlsBlade('NREL5MW.xlsx');
blade.ispan = linspace(blade.span(1),blade.span(end),51);

afdb = [blade.stations.airfoil];
afdb.resample(400,'cosine');  % higher point density for Plot3D output
switch 1
    case 1
        % attempt to make a uniform TE thickness
        chord = [blade.stations.chord];
        for k=3:numel(afdb)
            thickTE = afdb(k).thickness(end);
            afdb(k).adjustTE(0.02/chord(k)-thickTE);
            afdb(k).TEtype = 'flat';
        end
    case 2
        % give airfoils same percent TE thickness
        for k=3:numel(afdb)
            thickTE = afdb(k).thickness(end);
            afdb(k).adjustTE(0.02-thickTE);
            afdb(k).TEtype = 'flat';
        end
end

% write the Plot3D file
blade.updateGeometry
blade.writePlot3D('demo_Plot3D.p3d',[-.3, .3]);

% show the results
figure(1);
clf; plotp3d('demo_Plot3D.p3d');
axis equal

figure(2);
TEthickness = squeeze(blade.geometry(2    ,1:3,:)) - ...
              squeeze(blade.geometry(end-1,1:3,:));
TEthickness = sqrt(sum(TEthickness.^2,1));
span = squeeze(blade.geometry(2,3,:));
plot(span,TEthickness);