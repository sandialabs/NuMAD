scaleFactor=1e6; %Used to plot stress in MPa and strain in microstrain
close all;
figureNumbers=[4,5];
myYlabel='$y_3$ [mm]';
plotLocalFields(figureNumbers,myYlabel,myresult.element1950,scaleFactor,'k')
plotLocalFields(figureNumbers,myYlabel,myresult.element2558,scaleFactor,'r')
legend('LE Reinf. Element','Spar Cap element')