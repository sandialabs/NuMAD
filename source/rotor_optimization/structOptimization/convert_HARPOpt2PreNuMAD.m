load('blade.mat')

figure(1)
grid on
subplot(3,1,1)
plot(blade.span,blade.chord,'b-o')
grid on
subplot(3,1,2)
plot(blade.span,blade.degreestwist,'b-o')
grid on
subplot(3,1,3)
plot(blade.span,blade.percentthick,'b-o')
grid on

[~,maxc]=max(blade.chord);

data=[blade.span(maxc) blade.degreestwist(maxc) blade.chord(maxc) NaN];
% data=[data; blade.span(end) blade.degreestwist(end) blade.chord(end) NaN];

for i = 1:length(blade.stations)
    data=[data ; blade.stations(i).spanlocation, blade.stations(i).degreestwist, blade.stations(i).chord, blade.stations(i).percentthick];
end

data=sortrows(data,1);

fid=fopen('HARPOpt2PreNuMAD.csv','w+');
fprintf(fid,'%f, %f, %f, %f\n',data');
fclose(fid);

disp(blade.components(5).name)
blade.components(5).cp
disp('sparcapwidth')
blade.sparcapwidth
disp('teband')
blade.teband
disp('leband')
blade.leband
