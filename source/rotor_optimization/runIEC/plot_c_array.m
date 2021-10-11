% geezerplot

load carray

figure
subplot(1,2,1)
plot(carray(:,1),carray(:,6),'-o')
xlabel('Blade span (m)')
ylabel('c, Flapwise (m)')
subplot(1,2,2)
plot(carray(:,1),carray(:,5),'-o')
xlabel('Blade span (m)')
ylabel('c, Edgewise (m)')
