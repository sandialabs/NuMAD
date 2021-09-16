function coeffs=polyfitmodes(nargin)
%polyfitmodes  Read BMODES inputs and results to fit polynomials to
%              calculated mode shapes
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   polyfitmodes
%
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%      input  = description of inputs
%  
%      output = description of outputs
%

%===== CREDITS & CHANGELOG ================================================
% 2011.06.09  brr: initial release
% yyyy.mm.dd  initials: description

% calculate mode shape polynomial fit coefficients
bmodes=readBModesMain('bmodes.bmi');
bmout=readBModesOut('bmodes.out',length(bmodes.el_loc));

for i=1:length(bmout.freq)
   one=max(abs(bmout.tab{i}(end,[2 4])));
   bmout.tab{i}(:,[2 4])=bmout.tab{i}(:,[2 4])/one;
end

j=1;
k=1;

figure('Name','BModes Calculated Mode Shape Inspection')

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
title('Flapwise')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
title('Edgewise')
j=j+1;
k=k+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
j=j+1;
k=k+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
j=j+1;
k=k+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
j=j+1;
k=k+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
j=j+1;
k=k+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,2),'-o')
ylabel(['Mode ' num2str(k)])
j=j+1;

subplot(6,2,j)
plot(bmout.tab{k}(:,1),bmout.tab{k}(:,4),'-o')
j=j+1;
k=k+1;

disp('Enter mode number for 1st flapwise bending')
flp1=input('  (refer to Figure 1) Press Enter to choose Mode 1  ');
disp('Enter mode number for 2nd flapwise bending')
flp2=input('  (refer to Figure 1) Press Enter to choose Mode 3  ');
disp('Enter mode number for 1st edgewise bending')
edg1=input('  (refer to Figure 1) Press Enter to choose Mode 2  ');

if isempty(flp1);flp1=1;end
if isempty(flp2);flp2=3;end
if isempty(edg1);edg1=2;end

spnpts=bmout.tab{1}(:,1);
flp1pts=bmout.tab{flp1}(:,2);
flp2pts=bmout.tab{flp2}(:,2);
edg1pts=bmout.tab{edg1}(:,4);

if flp1pts(end)<0,flp1pts=flp1pts*-1;end
if flp2pts(end)<0,flp2pts=flp2pts*-1;end
if edg1pts(end)<0,edg1pts=edg1pts*-1;end

% fit the data with polynomial
A=[spnpts.^2 spnpts.^3 spnpts.^4 spnpts.^5 spnpts.^6];
[U,S,V]=svd(A);
x=(0:0.05:1);

figure('Name','Chosen Modes')

B=flp1pts;
p=(U*S*V')\B;
pflp1=[p(end:-1:1)' 0 0];
pflp1=pflp1/sum(pflp1);
y=polyval(pflp1,x);
subplot(1,3,1)
plot(spnpts,flp1pts/flp1pts(end),'-o',x,y,'k-d')
legend('Modal displacements','Polynomial Fit','Location','Northwest')
xlabel('Flap1')

B=flp2pts;
p=(U*S*V')\B;
pflp2=[p(end:-1:1)' 0 0];
pflp2=pflp2/sum(pflp2);
y=polyval(pflp2,x);
subplot(1,3,2)
plot(spnpts,flp2pts/flp2pts(end),'-o',x,y,'k-d')
xlabel('Flap2')
title('Chosen Modes')

B=edg1pts;
p=(U*S*V')\B;
pedg1=[p(end:-1:1)' 0 0];
pedg1=pedg1/sum(pedg1);
y=polyval(pedg1,x);
subplot(1,3,3)
plot(spnpts,edg1pts/edg1pts(end),'-o',x,y,'k-d')
xlabel('Edge1')

blade.BldFl1Sh=[0 pflp1(5:-1:1)];  % ignore the x^0 term
blade.BldFl2Sh=[0 pflp2(5:-1:1)];  % ignore the x^0 term
blade.BldEdgSh=[0 pedg1(5:-1:1)];  % ignore the x^0 term

save polyfitdata blade

coeffs=blade;

end

