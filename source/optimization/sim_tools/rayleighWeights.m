function wts = rayleighWeights(wsvector,avgws)
%RAYLEIGHWEIGHTS  Calculate the fractions of PDF for an array of wind speed bins
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   wts = rayleighWeights(wsvector,avgws)
%   Detailed function description describing inputs, outputs, and other necessary information.
%
%      wsvector = vector of bin centers at which to compute weights
%      avgws = average of the Rayleigh distribution function
%  
%      wts = vector of Rayleigh PDF fractions
%

ws=wsvector;

% define bin edges based on wind speeds in rccdata
binwidth=ws(2)-ws(1);
binedges=[ws(1)-binwidth/2 ws+binwidth/2];
% define wind bins
for j=1:length(binedges)-1
    windbins(j,:)=[binedges(j) binedges(j+1)];
end
% Calculate weights for each bin according to Rayleigh distribution
sig=avgws/sqrt(pi/2);
% find PDF and CDF of Rayleigh distribution
pdf=windbins.*exp(-windbins.^2/(2*sig^2))/sig^2;
cdf=1-exp(-windbins.^2/(2*sig^2));
% calculate weights
wts=diff(cdf,1,2);
% show sum of weights (should be close to one, and not greater than one)
disp(' ');
disp(['Sum of the Rayleigh weights is ' num2str(sum(wts))]);

figure
plot(ws,ws.*exp(-ws.^2/(2*sig^2))/sig^2)
ylabel('PDF')
xlabel('Wind Speed')

end