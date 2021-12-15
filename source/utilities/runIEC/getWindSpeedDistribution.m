function [wt,rccdata]=getWindSpeedDistribution(avgws)
    load('rccdata.mat','rccdata')
    % determine the wind speeds that are saved in rccdata
    for w=1:size(rccdata,2)
        ws(w)=rccdata{1,w}.windspeed;
    end

    % check to make sure wind speeds are spaced evenly
    if std(diff(ws))~=0
        disp(num2str(ws))
        error('Your windspeeds are not spaced evenly');
    end

    % define bin edges based on wind speeds in rccdata
    binwidth=ws(2)-ws(1);
    binedges=[ws(1)-binwidth/2 ws+binwidth/2];
    % define wind bins
    for jj=1:length(binedges)-1
        windbins(jj,:)=[binedges(jj) binedges(jj+1)];
    end
    % Calculate weights for each bin according to Rayleigh distribution
    sig=avgws/sqrt(pi/2);
    % find PDF and CDF of Rayleigh distribution
%     pdf=windbins.*exp(-windbins.^2/(2*sig^2))/sig^2;
    cdf=1-exp(-windbins.^2/(2*sig^2)); 
    % calculate weights
    wt=diff(cdf,1,2);
    % show sum of weights (should be close to one, and not greater than one)
    disp(' ');
    disp(['Sum of the Rayleigh weights is ' num2str(sum(wt))]);
    disp(' ');

    % check to make sure that the number of weights is equal to the number of
    % simulated wind speeds in rccdata
    if length(wt)~=size(rccdata,2)
        error('The number of Rayleigh weights is not equal to the number of wind speeds contained in the rccdata.mat file');
    end