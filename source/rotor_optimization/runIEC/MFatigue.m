function out = MFatigue(cs,EIs,matData,simtime,avgws,Yr,sf)
% MFATIGUE  Calculates fatigue damage using rainflow cycle count data
%                          Under Construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% out = MFatigue(cs,EIs,matData,simtime,avgws,Yr,sf)
%
% out = NxM matrix of Miner's fatigue damage ratios
% N=number of analysis locations, or virtual sensor channels
% M=number of materials to analyze
%
% Required Inputs:
% cs = 1xN vector of skin offset values from the neutral axis; order corresponds to label sequence used for rows in rccdata [m]
% EIs= 1xN vector of section EI values; order corresponds to label sequence used for rows in rccdata [Nm^2]
% matData = data structure of material properties
%     matData(i).Name, material name string
%     matData(i).E, Elastic Moduli, E [Pa]
%     matData(i).b, fatigue exponent, b [-]
%     matData(i).C, single cycle material tensile strength, C [Pa]
% simtime = scalar; the time represented by the rainflow count data found in rccdata, seconds; this is used to calculate a cycle accumulation rate
% avgws = scalar; average wind speed for Rayleigh wind distribution, m/s
% Yr = scalar; number of years for extrapolation of the simulated cycle rate, years
% sf = scalar; total safety factor for application to stresses
%
% DAMAGE MODEL - Miner's Rule
% Damage=sum( ni/Nf(psf_f*psf_m*Si) ) <= 1.0
% The model predicts fatigue failure when the damage exceeds 1.0
% ni=number of cycles at stress amplitude Si (peak-to-peak)
% psf_f and psf_m are partial factors of safety for loads and materials, respectively
% Nf=number of cycles to failure
% Nf is determined based on the material fatigue model that is implemented.  The size of matData determines the material model that is used by this function.  The function includes the following material model options:
% * Constant R Goodman Fit
% * 2-parameter Goodman Fit (under construction)
% * more models coming later…
%
% CONSTANT-R GOODMAN FIT
% R=min stress/max stress
% S=C*Nf^(-m)                Eq.(1)
% With C=So, or single cycle strength
% Normalize this equation using the single cycle strength, and rearrange
% Log(S/S_o)=log(C/S_o)-m*Log(Nf)
% Now it is easier to see how m is the slope of the straight line fit to S-N data with y-intercept at zero on a log-log plot of the above equation.
% Return to Eq.(1) and rearrange to calculate Nf,
% Nf=(1/C*S)^(-b)            Eq.(2)
% Where b=1/m
%
% 2-PARAMETER GOODMAN FIT (under construction)
% <text here>
%
% CALCULATION OF BENDING STRESS, S
% The turbine blade is assumed to be a simple beam with flap and edgewise stiffness provided in EIs.
% Axial strain resulting from application of an external moment is calculated with the following equation
% strain=M*c/EI
% M  = moment about the neutral axis
% c  = perpendicular distance to the neutral axis
% I  = the second moment of inertia about the neutral axis
% assume linear elastic material,
% S = strain*E
% Then,
% S=M*c/EI*E
%

% Matlab file 'rccdata.mat' is created by running the aeroelastic
% simulation scripts followed by Crunch analysis
load rccdata.mat

% determine the wind speeds that are saved in rccdata
for w=1:size(rccdata,2)
    ws(w)=rccdata{1,w}.windspeed;
end

% check to make sure wind speeds are spaced evenly
if ~isempty(find(diff(diff(ws))~=0))
    disp(num2str(ws))
    error('Your windspeeds are not spaced evenly');
end

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
wt=diff(cdf,1,2)
% show sum of weights (should be close to one, and not greater than one)
disp(' ');
disp(['Sum of the Rayleigh weights is ' num2str(sum(wt))]);
disp(' ');

% check to make sure that the number of weights is equal to the number of
% simulated wind speeds in rccdata
if length(wt)~=size(rccdata,2)
    error('The number of Rayleigh weights is not equal to the number of wind speeds contained in the rccdata.mat file');
end


% For each material and each channel, find the total cumulative fatigue
%  damage due to the entire range of wind speeds
for mat=1:length(matData)
    E=matData(mat).E;
    b=matData(mat).b;
    C=matData(mat).C;
    for ch=1:size(rccdata,1)
        % reset D for each new channel
        % D = sum( ni/Nfi, i = 1:total number of cycles)
        D=0;
        mxstr=0;
        % determine the blade span location of the channel
        if ~isempty(strfind(rccdata{ch,1}.label,'Root'))
            chSpan = 1; % root location along blade span        
        elseif ~isempty(strfind(rccdata{ch,1}.label,'Spn'))
            chSpan = str2num(rccdata{ch,1}.label(4));
        else 
            error('MFatigue does not recognize the variable channel in rccdata')
        end
        
        for w=1:size(rccdata,2)
            % put data from rccdata structure for this channel and wind speed into a temporary variable, data
            data=rccdata{ch,w};
            % Find total fatigue damage due to operation over entire wind
            %  speed range
            for j=1:length(data.cycles)
                M=data.amplitudes(j);
                c=cs(chSpan);
                EI=EIs(chSpan);
                % Calculate strain values based on cyclic moments from simulation
                %  Assumptions:
                %    om=stress=M*c/I
                %    ep=strain=M*c/(E*I)
                tmp=M*c/EI;
                if tmp>mxstr,mxstr=tmp;end
                S=M*c/EI*E;
                Nf=(1/C*S*sf)^-b;
                n=data.cycles(j)/simtime * (60*60*24*365.24*Yr) * wt(w);  % cycles per second times number of seconds in Yr years times Rayleigh weight factor
                D=D+n/Nf;
            end
        end
        % Display results on screen
%         disp(sprintf('%s (E=%3.1fGPa,C=%3.0fMPa,b=%3.1f), %s (c=%5.4fm,EI=%4.2fGPa) %4.3g',matData(mat).Name,E/1e9,C/1e6,b,data.label,cs(ch),EIs(ch)/1e9,D))
        out(ch,mat)=D;
    end
end

xlsfn='IECDLC_1p2_F.xlsx';
xlsData={''};
for i=1:length(matData)  % step across columns for various materials
    xlsData=[xlsData, {matData(i).Name}];
end
for ch=1:size(rccdata,1) % step down rows for all the computed gage locations
    tmp={rccdata{ch,1}.label};
    for mat=1:length(matData)
        tmp=[tmp out(ch,mat)];
    end
    xlsData=[xlsData; tmp];
end

try xlswrite(xlsfn,xlsData,1);
catch keyboard, end

end





