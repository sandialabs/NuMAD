function out = MFatigue(wt,rccdata,cs,EIs,params)
%% MFATIGUE  Calculates fatigue damage using rainflow cycle count data
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

Yr=params.designLife;
simtime=params.simtime;
N_save = [];
%% For each material and each channel, find the total cumulative fatigue
%  damage due to the entire range of wind speeds

% TODO remove material properties from the DLCDef. Import from blade object
% and/or GUI.
matData = params.matData;
for mat=1:length(matData)
    E=matData(mat).ex; % (Pa)
    m=matData(mat).m; % fatigue slope exponent
    uts=matData(mat).uts;  % assumes ultimate strength format S = [UCS UTS] (Pa)
    ucs=matData(mat).ucs;    % ultimate compressive strength (equal to UTS if S = [UTS])
    %ble: ultimate strains could be added in material properties or
    %calculated as done below
    eps_uts = uts/E;
    eps_ucs = ucs/E;
    % material reduction factors are applied to the characteristic strength
    gamma_ms = matData(mat).gamma_ms;
    gamma_mf = matData(mat).gamma_mf;
    FSloads = 1.35; % from IEC design standard for DLC 1.2
    
    for ch=1:size(rccdata,1)                
        % reset D for each new channel
        % D = sum( ni/Nfi, i = 1:total number of cycles)
        Dminer=0;
        
        % determine the blade span location of the channel
        if contains(rccdata{ch,1}.label,'Root')
            chSpan = 1; % root location along blade span        
        elseif contains(rccdata{ch,1}.label,'Spn')
            chSpan = str2num(rccdata{ch,1}.label(4))+1;
        else 
            error('MFatigue does not recognize the variable channel in rccdata')
        end
        
        for w=1:size(rccdata,2)
            % put data from rccdata structure for this channel and wind speed into a temporary variable, data
            data=rccdata{ch,w};
            % save the ultimate strain for the cycle respecting magnitude
            mean_negative = data.means<0;
            eps_u = repmat(eps_uts,size(data.means));
            eps_u(mean_negative) = eps_ucs;
            
            % determine the scaling for bending moments (ignore forces)
            if contains(data.label,'M') && contains(data.label,'x')
                % EI and neutral axis location, c, for edgewise moment (Mxb in FAST)
                cy = cs(chSpan,1);
                EIedge = EIs(chSpan,1);
                c=cy; EI=EIedge;
            elseif contains(data.label,'M') && contains(data.label,'y')
                % EI and neutral axis location, c, for flapwise moment (Myb in FAST)
                cx = cs(chSpan,2);
                EIflap = EIs(chSpan,2);
                c=cx;EI=EIflap;
            else % forces, torsional moments, and some calculated channels
                c=0; % don't report fatigue damage, not currently calculated
                EI=1;
            end            
                        
            % Find total fatigue damage due to operation over entire wind
            %  speed range
            for jj=1:length(data.cycles)
                
                if contains(data.label,'eps')
                    % calculate failure based on the strain state (combined loads)
                    eps_a = data.amplitudes(jj);
                    eps_m = data.means(jj);
                    
                else % calculate failure based on the stress state
                    M_a = data.amplitudes(jj);
                    M_m = data.means(jj);
                    
                    % Calculate strain values based on cyclic moments from simulation
                    %  Assumptions:
                    %    om=stress=M*c/I
                    %    eps=strain=M*c/(E*I)
                    eps_a = M_a*c/EI;
                    eps_m = M_m*c/EI;
                    
                    sigma_a = M_a*c/EI * E;
                    sigma_m = M_m*c/EI * E;
                    
                    % check data by plotting the strain levels
                    if contains(data.label,'M') && contains(data.label,'y')
                        eps_My_a(jj,1) = eps_a;
                        eps_My_m(jj,1) = eps_m;
                    end
                end
                
                
                
                % Determine the maximum number of cycles for failure based
                % on available fatigue failure criterion or from data
                switch params.fatigueCriterion
                    case 'Goodman'
                        % Use load amplitude only or use the equivalent load
                        switch params.fatigueStress
                            case 'Amplitude Only'
                                eps = eps_a;
                                sigma = sigma_a;
                            case 'Equivalent'
                                % eps_u is the ultimate strain of the material
                                eps = eps_a/(1-abs(eps_m/eps_u(jj)));
%                                 sigma = sigma_a/(1-sigma_m/sigma_u);
                        end
                        % Goodman line for failure
                        Nf=((eps*gamma_mf*FSloads)/abs(eps_u(jj)))^-m;
                    case 'Shifted Goodman'
                        % add later
                        Nf=((eps_uts+abs(eps_ucs)-abs(2*eps_m*gamma_ms*FSloads-eps_uts+abs(eps_ucs)))/(2*eps_a*gamma_mf*FSloads))^m;
                    case 'Fatigue Data'
                        % add 5 or 7 point fatigue data and calculations
                end
                
                % error checks
                if Nf < 0
                    Nf=1;
                    warning('check calculations for failure cycles')
                end
                
                % use Miner's rule for fatigue damage calculation 
                % NOTE: failure assumed for D=1
                n=data.cycles(jj)/simtime * (60*60*24*365.24*Yr) * wt(w);  % cycles per second times number of seconds in Yr years times Rayleigh weight factor
                Dminer=Dminer+n/Nf;
            end
            
% %             if 0%contains(rccdata{ch,1}.label,'eps') && w == 4 && mat == 1
% %                 disp(mat)
% %                 disp(w)
% %                 pause(0.2)
% %                 figure
% %                 plot(data.means.*1e6, data.amplitudes.*1e6, 'x'); grid on
% %                 title(['Channel: ' data.label ', wind speed = ' num2str(data.windspeed) ' m/s'])
% %                 ylabel('mean strain')
% %                 xlabel('alternating strain')
% %             end
% %             
% %             if 0%contains(data.label,'M') && contains(data.label,'y') && w == 4 && mat == 1              
% %                 figure
% %                 plot(eps_My_m.*1e6, eps_My_a.*1e6, 'rx'); grid on
% %                 title(['Channel: ' data.label ', wind speed = ' num2str(data.windspeed) ' m/s'])
% %                 ylabel('mean strain')
% %                 xlabel('alternating strain')
% %             end
        end
        % Display results on screen
%         fprintf('%s (E=%3.1fGPa,UTS=%3.0fMPa,UCS=%3.0fMPa,b=%3.1f), %s (c=%5.4fm,EI=%4.2fGPa) %4.3g\n',matData(mat).Name,E/1e9,uts/1e6,ucs/1e6,m,data.label,cs(chSpan),EIs(chSpan)/1e9,Dminer)
        out(ch,mat)=Dminer;
    end
end

xlsData={''};
for i=1:length(matData)  % step across columns for various materials
    xlsData=[xlsData, {matData(i).name}];
end
for ch=1:size(rccdata,1) % step down rows for all the computed gage locations
    tmp={rccdata{ch,1}.label};
    for mat=1:length(matData)
        tmp=[tmp out(ch,mat)];
    end
    xlsData=[xlsData; tmp];
end

csvfn = 'IECDLC_1p2_F.csv';
if exist(csvfn, 'file')
    delete(csvfn);
end

try % save the file as a .csv to prevent errors with saving in Excel
    saveFmt1 = [repmat('%s,',1,length(matData)+1) '\n'];
    saveFmt2 = ['%s,' repmat('%.4e,',1,length(matData)) '\n'];
    fid = fopen(csvfn,'w');
    fprintf(fid,saveFmt1,xlsData{1,:});    
    for ii = 2:size(xlsData,1)
        fprintf(fid,saveFmt2,xlsData{ii,:});
    end
    fclose(fid);
catch
    error('file save issue')
end

end





