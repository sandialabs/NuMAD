function output = sortExtremes(params,output,CaseName,outName,cs,EIs,EA_normal)

i=0;

% base subset of loads analysis (when params.fullLoads = 0)
i=i+1;
str(i).chanLabels={'RootMyb1','RootMyb2','RootMyb3'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MaxBladeFlapBendingMoment';
str(i).minOrMax='max';

i=i+1;
str(i).chanLabels={'RootMyb1','RootMyb2','RootMyb3'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinBladeFlapBendingMoment';
str(i).minOrMax='min';

% i=i+1;
% str(i).chanLabels={'RootMyc1','RootMyc2','RootMyc3'};
% str(i).gains=ones(length(str(i).chanLabels),1);
% str(i).fieldName='MaxRotorPlaneFlapBendingMoment';
% str(i).minOrMax='max';
% 
% i=i+1;
% str(i).chanLabels={'RootMyc1','RootMyc2','RootMyc3'};
% str(i).gains=ones(length(str(i).chanLabels),1);
% str(i).fieldName='MinRotorPlaneFlapBendingMoment';
% str(i).minOrMax='min';

% <<<<<<<<<<<<<<<<<<<<<< Blade Strain 1 <<<<<<<<<<<<<<<<<<<<<<
i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9}];
% calculate c/EI for each gage location
str(i).gains=cs(:,2)./EIs(:,2)*1e6*1000;  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MaxFlapStrain';
str(i).minOrMax='max';
i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9}];
% calculate c/EI for each gage location
str(i).gains=cs(:,2)./EIs(:,2)*1e6*1000;  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MinFlapStrain';
str(i).minOrMax='min';

i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMxb1'} params.bladeGageLabels_MLx{1:9}];
str(i).gains=cs(:,1)./EIs(:,1)*1e6*1000;  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MaxEdgeStrain';
str(i).minOrMax='max';
i=i+1; %warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMxb1'} params.bladeGageLabels_MLx{1:9}];
str(i).gains=cs(:,1)./EIs(:,1)*1e6*1000;  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MinEdgeStrain';
str(i).minOrMax='min';


% <<<<<<<<<<<<<<<<<<<<<< Blade Strain 1 <<<<<<<<<<<<<<<<<<<<<<
i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9} {'RootFzb1'} params.bladeGageLabels_FLz{1:9}];
% calculate c/EI for each gage location
str(i).gains=[cs(:,2)./EIs(:,2)*1e6*1000; 1000*1e6./EA_normal]';  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MaxHPStrain';
str(i).minOrMax='combined strain max';

i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9} {'RootFzb1'} params.bladeGageLabels_FLz{1:9}];
% calculate c/EI for each gage location
str(i).gains=[cs(:,2)./EIs(:,2)*1e6*1000 1000*1e6./EA_normal];  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MinHPStrain';
str(i).minOrMax='combined strain min';

i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9} {'RootFzb1'} params.bladeGageLabels_FLz{1:9}];
% calculate c/EI for each gage location
str(i).gains=[cs(:,2)./EIs(:,2)*1e6*1000; 1000*1e6./EA_normal];  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MaxLPStrain';
str(i).minOrMax='combined strain max';

i=i+1; warning('Only blade 1 is being calculated')
str(i).chanLabels=[{'RootMyb1'} params.bladeGageLabels_MLy{1:9} {'RootFzb1'} params.bladeGageLabels_FLz{1:9}];
% calculate c/EI for each gage location
str(i).gains=[cs(:,2)./EIs(:,2)*1e6*1000; 1000*1e6./EA_normal];  % convert from strain to microstrain and convert from kN to N
str(i).fieldName='MinLPStrain';
str(i).minOrMax='combined strain min';




% <<<<<<<<<<<<<<<<<<<<<< Blade Deflections <<<<<<<<<<<<<<<<<<<<<<
i=i+1;
str(i).chanLabels={'OoPDefl1','OoPDefl2','OoPDefl3'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MaxOoPDefl';
str(i).minOrMax='max';

i=i+1;
str(i).chanLabels={'OoPDefl1','OoPDefl2','OoPDefl3'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinOoPDefl';
str(i).minOrMax='min';

i=i+1;
str(i).chanLabels={'TipClrnc1','TipClrnc2','TipClrnc3'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinTipClrnc';
str(i).minOrMax='min';

i=i+1;
str(i).chanLabels={'LSShftFxs'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MaxRotorTh';
str(i).minOrMax='max';

i=i+1;
str(i).chanLabels={'LSShftFxs'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinRotorTh';
str(i).minOrMax='min';

% <<<<<<<<<<<<<<<<<<<<< Nacelle IMU Acceleration <<<<<<<<<<<<<<<<<<<<<
i=i+1;
str(i).chanLabels={'NcIMUTAxs'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MaxNacIMUAxs';
str(i).minOrMax='max';
i=i+1;
str(i).chanLabels={'NcIMUTAxs'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinNacIMUAxs';
str(i).minOrMax='min';

i=i+1;
str(i).chanLabels={'NcIMUTAys'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MaxNacIMUAys';
str(i).minOrMax='max';
i=i+1;
str(i).chanLabels={'NcIMUTAys'};
str(i).gains=ones(length(str(i).chanLabels),1);
str(i).fieldName='MinNacIMUAys';
str(i).minOrMax='min';


if params.momentMaxRotation < 180 % calculate max/min moments in rotated coordinates from the primary flap/edge directions
    
    % set the rotation angle, theta, for calculation of resultant moments.
    thetaMomentRotation = 0:params.momentMaxRotation:180; % [deg]
    thetaMomentRotation(thetaMomentRotation==180)=[]; % remove repetitive 180 deg rotation
    
    % Calculate the resultant moment maxima and minima for the root location.
    i=i+1;
    str(i).chanLabels={'RootMxb1' 'RootMyb1'}; % Mx must be listed first, then My
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRes_RootMrb1';
    str(i).minOrMax='resultant max';
    i=i+1;
    str(i).chanLabels={'RootMxb1' 'RootMyb1'}; % Mx must be listed first, then My
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRes_RootMrb1';
    str(i).minOrMax='resultant min';
    % Calculate the resultant moment maxima and minima for the strain gage locations.
    for ii = 1:length(params.bladeGageLabels_MLx)
        i=i+1;
        str(i).chanLabels=[params.bladeGageLabels_MLx(ii) params.bladeGageLabels_MLy(ii)]; % Mx must be listed first, then My
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['MaxRes_' params.bladeGageLabels_MLx{ii}(1:6) 'r' params.bladeGageLabels_MLx{ii}(end-1:end)],'L','');
        str(i).minOrMax='resultant max';
        i=i+1;
        str(i).chanLabels=[params.bladeGageLabels_MLx(ii) params.bladeGageLabels_MLy(ii)]; % Mx must be listed first, then My
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['MinRes_' params.bladeGageLabels_MLx{ii}(1:6) 'r' params.bladeGageLabels_MLx{ii}(end-1:end)],'L','');
        str(i).minOrMax='resultant min';
    end 
else
    thetaMomentRotation = nan;
    
    for ii = 1:length(params.bladeGageLabels_MLx)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLx{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_MLx{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLx{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_MLx{ii}],'L','');
        str(i).minOrMax='min';
    end
    
    for ii = 1:length(params.bladeGageLabels_MLy)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_MLy{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_MLy{ii}],'L','');
        str(i).minOrMax='min';
    end    
end


%% ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
if params.fullLoads   % perform the full loads analysis
    %%
%     i=i+1;
%     str(i).chanLabels={'RootFxc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFxc1';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFxc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFxc1';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootFxc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFxc2';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFxc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFxc2';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootFxc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFxc3';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFxc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFxc3';
%     str(i).minOrMax='min';
%     
%     %%
%     i=i+1;
%     str(i).chanLabels={'RootFyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFyc1';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFyc1';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootFyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFyc2';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFyc2';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootFyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootFyc3';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootFyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootFyc3';
%     str(i).minOrMax='min';
    
    
    %% ble <<<<<<<<<<<<<<<<<<<<<<<<<<
%     i=i+1;
%     str(i).chanLabels={'RootFxc1' 'RootFyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootFxc1Fyc1';
%     str(i).minOrMax='amplitude';
%     i=i+1;
%     str(i).chanLabels={'RootFxc2' 'RootFyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootFxc2Fyc2';
%     str(i).minOrMax='amplitude';
%     i=i+1;
%     str(i).chanLabels={'RootFxc3' 'RootFyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootFxc3Fyc3';
%     str(i).minOrMax='amplitude';
        
%     '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<';
%     i=i+1;
%     str(i).chanLabels={'RootMxc1' 'RootMyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootMxc1Myc1';
%     str(i).minOrMax='amplitude';
%     i=i+1;
%     str(i).chanLabels={'RootMxc2' 'RootMyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootMxc2Myc2';
%     str(i).minOrMax='amplitude';
%     i=i+1;
%     str(i).chanLabels={'RootMxc3' 'RootMyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='AmplRootMxc3Myc3';
%     str(i).minOrMax='amplitude';
    
    '<<<<<<<<<<<<<<<<<< BLADE ROOT BEARING <<<<<<<<<<<<<<<<<<<<<<<<<';
    i=i+1;
    str(i).chanLabels={'RootMxb1' 'RootMyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootMxb1Myb1';
    str(i).minOrMax='amplitude';
    i=i+1;
    str(i).chanLabels={'RootMxb2' 'RootMyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootMxb2Myb2';
    str(i).minOrMax='amplitude';
    i=i+1;
    str(i).chanLabels={'RootMxb3' 'RootMyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootMxb3Myb3';
    str(i).minOrMax='amplitude';
    
    '<<<<<<<<<<<<<<< BLADE ROOT BEARING - SHEAR FORCES <<<<<<<<<<<<<<<';
    i=i+1;
    str(i).chanLabels={'RootFxb1' 'RootFyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootFxb1Fyb1';
    str(i).minOrMax='amplitude';
    i=i+1;
    str(i).chanLabels={'RootFxb2' 'RootFyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootFxb2Fyb2';
    str(i).minOrMax='amplitude';
    i=i+1;
    str(i).chanLabels={'RootFxb3' 'RootFyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplRootFxb3Fyb3';
    str(i).minOrMax='amplitude';
    
    '<<<<<<<<<<<<<<<<<< YAW BEARING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<';       
    i=i+1;
    str(i).chanLabels={'YawBrMxn' 'YawBrMyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplYawBrMxnMyn';
    str(i).minOrMax='amplitude';
    
    i=i+1;
    str(i).chanLabels={'YawBrFxn' 'YawBrFyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplYawBrFxnFyn';
    str(i).minOrMax='amplitude';    
    
    '>>>>>>>>>>>>>>>>>>>> LOW-SPEED SHAFT >>>>>>>>>>>>>>>>>>>>>>>>>>>';
    i=i+1;
    str(i).chanLabels={'LSShftFys' 'LSShftFzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplLSShftFysFzs';
    str(i).minOrMax='amplitude';
    
    i=i+1;
    str(i).chanLabels={'LSSTipMys' 'LSSTipMzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplLSSTipMysMzs';
    str(i).minOrMax='amplitude';
    
    '>>>>>>>>>>>>>>>>>>>>> FOUNDATION/TOWER BASE >>>>>>>>>>>>>>>>>>>>>>>>>>>';
    i=i+1;
    str(i).chanLabels={'TwrBsMxt' 'TwrBsMyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplTwrBsMxtMyt';
    str(i).minOrMax='amplitude';   
    
    i=i+1;
    str(i).chanLabels={'TwrBsFxt' 'TwrBsFyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='AmplTwrBsFxtFyt';
    str(i).minOrMax='amplitude';
    
    %% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    i=i+1;
    str(i).chanLabels={'RootFxb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFxb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFxb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFxb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFxb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFxb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFxb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFxb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFxb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFxb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFxb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFxb3';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RootFyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFyb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFyb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFyb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFyb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFyb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFyb3';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RootFzb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFzb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFzb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFzb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFzb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFzb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFzb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFzb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootFzb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootFzb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootFzb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootFzb3';
    str(i).minOrMax='min';
    
%     %%
%     i=i+1;
%     str(i).chanLabels={'RootMxc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMxc1';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMxc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMxc1';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootMxc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMxc2';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMxc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMxc2';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootMxc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMxc3';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMxc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMxc3';
%     str(i).minOrMax='min';
%     
%     %%
%     i=i+1;
%     str(i).chanLabels={'RootMyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMyc1';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMyc1'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMyc1';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootMyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMyc2';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMyc2'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMyc2';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'RootMyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxRootMyc3';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'RootMyc3'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinRootMyc3';
%     str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RootMxb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMxb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMxb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMxb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMxb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMxb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMxb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMxb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMxb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMxb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMxb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMxb3';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RootMyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMyb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMyb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMyb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMyb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMyb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMyb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMyb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMyb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMyb3';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RootMzb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMzb1';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMzb1'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMzb1';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMzb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMzb2';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMzb2'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMzb2';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'RootMzb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRootMzb3';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RootMzb3'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRootMzb3';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'YawBrFxn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrFxn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrFxn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrFxn';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'YawBrFyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrFyn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrFyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrFyn';
    str(i).minOrMax='min';
        
    i=i+1;
    str(i).chanLabels={'YawBrFzn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrFzn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrFzn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrFzn';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'YawBrFxp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrFxp';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrFxp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrFxp';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'YawBrFyp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrFyp';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrFyp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrFyp';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'YawBrMxn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrMxn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrMxn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrMxn';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'YawBrMyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrMyn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrMyn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrMyn';
    str(i).minOrMax='min';
        
    i=i+1;
    str(i).chanLabels={'YawBrMzn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrMzn';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrMzn'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrMzn';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'YawBrMxp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrMxp';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrMxp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrMxp';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'YawBrMyp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxYawBrMyp';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'YawBrMyp'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinYawBrMyp';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'RotThrust'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxRotThrust';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'RotThrust'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinRotThrust';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'LSShftTq'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSShftTq';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSShftTq'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSShftTq';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'LSShftFxs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSShftFxs';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSShftFxs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSShftFxs';
    str(i).minOrMax='min';

    
    %%
    i=i+1;
    str(i).chanLabels={'LSShftFys'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSShftFys';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSShftFys'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSShftFys';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'LSShftFzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSShftFzs';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSShftFzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSShftFzs';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'LSShftMxs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSShftMxs';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSShftMxs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSShftMxs';
    str(i).minOrMax='min';
    
    %%
%     i=i+1;
%     str(i).chanLabels={'LSShftFya'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxLSShftFya';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'LSShftFya'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinLSShftFya';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'LSShftFza'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxLSShftFza';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'LSShftFza'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinLSShftFza';
%     str(i).minOrMax='min';
%     i=i+1;
%     str(i).chanLabels={'LSSTipMya'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxLSSTipMya';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'LSSTipMya'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinLSSTipMya';
%     str(i).minOrMax='min';
%     
%     i=i+1;
%     str(i).chanLabels={'LSSTipMza'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MaxLSSTipMza';
%     str(i).minOrMax='max';
%     i=i+1;
%     str(i).chanLabels={'LSSTipMza'};
%     str(i).gains=ones(length(str(i).chanLabels),1);
%     str(i).fieldName='MinLSSTipMza';
%     str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'LSSTipMys'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSSTipMys';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSSTipMys'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSSTipMys';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'LSSTipMzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxLSSTipMzs';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'LSSTipMzs'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinLSSTipMzs';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'TwrBsFxt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsFxt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsFxt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsFxt';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'TwrBsFyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsFyt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsFyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsFyt';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'TwrBsFzt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsFzt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsFzt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsFzt';
    str(i).minOrMax='min';
    
    %%
    i=i+1;
    str(i).chanLabels={'TwrBsMxt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsMxt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsMxt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsMxt';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'TwrBsMyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsMyt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsMyt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsMyt';
    str(i).minOrMax='min';
    
    i=i+1;
    str(i).chanLabels={'TwrBsMzt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MaxTwrBsMzt';
    str(i).minOrMax='max';
    i=i+1;
    str(i).chanLabels={'TwrBsMzt'};
    str(i).gains=ones(length(str(i).chanLabels),1);
    str(i).fieldName='MinTwrBsMzt';
    str(i).minOrMax='min';
    
    %%    
    for ii = 1:length(params.bladeGageLabels_FLx)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLx{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_FLx{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLx{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_FLx{ii}],'L','');
        str(i).minOrMax='min';
    end
    
    for ii = 1:length(params.bladeGageLabels_FLy)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_FLy{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_FLy{ii}],'L','');
        str(i).minOrMax='min';
    end
    
    for ii = 1:length(params.bladeGageLabels_FLz)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLz{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_FLz{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLz{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_FLz{ii}],'L','');
        str(i).minOrMax='min';
    end
    
    for ii = 1:length(params.bladeGageLabels_FLx)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_FLx{ii} params.bladeGageLabels_FLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Ampl' params.bladeGageLabels_FLx{ii} params.bladeGageLabels_FLy{ii}(5:end)],'L','');
        str(i).minOrMax='amplitude';
    end
        
    '-------------------------------------------------------';
    
    for ii = 1:length(params.bladeGageLabels_MLz)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLz{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Max' params.bladeGageLabels_MLz{ii}],'L','');
        str(i).minOrMax='max';
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLz{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Min' params.bladeGageLabels_MLz{ii}],'L','');
        str(i).minOrMax='min';
    end
    
    for ii = 1:length(params.bladeGageLabels_MLx)
        i=i+1;
        str(i).chanLabels={params.bladeGageLabels_MLx{ii} params.bladeGageLabels_MLy{ii}};
        str(i).gains=ones(length(str(i).chanLabels),1);
        str(i).fieldName=strrep(['Ampl' params.bladeGageLabels_MLx{ii} params.bladeGageLabels_MLy{ii}(5:end)],'L','');
        str(i).minOrMax='amplitude';
    end
end

% ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

%%
if strcmp(CaseName,'IECDLC1p1NTM')
    % extrapolate to 50 year loads using other method
    output=findExtremes_1p1peaks(outName,str,thetaMomentRotation,params.bladeGageCoordinateRotation,CaseName,output,0);
else
    % identify the time series maxima from FAST output files
    output=findExtremes(outName,str,thetaMomentRotation,params.bladeGageCoordinateRotation,CaseName,output,0);
end

end