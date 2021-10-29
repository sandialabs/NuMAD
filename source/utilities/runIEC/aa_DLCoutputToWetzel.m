clear all
clc


fn_inner = 'IECDLC_Results_A12S0_10kgInstrumentation_noTipMass_gageSet1inner_0204.xlsx';
fn_outer = 'IECDLC_Results_A12S0_10kgInstrumentation_noTipMass_gageSet2outer_0204.xlsx';

[innerGage] = readDLCoutput(fn_inner);
[outerGage] = readDLCoutput(fn_outer);

% Calculate the number of gage stations for the inner set
condition = 0; ctr = 1;
while condition == 0;
    if ~isfield(innerGage,['MaxRes_Spn' num2str(ctr) 'MLrb1_0deg'])
       condition = 1; innerGageStations = ctr-1;
    end
    ctr = ctr+1;
end
% Calculate the number of gage stations for the outer set
condition = 0; ctr = 1;
while condition == 0;
    if ~isfield(outerGage,['MaxRes_Spn' num2str(ctr) 'MLrb1_0deg'])
       condition = 1; outerGageStations = ctr-1;
    end
    ctr = ctr+1;
end


%% ************************************************************************
% This section saves the Wetzel design variables only (Fx,Fy,Fz,Mx,My,Mz) 
% and (Fxy, Mxy) in one structure
% *************************************************************************

GageStations = {};

% Save the Minimum Resultant Forces in a new structure to combine the inner
% and outer gage station locations

% ----------- Root Forces/Moments ------------ %
rootNamesMax = {'MaxRootFxb1'; 'MaxRootFyb1'; 'MaxRootFzb1'; 'MaxRootMxb1';...
    'MaxRootMyb1'; 'MaxRootMzb1'};
for ii = 1:length(rootNamesMax)
    GageStations.(rootNamesMax{ii}) = innerGage.(rootNamesMax{ii});
end
rootNamesMin = {'MinRootFxb1'; 'MinRootFyb1'; 'MinRootFzb1'; 'MinRootMxb1';...
    'MinRootMyb1'; 'MinRootMzb1'};
for ii = 1:length(rootNamesMin)
    GageStations.(rootNamesMin{ii}) = innerGage.(rootNamesMin{ii});
end
rootNamesAmpl = {'AmplRootFxb1Fyb1'; 'AmplRootMxb1Myb1'};
for ii = 1:length(rootNamesAmpl)
    GageStations.(rootNamesAmpl{ii}) = innerGage.(rootNamesAmpl{ii});
end
% Save the Resultant Moments in a new structure to combine the inner and 
% outer gage station locations
% ------------ Mr -------------- %
strVar = {}; orig_strVar = {};
for jj = 0:15:165
    strVar{end+1} = ['MaxRes_RootMrb1_' num2str(jj) 'deg'];
    GageStations.(strVar{end}) = innerGage.(strVar{end});
end
for jj = 0:15:165
    strVar{end+1} = ['MinRes_RootMrb1_' num2str(jj) 'deg'];
    GageStations.(strVar{end}) = innerGage.(strVar{end});
end

% ------------ Fx -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Max and minimum Force x - gage set 1
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLxb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLxb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Max and minimum Force x - gage set 2        
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLxb1'];
        orig_strVar{end+1} = ['MaxSpn' num2str(ii-innerGageStations) 'FLxb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLxb1'];
        orig_strVar{end+1} = ['MinSpn' num2str(ii-innerGageStations) 'FLxb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end
% ------------ Fy -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Max and minimum Force y - gage set 1
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLyb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLyb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Max and minimum Force y - gage set 2        
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLyb1'];
        orig_strVar{end+1} = ['MaxSpn' num2str(ii-innerGageStations) 'FLyb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLyb1'];
        orig_strVar{end+1} = ['MinSpn' num2str(ii-innerGageStations) 'FLyb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end
% ------------ Fz -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Max and minimum Force z - gage set 1
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLzb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLzb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Max and minimum Force z - gage set 2        
        strVar{end+1} = ['MaxSpn' num2str(ii) 'FLzb1'];
        orig_strVar{end+1} = ['MaxSpn' num2str(ii-innerGageStations) 'FLzb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'FLzb1'];
        orig_strVar{end+1} = ['MinSpn' num2str(ii-innerGageStations) 'FLzb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end
% ------------ Ampl:Fxy -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Amplitude of Force xy - gage set 1
        strVar{end+1} = ['AmplSpn' num2str(ii) 'FLxb1FLyb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Amplitude of Force xy - gage set 2
        strVar{end+1} = ['AmplSpn' num2str(ii) 'FLxb1FLyb1'];
        orig_strVar{end+1} = ['AmplSpn' num2str(ii-innerGageStations) 'FLxb1FLyb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end

% Save the Resultant Moments in a new structure to combine the inner and 
% outer gage station locations

% ------------ Mr -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        for jj = 0:15:165
            strVar{end+1} = ['MaxRes_Spn' num2str(ii) 'MLrb1_' num2str(jj) 'deg'];
            GageStations.(strVar{end}) = innerGage.(strVar{end});
        end
        for jj = 0:15:165
            strVar{end+1} = ['MinRes_Spn' num2str(ii) 'MLrb1_' num2str(jj) 'deg'];
            GageStations.(strVar{end}) = innerGage.(strVar{end});
        end
    else
        for jj = 0:15:165
            strVar{end+1} = ['MaxRes_Spn' num2str(ii) 'MLrb1_' num2str(jj) 'deg'];
            orig_strVar{end+1} = ['MaxRes_Spn' num2str(ii-innerGageStations) 'MLrb1_' num2str(jj) 'deg'];
            GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        end
        for jj = 0:15:165
            strVar{end+1} = ['MinRes_Spn' num2str(ii) 'MLrb1_' num2str(jj) 'deg'];
            orig_strVar{end+1} = ['MinRes_Spn' num2str(ii-innerGageStations) 'MLrb1_' num2str(jj) 'deg'];
            GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        end
    end
end
% ------------ Mz -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Max and minimum Force z - gage set 1
        strVar{end+1} = ['MaxSpn' num2str(ii) 'MLzb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'MLzb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Max and minimum Force z - gage set 2        
        strVar{end+1} = ['MaxSpn' num2str(ii) 'MLzb1'];
        orig_strVar{end+1} = ['MaxSpn' num2str(ii-innerGageStations) 'MLzb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
        strVar{end+1} = ['MinSpn' num2str(ii) 'MLzb1'];
        orig_strVar{end+1} = ['MinSpn' num2str(ii-innerGageStations) 'MLzb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end
% ------------ Ampl:Mxy -------------- %
strVar = {}; orig_strVar = {};
for ii = 1:innerGageStations+outerGageStations
    if ii <= innerGageStations
        % Amplitude of Force xy - gage set 1
        strVar{end+1} = ['AmplSpn' num2str(ii) 'MLxb1MLyb1'];
        GageStations.(strVar{end}) = innerGage.(strVar{end});
    else
        % Amplitude of Force xy - gage set 2
        strVar{end+1} = ['AmplSpn' num2str(ii) 'MLxb1MLyb1'];
        orig_strVar{end+1} = ['AmplSpn' num2str(ii-innerGageStations) 'MLxb1MLyb1'];
        GageStations.(strVar{end}) = outerGage.(orig_strVar{end});
    end
end


%% ************************************************************************
% This section saves the Wetzel design variables only (Fx,Fy,Fz,Mx,My,Mz) 
% and (Fxy, Mxy) in one structure
% *************************************************************************

% perform for the saved gage station information
names = fieldnames(GageStations);
ExtremeValGage = {};
for ii = 1:length(names)
    if strcmp(names{ii}(1:3),'Min')
        % calculate the table minimum
        [val, ind] = min(GageStations.(names{ii}).data);
        ExtremeValGage.(names{ii}) = GageStations.(names{ii})(ind,:);
    elseif strcmp(names{ii}(1:3),'Max') || strcmp(names{ii}(1:4),'Ampl')
        % calculate the table maximum
        [val, ind] = max(GageStations.(names{ii}).data);
        ExtremeValGage.(names{ii}) = GageStations.(names{ii})(ind,:);
    end
end
writeDLCoutput(ExtremeValGage, 'ExtremeValues_GageStations.xlsx');


% perform for the inner gage stations only
names = fieldnames(innerGage);
innerGage_ExtremeValGage = {};
for ii = 1:length(names)
    if strcmp(names{ii}(1:3),'Min')
        % calculate the table minimum
        [val, ind] = min(innerGage.(names{ii}).data);
        innerGage_ExtremeValGage.(names{ii}) = innerGage.(names{ii})(ind,:);
    elseif strcmp(names{ii}(1:3),'Max') || strcmp(names{ii}(1:4),'Ampl')
        % calculate the table maximum
        [val, ind] = max(innerGage.(names{ii}).data);
        innerGage_ExtremeValGage.(names{ii}) = innerGage.(names{ii})(ind,:);
    end
end
% writeDLCoutput(innerGage_ExtremeValGage, 'ExtremeValues_InnerGage.xlsx');

% perform for the outer gage stations only
names = fieldnames(outerGage);
outerGage_ExtremeValGage = {};
for ii = 1:length(names)
    if strcmp(names{ii}(1:3),'Min')
        % calculate the table minimum
        [val, ind] = min(outerGage.(names{ii}).data);
        outerGage_ExtremeValGage.(names{ii}) = outerGage.(names{ii})(ind,:);
    elseif strcmp(names{ii}(1:3),'Max') || strcmp(names{ii}(1:4),'Ampl')
        % calculate the table maximum
        [val, ind] = max(outerGage.(names{ii}).data);
        outerGage_ExtremeValGage.(names{ii}) = outerGage.(names{ii})(ind,:);
    end
end
% writeDLCoutput(outerGage_ExtremeValGage, 'ExtremeValues_OuterGage.xlsx');


%% ************************************************************************
% I. Finds the maximum Fx, ..., Mx ... values, and 6-component set of (Fi,Mi)
% *************************************************************************

% create the span location matrices for maximum/minimum forces and moments
namesGage = fieldnames(ExtremeValGage);

maxTbl = {};
for ii = 1:length(rootNamesMax)
    maxTbl = [maxTbl; ExtremeValGage.(rootNamesMax{ii})];
    columnVarMax{ii,1} = rootNamesMax{ii};
end
maxTbl.variable = columnVarMax

minTbl = {};
for ii = 1:length(rootNamesMin)
    minTbl = [minTbl; ExtremeValGage.(rootNamesMin{ii})];
    columnVarMin{ii,1} = rootNamesMin{ii};
end
minTbl.variable = columnVarMin

amplTbl = {};
for ii = 1:length(rootNamesAmpl)
    amplTbl = [amplTbl; ExtremeValGage.(rootNamesAmpl{ii})];
    columnVarAmpl{ii,1} = rootNamesAmpl{ii};
end
amplTbl.variable = columnVarAmpl

% create a table that summarizes the maximum/minimum values at each
% span location
outputMaxMin.('Root_Max') = maxTbl;
outputMaxMin.('Root_Min') = minTbl;
outputMaxMin.('Root_Ampl') = amplTbl;

% Blade gage spanwise location sensors
for spn = 1:innerGageStations+outerGageStations
    spnNamesMax = {['MaxSpn' num2str(spn) 'FLxb1']; ['MaxSpn' num2str(spn) 'FLyb1']; ...
        ['MaxSpn' num2str(spn) 'FLzb1']; ['MaxRes_Spn' num2str(spn) 'MLrb1_0deg'];...
        ['MaxRes_Spn' num2str(spn) 'MLrb1_90deg']; ['MaxSpn' num2str(spn) 'MLzb1']}
    
    spnNamesMin = {['MinSpn' num2str(spn) 'FLxb1']; ['MinSpn' num2str(spn) 'FLyb1']; ...
        ['MinSpn' num2str(spn) 'FLzb1']; ['MinRes_Spn' num2str(spn) 'MLrb1_0deg'];...
        ['MinRes_Spn' num2str(spn) 'MLrb1_90deg']; ['MinSpn' num2str(spn) 'MLzb1']}
    spnNamesAmpl = {['AmplSpn' num2str(spn) 'FLxb1FLyb1']; ['AmplSpn' num2str(spn) 'MLxb1MLyb1']};
    
    maxTbl = {};
    for ii = 1:length(spnNamesMax)
        maxTbl = [maxTbl; ExtremeValGage.(spnNamesMax{ii})];
        columnVarMax{ii,1} = spnNamesMax{ii};
    end
    maxTbl.variable = columnVarMax
    
    minTbl = {};
    for ii = 1:length(spnNamesMin)
        minTbl = [minTbl; ExtremeValGage.(spnNamesMin{ii})];
        columnVarMin{ii,1} = spnNamesMin{ii};
    end
    minTbl.variable = columnVarMin
    
    amplTbl = {};
    for ii = 1:length(spnNamesAmpl)
        amplTbl = [amplTbl; ExtremeValGage.(spnNamesAmpl{ii})];
        columnVarAmpl{ii,1} = spnNamesAmpl{ii};
    end
    amplTbl.variable = columnVarAmpl
    
    % create a table that summarizes the maximum/minimum values at each
    % span location
    outputMaxMin.(['Spn' num2str(spn) '_Max']) = maxTbl;
    outputMaxMin.(['Spn' num2str(spn) '_Min']) = minTbl;
    outputMaxMin.(['Spn' num2str(spn) '_Ampl']) = amplTbl;    
end

writeDLCoutput(outputMaxMin, 'SpanStation_MaxMin.xlsx', '1');




%% create the span location matrices for maximum/minimum forces and moments
% with the 6-component set
clear finalSave

% create the root 6x6 component set
gagefolder = 'out-gageSet1inner';
spnAct = spn;
% -------------- Max Root Forces/Moments -------------- %
clear FMxyz
MaxOrMin = 'Max'
for ii = 1:length(outputMaxMin.(['Root_' MaxOrMin]).data)
    fname = outputMaxMin.(['Root_' MaxOrMin]).file{ii};
    ftime = outputMaxMin.(['Root_' MaxOrMin]).time(ii);
    compLabels = {'RootFxb1'; 'RootFyb1'; 'RootFzb1'; 'RootMxb1'; 'RootMyb1'; 'RootMzb1'}
    out = loadFASTOutData([gagefolder fname(4:end)]);
    % set up the 6x6 matrix of forces and moments - diagonal is the
    % maximum for each force/moment.
    for jj = 1:length(compLabels)
        i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
        FMxyz(ii,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
    end
    finalSave.(['Root_' MaxOrMin]) = array2table(FMxyz,'VariableNames', compLabels)
end
finalSave.(['Root_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fx'] [MaxOrMin ' Fy']...
    [MaxOrMin ' Fz'] [MaxOrMin ' Mx'] [MaxOrMin ' My'] [MaxOrMin ' Mz']}
% -------------- Min Root Forces/Moments -------------- %
clear FMxyz
MaxOrMin = 'Min'
for ii = 1:length(outputMaxMin.(['Root_' MaxOrMin]).data)
    fname = outputMaxMin.(['Root_' MaxOrMin]).file{ii};
    ftime = outputMaxMin.(['Root_' MaxOrMin]).time(ii);
    compLabels = {'RootFxb1'; 'RootFyb1'; 'RootFzb1'; 'RootMxb1'; 'RootMyb1'; 'RootMzb1'}
    out = loadFASTOutData([gagefolder fname(4:end)]);
    % set up the 6x6 matrix of forces and moments - diagonal is the
    % maximum for each force/moment.
    for jj = 1:length(compLabels)
        i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
        FMxyz(ii,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
    end
    finalSave.(['Root_' MaxOrMin]) = array2table(FMxyz,'VariableNames', compLabels)
end
finalSave.(['Root_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fx'] [MaxOrMin ' Fy']...
    [MaxOrMin ' Fz'] [MaxOrMin ' Mx'] [MaxOrMin ' My'] [MaxOrMin ' Mz']}
% -------------- Ampl Root Forces/Moments -------------- %
clear FMxyz
MaxOrMin = 'Ampl'
for spn = 1:innerGageStations+outerGageStations
    finalSave.(['Root_' MaxOrMin]) = outputMaxMin.(['Root_' MaxOrMin])
    finalSave.(['Root_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fxy'] [MaxOrMin ' Mxy']}
end
% --------------- Blade Gage Stations ---------------- %
% calculate the matrix for the maximum values.
MaxOrMin = 'Max'
for spn = 1:innerGageStations+outerGageStations 
    if spn <= innerGageStations
        gagefolder = 'out-gageSet1inner';
        spnAct = spn;
    else
        gagefolder = 'out-gageSet2outer';
        spnAct = spn - innerGageStations;
    end    
    clear FMxyz
    for ii = 1:length(outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).data)                  
        fname = outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).file{ii};
        ftime = outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).time(ii);        
        compLabels = {['Spn' num2str(spnAct) 'FLxb1']; ['Spn' num2str(spnAct) 'FLyb1'];...
            ['Spn' num2str(spnAct) 'FLzb1']; ['Spn' num2str(spnAct) 'MLxb1']; ...
            ['Spn' num2str(spnAct) 'MLyb1']; ['Spn' num2str(spnAct) 'MLzb1']};                  
        out = loadFASTOutData([gagefolder fname(4:end)]);        
        % set up the 6x6 matrix of forces and moments - diagonal is the
        % maximum for each force/moment.
        for jj = 1:length(compLabels)
            i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
            FMxyz(ii,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));            
        end                
        finalSave.(['Spn' num2str(spn) '_' MaxOrMin]) = array2table(FMxyz,'VariableNames', compLabels)         
    end
    finalSave.(['Spn' num2str(spn) '_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fx'] [MaxOrMin ' Fy']...
        [MaxOrMin ' Fz'] [MaxOrMin ' Mx'] [MaxOrMin ' My'] [MaxOrMin ' Mz']}
end

% calculate the matrix for the minimum values.
MaxOrMin = 'Min'
for spn = 1:innerGageStations+outerGageStations 
    if spn <= innerGageStations
        gagefolder = 'out-gageSet1inner';
        spnAct = spn;
    else
        gagefolder = 'out-gageSet2outer';
        spnAct = spn - innerGageStations;
    end  
    clear FMxyz
    for ii = 1:length(outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).data)                  
        fname = outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).file{ii};
        ftime = outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin]).time(ii);        
        compLabels = {['Spn' num2str(spnAct) 'FLxb1']; ['Spn' num2str(spnAct) 'FLyb1'];...
            ['Spn' num2str(spnAct) 'FLzb1']; ['Spn' num2str(spnAct) 'MLxb1']; ...
            ['Spn' num2str(spnAct) 'MLyb1']; ['Spn' num2str(spnAct) 'MLzb1']};                  
        out = loadFASTOutData([gagefolder fname(4:end)]);        
        % set up the 6x6 matrix of forces and moments - diagonal is the
        % maximum for each force/moment.
        for jj = 1:length(compLabels)
            i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
            FMxyz(ii,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));            
        end                
        finalSave.(['Spn' num2str(spn) '_' MaxOrMin]) = array2table(FMxyz,'VariableNames', compLabels)         
    end
    finalSave.(['Spn' num2str(spn) '_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fx'] [MaxOrMin ' Fy']...
        [MaxOrMin ' Fz'] [MaxOrMin ' Mx'] [MaxOrMin ' My'] [MaxOrMin ' Mz']}
end

% calculate the matrix for the minimum values.
MaxOrMin = 'Ampl'
for spn = 1:innerGageStations+outerGageStations
    finalSave.(['Spn' num2str(spn) '_' MaxOrMin]) = outputMaxMin.(['Spn' num2str(spn) '_' MaxOrMin])
    finalSave.(['Spn' num2str(spn) '_' MaxOrMin]).Properties.RowNames = {[MaxOrMin ' Fxy'] [MaxOrMin ' Mxy']}
end


%% Save the final output

filename = '1_ComponentSets_MaxMin.xlsx';

% add the root sheet
writetable(finalSave.('Root_Max'), filename, 'Sheet', ...
    'Root', 'Range', 'A1:G7','WriteRowNames',1)
writetable(finalSave.('Root_Min'), filename, 'Sheet', ...
    'Root', 'Range', 'A9:G15','WriteRowNames',1)
writetable(finalSave.('Root_Ampl'), filename, 'Sheet', ...
    'Root', 'Range', 'A17:G19','WriteRowNames',1)

% add a sheet for each spanwise blade gage station
for spn = 1:innerGageStations+outerGageStations
    writetable(finalSave.(['Spn' num2str(spn) '_Max']), filename, 'Sheet', ...
        ['Spn' num2str(spn)], 'Range', 'A1:G7','WriteRowNames',1)
    writetable(finalSave.(['Spn' num2str(spn) '_Min']), filename, 'Sheet', ...
        ['Spn' num2str(spn)], 'Range', 'A9:G15','WriteRowNames',1)
    writetable(finalSave.(['Spn' num2str(spn) '_Ampl']), filename, 'Sheet', ...
        ['Spn' num2str(spn)], 'Range', 'A17:G19','WriteRowNames',1)
end



%% ************************************************************************
% II. This section saves the resultant moment maximum and minimum data
% *************************************************************************
theta = 0:15:165;
clear ResultantMomentSave

% ---------------- Root Resultant Moments ------------------ %
% Save the maximum resultant moments
MaxOrMin = 'Max';
tbl = {}; columnVar = {}
for ii = theta
    tbl = [tbl; ExtremeValGage.([MaxOrMin 'Res_RootMrb1_' num2str(ii) 'deg'])];
    columnVar{end+1,1} = [MaxOrMin 'Res_RootMrb1_' num2str(ii) 'deg'];
end
ResultantMomentSave.([MaxOrMin 'Res_Root']) = tbl
ResultantMomentSave.([MaxOrMin 'Res_Root']).Properties.RowNames = columnVar;   
    ResultantMomentSave.([MaxOrMin 'Res_Root']).Properties.VariableNames{2} = ...
        [MaxOrMin 'Res_RootMrb1'];

% Save the minimum resultant moments
MaxOrMin = 'Min';
tbl = {}; columnVar = {}
for ii = theta
    tbl = [tbl; ExtremeValGage.([MaxOrMin 'Res_RootMrb1_' num2str(ii) 'deg'])];
    columnVar{end+1,1} = [MaxOrMin 'Res_RootMrb1_' num2str(ii) 'deg'];
end
ResultantMomentSave.([MaxOrMin 'Res_Root']) = tbl
ResultantMomentSave.([MaxOrMin 'Res_Root']).Properties.RowNames = columnVar;
ResultantMomentSave.([MaxOrMin 'Res_Root']).Properties.VariableNames{2} = ...
    [MaxOrMin 'Res_RootMrb1'];

% ---------------- Gage Resultant Moments ------------------ %
for spn = 1:innerGageStations+outerGageStations        
    % Save the maximum resultant moments
    MaxOrMin = 'Max';
    tbl = {}; columnVar = {}
    for ii = theta        
        tbl = [tbl; ExtremeValGage.([MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1_' num2str(ii) 'deg'])];
        columnVar{end+1,1} = [MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1_' num2str(ii) 'deg'];
    end
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]) = tbl
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]).Properties.RowNames = columnVar;   
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]).Properties.VariableNames{2} = ...
        [MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1'];
    
    % Save the minimum resultant moments
    MaxOrMin = 'Min';
    tbl = {}; columnVar = {}
    for ii = theta        
        tbl = [tbl; ExtremeValGage.([MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1_' num2str(ii) 'deg'])];
        columnVar{end+1,1} = [MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1_' num2str(ii) 'deg'];
    end
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]) = tbl
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]).Properties.RowNames = columnVar;   
    ResultantMomentSave.([MaxOrMin 'Res_Spn' num2str(spn)]).Properties.VariableNames{2} = ...
        [MaxOrMin 'Res_Spn' num2str(spn) 'MLrb1'];
end

writeDLCoutput(ResultantMomentSave, '2_ResulantMomentStationsResults.xlsx','2')




%% ************************************************************************
% V. This section saves the forces during the moment of max blade tip
% deflection
% *************************************************************************

% ------------------ Maximum Tip Deflection -------------------- %
innerGage.MaxOoPDefl
innerGage_ExtremeValGage.MaxOoPDefl
innerGage_ExtremeValGage.MinOoPDefl

outerGage.MaxOoPDefl
outerGage_ExtremeValGage.MaxOoPDefl
outerGage_ExtremeValGage.MinOoPDefl

% calculate the matrix for the maximum out of plane deflection
clear FMxyz; columnVector = {};
fname = innerGage_ExtremeValGage.MaxOoPDefl.file{1};
ftime = innerGage_ExtremeValGage.MaxOoPDefl.time(1);

out = loadFASTOutData(['out-gageSet1inner' fname(4:end)]);
% ---------------- root loads ------------------ %
compLabels = {'RootFxb1'; 'RootFyb1'; 'RootFzb1'; 'RootMxb1'; 'RootMyb1'; 'RootMzb1'};
columnVector{end+1} = 'Root';
for jj = 1:length(compLabels)
    i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
    FMxyz(1,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
end

for spn = 1:innerGageStations+outerGageStations
    if spn <= innerGageStations
        fname = innerGage_ExtremeValGage.MaxOoPDefl.file{1};
        ftime = innerGage_ExtremeValGage.MaxOoPDefl.time(1);
        gageFolder = 'out-gageSet1inner';
        spnAct = spn;
    else
        fname = outerGage_ExtremeValGage.MaxOoPDefl.file{1};
        ftime = outerGage_ExtremeValGage.MaxOoPDefl.time(1);
        gageFolder = 'out-gageSet2outer';
        spnAct = spn - innerGageStations;
    end    
    out = loadFASTOutData([gageFolder fname(4:end)]);
    
    compLabels = {['Spn' num2str(spnAct) 'FLxb1']; ['Spn' num2str(spnAct) 'FLyb1'];...
        ['Spn' num2str(spnAct) 'FLzb1']; ['Spn' num2str(spnAct) 'MLxb1']; ...
        ['Spn' num2str(spnAct) 'MLyb1']; ['Spn' num2str(spnAct) 'MLzb1']};
    columnVector{end+1} = ['Spn_' num2str(spn)];
    
    % set up the 6x6 matrix of forces and moments - diagonal is the
    % maximum for each force/moment.
    for jj = 1:length(compLabels)
        i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
        FMxyz(spn+1,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
    end
end
rowNames = {'Fx' 'Fy' 'Fz' 'Mx' 'My' 'Mz'};
DeflLoads_Final.FMxyz_maxDefl = array2table(FMxyz,'VariableNames', rowNames)
DeflLoads_Final.FMxyz_maxDefl.Properties.RowNames = columnVector


% ------------------ Minimum Tip Clearance -------------------- %
innerGage.MinTipClrnc
innerGage_ExtremeValGage.MinTipClrnc

outerGage.MinTipClrnc
outerGage_ExtremeValGage.MinTipClrnc

% calculate the matrix for the minimum tip clearance
clear FMxyz; columnVector = {};
fname = innerGage_ExtremeValGage.MinTipClrnc.file{1};
ftime = innerGage_ExtremeValGage.MinTipClrnc.time(1);

out = loadFASTOutData(['out-gageSet1inner' fname(4:end)]);
% ---------------- root loads ------------------ %
compLabels = {'RootFxb1'; 'RootFyb1'; 'RootFzb1'; 'RootMxb1'; 'RootMyb1'; 'RootMzb1'};
columnVector{end+1} = 'Root';
for jj = 1:length(compLabels)
    i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
    FMxyz(1,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
end

for spn = 1:innerGageStations+outerGageStations
    if spn <= innerGageStations
        fname = innerGage_ExtremeValGage.MinTipClrnc.file{1};
        ftime = innerGage_ExtremeValGage.MinTipClrnc.time(1);
        gageFolder = 'out-gageSet1inner';
        spnAct = spn;
    else
        fname = outerGage_ExtremeValGage.MinTipClrnc.file{1};
        ftime = outerGage_ExtremeValGage.MinTipClrnc.time(1);
        gageFolder = 'out-gageSet2outer';
        spnAct = spn - innerGageStations;
    end    
    out = loadFASTOutData([gageFolder fname(4:end)]);
    
    compLabels = {['Spn' num2str(spnAct) 'FLxb1']; ['Spn' num2str(spnAct) 'FLyb1'];...
        ['Spn' num2str(spnAct) 'FLzb1']; ['Spn' num2str(spnAct) 'MLxb1']; ...
        ['Spn' num2str(spnAct) 'MLyb1']; ['Spn' num2str(spnAct) 'MLzb1']};
    columnVector{end+1} = ['Spn_' num2str(spn)];
    
    % set up the 6x6 matrix of forces and moments - diagonal is the
    % maximum for each force/moment.
    for jj = 1:length(compLabels)
        i_ftime = find(out.data(:, strcmp('Time',out.list)) >= ftime,1);
        FMxyz(spn+1,jj) = out.data(i_ftime, strcmp(compLabels(jj),out.list));
    end
end
rowNames = {'Fx' 'Fy' 'Fz' 'Mx' 'My' 'Mz'};
DeflLoads_Final.FMxyz_minTipClrnc = array2table(FMxyz,'VariableNames', rowNames)
DeflLoads_Final.FMxyz_minTipClrnc.Properties.RowNames = columnVector


writeDLCoutput(DeflLoads_Final, '5_MaximumDeflectionLoadsParked.xlsx','3')







