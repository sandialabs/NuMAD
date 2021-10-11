function output=findExtremes_1p1peaks(outNames,str,thetaMomentRotation,gageCoordinateRotation,caseName,output,plotFlag)
outputId=length(getfield(output,str(1).fieldName))+1;

for i=1:length(str)
    switch lower(str(i).minOrMax)  % initialize the data in the structure
        case 'max'
            output=setfield(output,str(i).fieldName,{outputId},'data',nan);
        case 'min'
            output=setfield(output,str(i).fieldName,{outputId},'data',nan);
        case 'amplitude'
            output=setfield(output,str(i).fieldName,{outputId},'data',nan);
        case 'resultant max'
            for theta = thetaMomentRotation
                thStr = ['_' num2str(theta) 'deg'];
                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',nan);
            end
        case 'resultant min'
            for theta = thetaMomentRotation
                thStr = ['_' num2str(theta) 'deg'];
                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',nan);
            end
        case 'combined strain max'
            output=setfield(output,str(i).fieldName,{outputId},'data',nan);
        case 'combined strain min'
            output=setfield(output,str(i).fieldName,{outputId},'data',nan);
        otherwise
            error('error in findExtremes.m')
    end
end

for n=1:size(outNames,1)
    out=loadFASTOutData(outNames{n});  % load a file from the list
    
    % correct for spanwise gage channels which are in the local blade
    % coordinate system oriented along the principal axis
    for gg = 1:length(gageCoordinateRotation)
        % determine the rotation matrix
        theta = gageCoordinateRotation(gg);
        rotMat = [cosd(theta) sind(theta); -sind(theta) cosd(theta)];
        
        % force transformation
        fxData = out.data(:,strcmp(['Spn' num2str(gg) 'FLxb1'],out.list));
        fyData = out.data(:,strcmp(['Spn' num2str(gg) 'FLyb1'],out.list));     
        Ftransform = rotMat*[fxData'; fyData'];
        out.data(:,strcmp(['Spn' num2str(gg) 'FLxb1'],out.list)) = Ftransform(1,:)';
        out.data(:,strcmp(['Spn' num2str(gg) 'FLyb1'],out.list)) = Ftransform(2,:)';
        out.tbl.(['Spn' num2str(gg) 'FLxb1']) = Ftransform(1,:)';
        out.tbl.(['Spn' num2str(gg) 'FLyb1']) = Ftransform(2,:)';
        
        % moment transformation
        mxData = out.data(:,strcmp(['Spn' num2str(gg) 'MLxb1'],out.list));
        myData = out.data(:,strcmp(['Spn' num2str(gg) 'MLyb1'],out.list));
        Mtransform = rotMat*[mxData'; myData'];
        out.data(:,strcmp(['Spn' num2str(gg) 'MLxb1'],out.list)) = Mtransform(1,:)';
        out.data(:,strcmp(['Spn' num2str(gg) 'MLyb1'],out.list)) = Mtransform(2,:)';
        out.tbl.(['Spn' num2str(gg) 'MLxb1']) = Mtransform(1,:)';
        out.tbl.(['Spn' num2str(gg) 'MLyb1']) = Mtransform(2,:)';
    end
    
    for i=1:length(str) % ble: number of entries in sortExtremes
        
        for ch=1:length(str(i).chanLabels)
            
            pointer=strcmp(str(i).chanLabels{ch},out.list); % pointer to proper column in FAST data output
            nums=out.data(:,pointer)*str(i).gains(ch); % ble: output from str{i}
            time=out.tbl.Time;
            
            switch lower(str(i).minOrMax)
                case 'max'
                    % identify the maximum peaks based on the IEC 1.1 standard
                    isMinimum = 0;
                    tmp = findPeaks1p1(nums,str(i).fieldName,isMinimum,time);                    
                    % keep track of each variable being used and associated peak data
                    if ch == 1, output.(str(i).fieldName).Chan = str(i).chanLabels{ch};
                    else, output.(str(i).fieldName).Chan = [output.(str(i).fieldName).Chan '; ' str(i).chanLabels{ch}];
                    end
                    % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                    output=setfield(output,str(i).fieldName,{outputId},'File','all');
                    output=setfield(output,str(i).fieldName,{outputId},'time','fifty year');
                    output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);                                        
                    % append the identified peaks to the data variable
                    if n == 1, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = tmp;
                    else, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = ...
                            [output.(str(i).fieldName).(['peaks_' num2str(ch)]) tmp]; 
                    end
                    
                case 'min'
                    % identify the minimum peaks based on the IEC 1.1 standard
                    isMinimum = 1;
                    tmp = findPeaks1p1(nums,str(i).fieldName,isMinimum,time);                    
                    % keep track of each variable being used and associated peak data
                    if ch == 1, output.(str(i).fieldName).Chan = str(i).chanLabels{ch};
                    else, output.(str(i).fieldName).Chan = [output.(str(i).fieldName).Chan '; ' str(i).chanLabels{ch}];
                    end
                    % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                    output=setfield(output,str(i).fieldName,{outputId},'File','all');
                    output=setfield(output,str(i).fieldName,{outputId},'time','fifty year');
                    output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);                    
                    % append the identified peaks to the data variable
                    if n == 1, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = tmp;
                    else, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = ...
                            [output.(str(i).fieldName).(['peaks_' num2str(ch)]) tmp]; 
                    end
                    
                case 'amplitude' 
                    if ch == 1 % only perform this sum of squares amplitude analysis once per 'str' entry
                        chString = '';
                        numsVector = zeros(length(out.data),length(str(i).chanLabels));
                         
                        for chVector=1:length(str(i).chanLabels)
                            pointer=strcmp(str(i).chanLabels{chVector},out.list);
                            numsVector(:,chVector)=(out.data(:,pointer)*str(i).gains(chVector)).^2;
                            chString = [chString ' ' str(i).chanLabels{chVector}];
                        end
                        numsAmplitude = sqrt(sum(numsVector,2)); 
                        
                        % identify the maximum amplitude peaks based on the IEC 1.1 standard
                        isMinimum = 0;
                        tmp = findPeaks1p1(numsAmplitude,str(i).fieldName,isMinimum,time);                        
                        output=setfield(output,str(i).fieldName,{outputId},'Chan', ['Amplitude of:' chString]);
                        % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                        output=setfield(output,str(i).fieldName,{outputId},'File','all');
                        output=setfield(output,str(i).fieldName,{outputId},'time','fifty year');
                        output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);
                        % append the identified peaks to the data variable
                        if n == 1, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = tmp;
                        else, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = ...
                                [output.(str(i).fieldName).(['peaks_' num2str(ch)]) tmp];
                        end
                    end
                    
                case 'resultant max'
                    if ch == 1 % only perform this sum of squares amplitude analysis once per 'str' entry                        
                        for theta = thetaMomentRotation  % Rotation angle [deg]
                            thStr = ['_' num2str(theta) 'deg'];
                            % calculate the resultant moment at incremental rotation angles, theta.
                            MxCosTh = out.data(:,strcmp(str(i).chanLabels(1),out.list)) .* cosd(theta);
                            MySinTh = out.data(:,strcmp(str(i).chanLabels(2),out.list)) .* sind(theta);
                            Mresultant = MxCosTh + MySinTh;
                            
                            % identify the maximum amplitude peaks based on the IEC 1.1 standard
                            isMinimum = 0;
                            tmp = findPeaks1p1(Mresultant,str(i).fieldName,isMinimum,time);
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'Chan', ['Resultant Moment: ' str(i).chanLabels{1} ' ' str(i).chanLabels{2}]);
                            % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'File','all');
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'time','fifty year');
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'Name',caseName);
                            % append the identified peaks to the data variable
                            if n == 1, output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) = tmp;
                            else, output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) = ...
                                    [output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) tmp];
                            end
                        end
                    end                                 
                    
                case 'resultant min'
                    if ch == 1 % only perform this sum of squares amplitude analysis once per 'str' entry                        
                        for theta = thetaMomentRotation  % Rotation angle [deg]
                            thStr = ['_' num2str(theta) 'deg'];
                            % calculate the resultant moment at incremental rotation angles, theta.
                            MxCosTh = out.data(:,strcmp(str(i).chanLabels(1),out.list)) .* cosd(theta);
                            MySinTh = out.data(:,strcmp(str(i).chanLabels(2),out.list)) .* sind(theta);
                            Mresultant = MxCosTh + MySinTh;
                            
                            % identify the minimum amplitude peaks based on the IEC 1.1 standard
                            isMinimum = 1;
                            tmp = findPeaks1p1(Mresultant,str(i).fieldName,isMinimum,time);
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'Chan', ['Resultant Moment: ' str(i).chanLabels{1} ' ' str(i).chanLabels{2}]);
                            % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'File','all');
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'time','fifty year');
                            output=setfield(output,[str(i).fieldName thStr],{outputId},'Name',caseName);
                            % append the identified peaks to the data variable
                            if n == 1, output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) = tmp;
                            else, output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) = ...
                                    [output.([str(i).fieldName thStr]).(['peaks_' num2str(ch)]) tmp];
                            end
                        end
                    end 
                    
                case 'combined strain max'
                    % used for calculating the combined strain from bending stress and centrifugal stress
                    Ngages = length(str(i).chanLabels)/2;
                    if ch <= Ngages % only perform the calculation on the total number of gages
                        % calculate strain due to bending moment (eps_b) and due to axial loading (eps_n)
                        eps_b = out.data(:,strcmp(str(i).chanLabels(ch),out.list)) .* str(i).gains(ch);
                        eps_n = out.data(:,strcmp(str(i).chanLabels(ch+Ngages),out.list)) .* str(i).gains(ch+Ngages);
                        % check that moment and force channels are correct
                        if strcmp(str(i).chanLabels{ch}(1:4),str(i).chanLabels{ch+Ngages}(1:4)) ~= 1
                            error('not using the same spanwise station')
                        elseif ~contains(str(i).chanLabels{ch+Ngages}, 'F') || ~contains(str(i).chanLabels{ch+Ngages}, 'z') || ~contains(str(i).chanLabels{ch}, 'M')
                            % ensure that the second channel is Fz for spanwise force and M for first channel
                            error('should be Fz for centrifugal force and moment for bending (flap or edge)')
                        elseif (contains(str(i).fieldName,'HP') || contains(str(i).fieldName,'LP')) && ~contains(str(i).chanLabels{ch}, 'y')
                            error('should be My for flapwise bending moment')
                        elseif (contains(str(i).fieldName,'TE') || contains(str(i).fieldName,'LE')) && ~contains(str(i).chanLabels{ch}, 'x')
                            error('should be Mx for edgewise bending moment')
                        end
                        % combine the bending and normal strain contributions based on coordinate system
                        if contains(str(i).fieldName,'HP')
                            % calculate combined strain on high pressure side of spar cap                            
                            eps_comb = eps_n + eps_b;
                        elseif contains(str(i).fieldName,'LP')
                            % calculate combined strain on low pressure side of spar cap
                            eps_comb = eps_n - eps_b;
                        elseif contains(str(i).fieldName,'TE')
                            % calculate combined strain on trailing edge reinforcement                           
                            eps_comb = eps_n + eps_b;
                        elseif contains(str(i).fieldName,'LE')
                            % calculate combined strain on leading edge reinforcement 
                            eps_comb = eps_n - eps_b;
                        else % cannot identify location on blade
                            error('field name must identify location on blade')
                        end
                        
                        % identify the maximum amplitude peaks based on the IEC 1.1 standard
                        isMinimum = 0;
                        tmp = findPeaks1p1(eps_comb,str(i).fieldName,isMinimum,time);
                        output=setfield(output,str(i).fieldName,{outputId},'Chan', ...
                            ['Combined Load: ' str(i).chanLabels{ch} ' ' str(i).chanLabels{ch+Ngages}]);
                        % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                        output=setfield(output,str(i).fieldName,{outputId},'File','all');
                        output=setfield(output,str(i).fieldName,{outputId},'time','fifty year');
                        output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);    
                        % append the identified peaks to the data variable
                        if n == 1, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = tmp;
                        else, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = ...
                                [output.(str(i).fieldName).(['peaks_' num2str(ch)]) tmp];
                        end                    
                    end
                    
                case 'combined strain min'
                    % used for calculating the combined strain from bending stress and centrifugal stress
                    Ngages = length(str(i).chanLabels)/2;
                    if ch <= Ngages % only perform the calculation on the total number of gages
                        % calculate strain due to bending moment (eps_b) and due to axial loading (eps_n)
                        eps_b = out.data(:,strcmp(str(i).chanLabels(ch),out.list)) .* str(i).gains(ch);
                        eps_n = out.data(:,strcmp(str(i).chanLabels(ch+Ngages),out.list)) .* str(i).gains(ch+Ngages);
                        % check that moment and force channels are correct
                        if strcmp(str(i).chanLabels{ch}(1:4),str(i).chanLabels{ch+Ngages}(1:4)) ~= 1
                            error('not using the same spanwise station')
                        elseif ~contains(str(i).chanLabels{ch+Ngages}, 'F') || ~contains(str(i).chanLabels{ch+Ngages}, 'z') || ~contains(str(i).chanLabels{ch}, 'M')
                            % ensure that the second channel is Fz for spanwise force and M for first channel
                            error('should be Fz for centrifugal force and moment for bending (flap or edge)')
                        elseif (contains(str(i).fieldName,'HP') || contains(str(i).fieldName,'LP')) && ~contains(str(i).chanLabels{ch}, 'y')
                            error('should be My for flapwise bending moment')
                        elseif (contains(str(i).fieldName,'TE') || contains(str(i).fieldName,'LE')) && ~contains(str(i).chanLabels{ch}, 'x')
                            error('should be Mx for edgewise bending moment')
                        end
                        % combine the bending and normal strain contributions based on coordinate system
                        if contains(str(i).fieldName,'HP')
                            % calculate combined strain on high pressure side of spar cap                            
                            eps_comb = eps_n + eps_b;
                        elseif contains(str(i).fieldName,'LP')
                            % calculate combined strain on low pressure side of spar cap
                            eps_comb = eps_n - eps_b;
                        elseif contains(str(i).fieldName,'TE')
                            % calculate combined strain on trailing edge reinforcement                           
                            eps_comb = eps_n + eps_b;
                        elseif contains(str(i).fieldName,'LE')
                            % calculate combined strain on leading edge reinforcement 
                            eps_comb = eps_n - eps_b;
                        else % cannot identify location on blade
                            error('field name must identify location on blade')
                        end
                        
                        % identify the minimum amplitude peaks based on the IEC 1.1 standard
                        isMinimum = 1;
                        tmp = findPeaks1p1(eps_comb,str(i).fieldName,isMinimum,time);
                        output=setfield(output,str(i).fieldName,{outputId},'Chan', ...
                            ['Combined Load: ' str(i).chanLabels{ch} ' ' str(i).chanLabels{ch+Ngages}]);
                        % IEC DLC 1.1 is a 50 year extrapolation of the time series simulation set
                        output=setfield(output,str(i).fieldName,{outputId},'File','all');
                        output=setfield(output,str(i).fieldName,{outputId},'time','fifty year');
                        output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);
                        % append the identified peaks to the data variable
                        if n == 1, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = tmp;
                        else, output.(str(i).fieldName).(['peaks_' num2str(ch)]) = ...
                                [output.(str(i).fieldName).(['peaks_' num2str(ch)]) tmp];
                        end
                    end
                        
                otherwise
                    error('set to min or max in sortExtremes.m')
            end
        end
    end
end

end


function [binMaxPeak] = findPeaks1p1(data,labels,isMinimum,time)

plotTimeSeries = 0;
dataMean = mean(data);
dataStdev = std(data);

% find the maximum/minimum value peaks for load extrapolation
if isMinimum == 1 % find maximum peaks (by inverting the data channel)
    gains = -1;
else % find minimum peaks
    gains = 1;
end
dataMax = gains.*data;
dataMax(dataMax<gains*dataMean+1.4*dataStdev) = nan;

index1 = find(~isnan(dataMax)==1);
index2 = find(diff(index1)>1);
index3 = index1(index2);
if ~isempty(index3)
    indexRangeMax = [[index1(1);index1(index2+1)] [index3;index1(end)]];
    
    for kk = 1:size(indexRangeMax,1)
        [binMaxPeak(kk), peakIndex(kk,1)] = max(dataMax(indexRangeMax(kk,1):indexRangeMax(kk,2)));
        binMaxPeak(kk) = gains.*binMaxPeak(kk);
    end
else
    binMaxPeak(1) = max(gains.*data);
    binMaxPeak(1) = gains.*binMaxPeak(1);
end

if plotTimeSeries
    figure(1001); set(gcf,'Position',[50 50 1000 800]); clf; hold on; grid on
    plot(time,data,'.-')
    try % plot peaks if they exist
        plot(time(index1),data(index1),'x')
        plot(time(indexRangeMax(:,1)+peakIndex-1),binMaxPeak,'o')
    catch % otherwise continue
    end
    title(labels)
end

end
