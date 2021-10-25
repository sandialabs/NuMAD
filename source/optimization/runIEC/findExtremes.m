function output=findExtremes(outNames,str,thetaMomentRotation,gageCoordinateRotation,caseName,output,plotFlag)

outputId=length(getfield(output,str(1).fieldName))+1;

for i=1:length(str)
    switch lower(str(i).minOrMax)  % initialize the data in the structure
        case 'max'
            output=setfield(output,str(i).fieldName,{outputId},'data',-Inf);
        case 'min'
            output=setfield(output,str(i).fieldName,{outputId},'data',Inf);
        case 'amplitude'
            output=setfield(output,str(i).fieldName,{outputId},'data',-Inf);
        case 'resultant max'
            for theta = thetaMomentRotation
                thStr = ['_' num2str(theta) 'deg'];
                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',-Inf);
            end
        case 'resultant min'
            for theta = thetaMomentRotation
                thStr = ['_' num2str(theta) 'deg'];
                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',Inf);
            end
        case 'combined strain max'
            output=setfield(output,str(i).fieldName,{outputId},'data',-Inf);
        case 'combined strain min'
            output=setfield(output,str(i).fieldName,{outputId},'data',Inf);
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

            switch lower(str(i).minOrMax)
                case 'max'
                    tmp=max(nums);
                    if tmp>getfield(output,str(i).fieldName,{length(getfield(output,str(i).fieldName))},'data')  % if current maximum is less than overall maximum
                        output=setfield(output,str(i).fieldName,{outputId},'data',tmp);
                        output=setfield(output,str(i).fieldName,{outputId},'File',outNames{n});
                        output=setfield(output,str(i).fieldName,{outputId},'Chan',str(i).chanLabels{ch});
                        timepointer=find(nums==tmp); % find pointer to time of max data point
                        output=setfield(output,str(i).fieldName,{outputId},'time',out.data(timepointer(1),1));
                        output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);                        
                    end
                    
                case 'min'
                    tmp=min(nums);
                    if tmp<getfield(output,str(i).fieldName,{length(getfield(output,str(i).fieldName))},'data')  % if current maximum is less than overall maximum
                        output=setfield(output,str(i).fieldName,{outputId},'data',tmp);
                        output=setfield(output,str(i).fieldName,{outputId},'File',outNames{n});
                        output=setfield(output,str(i).fieldName,{outputId},'Chan',str(i).chanLabels{ch});
                        timepointer=find(nums==tmp); % find pointer to time of max data point
                        output=setfield(output,str(i).fieldName,{outputId},'time',out.data(timepointer(1),1));
                        output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);                        
                    end
                    
                case 'amplitude'                    
                    if ch == 1 % only perform this sum of squares amplitude analysis once per 'str' entry
                        chString = ''; numsVector = []; numsAmplitude = []; 
                        numsVector = zeros(length(out.data),length(str(i).chanLabels));
                         
                        for chVector=1:length(str(i).chanLabels)
                            pointer=strcmp(str(i).chanLabels{chVector},out.list);
                            numsVector(:,chVector)=(out.data(:,pointer)*str(i).gains(chVector)).^2;
                            chString = [chString ' ' str(i).chanLabels{chVector}];
                        end
                        numsAmplitude = sqrt(sum(numsVector,2));
                        tmp = max(numsAmplitude);
                        if tmp>getfield(output,str(i).fieldName,{length(getfield(output,str(i).fieldName))},'data')  % if current maximum is less than overall maximum
                            output=setfield(output,str(i).fieldName,{outputId},'data',tmp);
                            output=setfield(output,str(i).fieldName,{outputId},'File',outNames{n});
                            output=setfield(output,str(i).fieldName,{outputId},'Chan', ['Amplitude of:' chString]);
                            timepointer=find(numsAmplitude==tmp); % find pointer to time of max data point
                            output=setfield(output,str(i).fieldName,{outputId},'time',out.data(timepointer(1),1));
                            output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);
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
                            tmp = max(Mresultant);
                            if tmp>getfield(output,[str(i).fieldName thStr],{length(getfield(output,[str(i).fieldName thStr]))},'data')  % if current maximum is greater than overall maximum
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',tmp);
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'File',outNames{n});
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'Chan', ['Resultant Moment: ' str(i).chanLabels{1} ' ' str(i).chanLabels{2}]);
                                timepointer=find(Mresultant==tmp); % find pointer to time of max data point
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'time',out.data(timepointer(1),1));
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'Name',caseName);
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
                            tmp = min(Mresultant);
                            if tmp<getfield(output,[str(i).fieldName thStr],{length(getfield(output,[str(i).fieldName thStr]))},'data')  % if current maximum is greater than overall maximum
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'data',tmp);
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'File',outNames{n});
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'Chan', ['Resultant Moment: ' str(i).chanLabels{1} ' ' str(i).chanLabels{2}]);
                                timepointer=find(Mresultant==tmp); % find pointer to time of max data point
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'time',out.data(timepointer(1),1));
                                output=setfield(output,[str(i).fieldName thStr],{outputId},'Name',caseName);
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
                        tmp=max(eps_comb);
                        if tmp>getfield(output,str(i).fieldName,{length(getfield(output,str(i).fieldName))},'data')  % if current maximum is less than overall maximum
                            output=setfield(output,str(i).fieldName,{outputId},'data',tmp);
                            output=setfield(output,str(i).fieldName,{outputId},'File',outNames{n});
                            output=setfield(output,str(i).fieldName,{outputId},'Chan',str(i).chanLabels{ch});
                            timepointer=find(eps_comb==tmp); % find pointer to time of max data point
                            output=setfield(output,str(i).fieldName,{outputId},'time',out.data(timepointer(1),1));
                            output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);
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
                        tmp=min(eps_comb);
                        if tmp<getfield(output,str(i).fieldName,{length(getfield(output,str(i).fieldName))},'data')  % if current maximum is less than overall maximum
                            output=setfield(output,str(i).fieldName,{outputId},'data',tmp);
                            output=setfield(output,str(i).fieldName,{outputId},'File',outNames{n});
                            output=setfield(output,str(i).fieldName,{outputId},'Chan',str(i).chanLabels{ch});
                            timepointer=find(eps_comb==tmp); % find pointer to time of max data point
                            output=setfield(output,str(i).fieldName,{outputId},'time',out.data(timepointer(1),1));
                            output=setfield(output,str(i).fieldName,{outputId},'Name',caseName);
                        end
                    end
                        
                otherwise
                    error('set to min or max in sortExtremes.m')
            end
        end
    end
end

if plotFlag % plot results
    for i=1:length(str)
        figure
        outfn=getfield(output,str(i).fieldName,{outputId},'File');
        out=loadFASTOutData(outfn);  % load a file from the list
        chLab=getfield(output,str(i).fieldName,{outputId},'Chan');
        pointer=strcmp(chLab,out.list); % pointer to proper column in FAST data output
        chId=strcmp(chLab,str(i).chanLabels);
        x1=out.data(:,1);
        y1=out.data(:,pointer);
        x2=getfield(output,str(i).fieldName,{outputId},'time') ;
        gain=str(i).gains(chId);
        y2=getfield(output,str(i).fieldName,{outputId},'data') / gain;
        plot(x1,y1,x2,y2,'o')
        xlabel(out.list{1})
        ylabel(out.list{pointer})
        fn=getfield(output,str(i).fieldName,{outputId},'File');
        title({caseName,strrep(fn,'_','\_')})
        fn=fn(5:end-4);
        print('-dpng',sprintf('figs/%s_%s.png',fn,out.list{pointer}));
    end
end

end
