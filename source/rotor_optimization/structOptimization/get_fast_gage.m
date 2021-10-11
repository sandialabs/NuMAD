function fast_gage=get_fast_gage(rotationAngle)
    % rotationAngle - Loads analysis resolution angle in degrees
    
    thetaMomentRotation = 0:rotationAngle:180; % [deg]
    thetaMomentRotation(thetaMomentRotation==180)=[]; % remove repetitive 180 deg rotation
    gage_labels={}; gage_maximum=[]; gage_theta=[];
    for bb = 1 % blade number(s)
        for tt = 1:length(thetaMomentRotation)
            for gg = 0:9 % number of gages
                % save the gage variable names
                if gg==0
                    % create the maximum moment vector
                    gage_labels{bb}{2*(tt-1)+1}{gg+1} = ['MaxRes_RootMrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                    gage_maximum{bb}{2*(tt-1)+1}(gg+1) = 1;
                    gage_theta{bb}{2*(tt-1)+1}(gg+1) = thetaMomentRotation(tt);
                    % create the minimum moment vector
                    gage_labels{bb}{2*(tt-1)+2}{gg+1}  = ['MinRes_RootMrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                    gage_maximum{bb}{2*(tt-1)+2}(gg+1) = -1;
                    gage_theta{bb}{2*(tt-1)+2}(gg+1) = thetaMomentRotation(tt);
                else
                    % create the maximum moment vector
                    gage_labels{bb}{2*(tt-1)+1}{gg+1} = ['MaxRes_Spn' num2str(gg) 'Mrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                    gage_maximum{bb}{2*(tt-1)+1}(gg+1) = 1;
                    gage_theta{bb}{2*(tt-1)+1}(gg+1) = thetaMomentRotation(tt);
                    % create the minimum moment vector
                    gage_labels{bb}{2*(tt-1)+2}{gg+1} = ['MinRes_Spn' num2str(gg) 'Mrb' num2str(bb) '_' num2str(thetaMomentRotation(tt)) 'deg'];
                    gage_maximum{bb}{2*(tt-1)+2}(gg+1) = -1;
                    gage_theta{bb}{2*(tt-1)+2}(gg+1) = thetaMomentRotation(tt);
                end
            end
        end
    end
    % save the information into the fast_gage structure: fast_gage.labels{bb}{tt}
    fast_gage.labels = gage_labels;
    fast_gage.maximum = gage_maximum;
    fast_gage.theta = gage_theta;