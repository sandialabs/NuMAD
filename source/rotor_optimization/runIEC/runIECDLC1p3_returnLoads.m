function [results]=runIECDLC1p3_returnLoads(params,output)

CaseName='IECDLC1p3ETM';

% ====================== MAIN PROCESS - SIMULATIONS =======================
% create a counter for parallel runs
for ii = 1:length(params.ws)
    ctr.w((ii-1)*params.numSeeds+1:ii*params.numSeeds) = params.ws(ii);
    ctr.s((ii-1)*params.numSeeds+1:ii*params.numSeeds) = 1:params.numSeeds;
end
ctr.w1 = ctr.w; ctr.s1 = ctr.s;
yawMat = repmat(params.yaw, params.numSeeds*length(params.ws), 1);
ctr.y = reshape(yawMat, 1, numel(yawMat));
ctr.w = repmat(ctr.w, 1, length(params.yaw));
ctr.s = repmat(ctr.s, 1, length(params.yaw));
% loop through "wind direction" as well
dirMat = repmat(params.wd, numel(yawMat), 1);
ctr.d = reshape(dirMat, 1, numel(dirMat));
ctr.y = repmat(ctr.y, 1, length(params.wd));
ctr.w = repmat(ctr.w, 1, length(params.wd));
ctr.s = repmat(ctr.s, 1, length(params.wd));


% Create a cell with the FAST output filenames
outName={};
for yaw = params.yaw
    for w=params.ws
        for s=1:params.numSeeds
            for d=params.wd
                outName{end+1}=[params.parDir 'out/' CaseName '_yaw' num2str(yaw) '_' ...
                    num2str(w) 'mps_seed' num2str(s) '_dir' num2str(d) '.out'];
            end
        end
    end
end


% read in the FAST output files and calculate channels as desired
temp = cell(1,length(outName));
parfor cc = 1:length(outName)
    out = loadFASTOutData(outName{cc});
    
    % 1 - compare prediction for estimated force from known moment
    Mxb1_meas = out.data(:,strcmp(out.list,'RootMxb1'));
    Myb1_meas = out.data(:,strcmp(out.list,'RootMyb1'));
    Mxb2_meas = out.data(:,strcmp(out.list,'RootMxb2'));
    Myb2_meas = out.data(:,strcmp(out.list,'RootMyb2'));
    Mxb3_meas = out.data(:,strcmp(out.list,'RootMxb3'));
    Myb3_meas = out.data(:,strcmp(out.list,'RootMyb3'));
    Mx_meas = [Mxb1_meas; Mxb2_meas; Mxb3_meas];
    My_meas = [Myb1_meas; Myb2_meas; Myb3_meas];
    
    Fxb1_meas = out.data(:,strcmp(out.list,'RootFxb1'));
    Fyb1_meas = out.data(:,strcmp(out.list,'RootFyb1'));
    Fxb2_meas = out.data(:,strcmp(out.list,'RootFxb2'));
    Fyb2_meas = out.data(:,strcmp(out.list,'RootFyb2'));
    Fxb3_meas = out.data(:,strcmp(out.list,'RootFxb3'));
    Fyb3_meas = out.data(:,strcmp(out.list,'RootFyb3'));
    Fx_meas = [Fxb1_meas; Fxb2_meas; Fxb3_meas];
    Fy_meas = [Fyb1_meas; Fyb2_meas; Fyb3_meas];
    
    Lblade = 13;
    Fx_hat = 2*My_meas/Lblade;
    Fy_hat = -2*Mx_meas/Lblade;
    % calculate errors
    temp{cc}.Fx_Err_kN = Fx_hat-Fx_meas;
    temp{cc}.Fy_Err_kN = Fy_hat-Fy_meas;
    temp{cc}.Fx_PercErr = (Fx_hat-Fx_meas)./Fx_meas.*100;
    temp{cc}.Fy_PercErr = (Fy_hat-Fy_meas)./Fy_meas.*100;
    
    if 1
        figure
        subplot 211
        plot(Fx_meas); hold on
        plot(Fx_hat,'x')
        legend('Fx','\overbar{F_x}')
        ylabel('Thrust Direction')
        subplot 212
        plot(Fx_meas); hold on
        plot(Fx_hat,'x')
        legend('Fy','\overbar{F_y}')
        xlabel('Torque Direction')
        
        figure
        subplot 211
        histogram(temp{cc}.Fx_Err_kN)
        ylabel('Fx Error (kN)')
        subplot 212
        histogram(temp{cc}.Fy_Err_kN)
        ylabel('Fy Error (kN)')
        
        figure
        subplot 211
        histogram(temp{cc}.Fx_PercErr)
        ylabel('Fx % Error')
        subplot 212
        histogram(temp{cc}.Fy_PercErr)
        ylabel('Fy % Error')
    end
end


% create the final max/min/ampl values from each of the simulations
names = fieldnames(temp{1});
for ii = 1:length(names) % loop through the channel variable
    tempVec = [];
    % save the data values in a temporary vector
    for jj = 1:length(temp)        
        tempVec = [tempVec; temp{jj}.(names{ii})];
    end
    
    results.(names{ii}) = tempVec;
    
end




end