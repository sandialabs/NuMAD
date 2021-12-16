%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef IECDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``IECDef``  A class definition for IEC input parameters
%
% Example: 
%
%	``IEC = IECDef(inputFile);``
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = 'public', GetAccess = 'public')  
        avgws                       % Average wind speed, Default = [] (m/s)
        BldGagNd=[1,2,3,4,5,6,7];   % Vector of length 7; blade gage nodes (corresponding to aerodyn nodes) for moment (strain) gages in FAST computations 
        bladeGageCoordinateRotation
        bladeGageLabels
        bladeGageLabels_MLx
        bladeGageLabels_MLy
        bladeGageLabels_MLz
        bladeGageLabels_FLx
        bladeGageLabels_FLy
        bladeGageLabels_FLz        
        combinedStrain = 1;         % perform fatigue calculations on combined bending and normal strain for spar, Default = 1
        seeds                       % random seeds
        simulate                    % On/Off flag, Options: 1 = call FAST and perform simulations, 0 = process existing data, Default = []
        Class=3;                    % Turbine Class, Options: 1, 2, 3, Default = 3
        delay=29;                   % Simulation delay, discardedsimulated data at the beginning of each simulation (turbulent and otherwise), Default = 29 (s)
        designLife = 30;            % Years of life 
        fastsim = 'fast';           % Fast simulation, Options: 'fast', 'fast simulink', 'adams', Default: 'fast'
        fatigueCriterion = 'Shifted Goodman';   % Fatigue Criteria, Options: 'Shifted Goodman'
        fatigueStress = 'Equivalent';  % stress to use for fatigue failure
        fstfn= 'NOT DEFINED';   	% FAST path, Default = "NOT DEFINED"
        fullLoads = 1;              % Perform full loads analysis, Options: 1 = On, 0 = Off, Default = 1
        gageSetCase = 'set1'        % Gage set case, Options= 'set1', 'set 2', Default = 'set1'
        lin=10;                     % Range of steady wind speeds for linearizations, Default = 10
        matData = [];               % List of material objects, with properties relavent to IEC fatigue analysis
        momentMaxRotation = 45;     % Angular discretization for coordinate rotation and maxima moment calculation (used for fatigue and ultimate) (deg)
        nondimFlag = true;          % Flag indicating non-dimensionality in plotCampbell() function
        numadfn='NOT DEFINED';      % NuMAD path, including extension, Default = "NOT DEFINED"
        NumGrid=10              	% Number of grid points in turbsim 4-D wind field, Default = 10
        numSeeds=6;                 % Number of seeds - number of 10-minute simulations - for turbulent simulations
        operatingPoints=[0 0 0];    % Operating Points: [Cutin RatedSpeed CutOut], Defaults = [0 0 0]
        parDir = '';                % Directory, Default = '';             
        ratedSpeed=12;              % Rated speed (rpm)
        sf_fat=0;                   % Total fatigue safety factor, Default = 0
        sf_uts=0;                   % Total ultimate strength safety factor, Default = 0
        sf_tow=0;                   % Total tower clearance safety factor, Default = 0
        SimTime=600;                % Total simulation time, for a given run Default = 600 (s)
        simtime=0;                  % Total simulation time for a set of seeds, for rain-cycle counting
        simulinkModel = 'NOT DEFINED';          % Simulink model file, Default = "NOT DEFINED"
        simulinkModelFolder = 'NOT DEFINED';    % Simulink model directory, Default = "NOT DEFINED"
        TurbClass='C';              % Turbulence Class, Options: 'A', 'B', 'C', Default = 'C'
        ws=[3:2:25];                % Range of mean wind speeds for turbulent simulations (m/s)
        wd=[180];                   % Range of "wind direction" bias for look-up table simulations - programmed as yaw position (units)
        yaw = [0];                  % Intentional yaw misalignment, degrees (for DLC 1.1)        
    end
    
    properties (SetAccess=private)                  
    end
    
    properties (Hidden, SetAccess=private)
    end
    
    properties (Hidden)
    end
    
    methods
        function obj = IECDef(fileName)
            % This method initilizes ``IECDef`` and creates a
            % ``IEC`` object.          
            % Returns
            % ------------
            %     IEC : obj
            %         IECDef object    
            %             
            
            if(nargin > 0)
                obj.readInput(fileName);
            end
        end
        
        function readInput(obj,fileName)
            fid = fopen(fileName);
            fl = fgetl(fid);
            currentMat = 0;
            while(fl ~= -1)
                if(contains(fl,'%'))
                    strList = split(fl,'%');
                    fl = strList{1};
                end
                fl = strip(fl);
                if(~isempty(fl))
                    strList = split(fl);
                    varName = strList{1};
                    strList = split(fl,varName);
                    varName = lower(varName);
                    inDat = strip(strList{2});
                    numList = str2num(inDat);
                    if(isempty(numList))
                        varType = 'char';
                    else
                        varType = 'num';
                    end
                    found = 0;
                    switch varName
                        case 'class'
                            if(strcmp(varType,'num'))
                                obj.Class = numList;
                                found = 1;
                            end
                        case 'designlife'
                            if(strcmp(varType,'num'))
                                obj.designLife = numList;
                                found = 1;
                            end
                        case 'fastfile'
                            if(strcmp(varType,'char'))
                                obj.fstfn = inDat;
                                found = 1;
                            end
                        case 'fastsimmethod'
                            if(strcmp(varType,'char'))
                                obj.fastsim = inDat;
                                found = 1;
                            end
                        case 'fatiguecriterion'
                            if(strcmp(varType,'char'))
                                obj.fatigueCriterion = inDat;
                                found = 1;
                            end
                        case 'fatiguesafetyfactor'
                            if(strcmp(varType,'num'))
                                obj.sf_fat = numList;
                                found = 1;
                            end
                        case 'fullloads'
                            if(strcmp(varType,'num'))
                                obj.fullLoads = numList;
                                found = 1;
                            end
                        case 'gagenodes'
                            if(strcmp(varType,'num'))
                                obj.BldGagNd = numList;
                                found = 1;
                            end
                        case 'gageSetCase'
                            if(strcmp(varType,'char'))
                                obj.gageSetCase = inDat;
                                found = 1;
                            end
                        case 'linearizationws'
                            if(strcmp(varType,'num'))
                                obj.lin = numList;
                                found = 1;
                            end
                        case 'maxloadangleincrement'
                            if(strcmp(varType,'num'))
                                obj.momentMaxRotation = numList;
                                found = 1;
                            end
                        case 'nondimensionalflag'
                            if(strcmp(varType,'char'))
                                if(strcmp(inDat,'false'))
                                    obj.nondimFlag = false;
                                else
                                    obj.nondimFlad = true;
                                end
                                found = 1;
                            end
                        case 'numadfile'
                            if(strcmp(varType,'char'))
                                obj.numadfn = inDat;
                                found = 1;
                            end
                        case 'numgridpoints'
                            if(strcmp(varType,'num'))
                                obj.NumGrid = numList;
                                found = 1;
                            end
                        case 'numseeds'
                            if(strcmp(varType,'num'))
                                obj.numSeeds = numList;
                                found = 1;
                            end
                        case 'operatingpoints'
                            if(strcmp(varType,'num'))
                                obj.operatingPoints = numList;
                                found = 1;
                            end
                        case 'parentdirectory'
                            if(strcmp(varType,'char'))
                                obj.parDir = inDat;
                                found = 1;
                            end
                        case 'ratedspeed'
                            if(strcmp(varType,'num'))
                                obj.ratedSpeed = numList;
                                found = 1;
                            end
                        case 'simtime'
                            if(strcmp(varType,'num'))
                                obj.SimTime = numList;
                                found = 1;
                            end
                        case 'simulinkmodel'
                            if(strcmp(varType,'char'))
                                obj.simulinkModel = inDat;
                                found = 1;
                            end
                        case 'simulinkmodelfolder'
                            if(strcmp(varType,'char'))
                                obj.simulinkModelFolder = inDat;
                                found = 1;
                            end
                        case 'startupdelay'
                            if(strcmp(varType,'num'))
                                obj.delay = numList;
                                found = 1;
                            end
                        case 'towersafetyfactor'
                            if(strcmp(varType,'num'))
                                obj.sf_tow = numList;
                                found = 1;
                            end
                        case 'turbulenceclass'
                            if(strcmp(varType,'char'))
                                obj.TurbClass = inDat;
                                found = 1;
                            end
                        case 'ultimatesafetyfactor'
                            if(strcmp(varType,'num'))
                                obj.sf_uts = numList;
                                found = 1;
                            end
                        case 'winddirections'
                            if(strcmp(varType,'num'))
                                obj.wd = numList;
                                found = 1;
                            end
                        case 'windspeeds'
                            if(strcmp(varType,'num'))
                                obj.ws = numList;
                                found = 1;
                            end
                        case 'yaw'
                            if(strcmp(varType,'num'))
                                obj.yaw = numList;
                                found = 1;
                            end
                        case 'fatiguematerial:'
                            newMat = MaterialDef;
                            obj.matData = [obj.matData,newMat];
                            currentMat = currentMat + 1;
                            found = 1;
                        case 'name'
                            if(strcmp(varType,'char'))
                                obj.matData(currentMat).name = inDat;
                                found = 1;
                            end
                        case 'modulus'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).ex = numList;
                                found = 1;
                            end
                        case 'fatigueslope'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).m = numList;
                                found = 1;
                            end
                        case 'tensilestrength'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).uts = numList;
                                found = 1;
                            end
                        case 'compressivestrength'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).ucs = numList;
                                found = 1;
                            end
                        case 'strengthreductionfact'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).gamma_mf = numList;
                                found = 1;
                            end
                        case 'shorttermreductionfact'
                            if(strcmp(varType,'num'))
                                obj.matData(currentMat).gamma_ms = numList;
                                found = 1;
                            end
                    end
                    if(found == 0)
                        wrnmsg = ['The line in the IEC input file containing "' varName '" was not recognized.  Either the input name or value type entered are invalid.  The line will be ignored.'];
                        warning(wrnmsg);
                    end
                end
                fl = fgetl(fid);
                if(isempty(fl))
                    fl = ' ';
                end
            end
            fclose(fid);
        end
        
        function checkInputs(obj)
            % This method checks NuMAD user inputs and generates error
            % messages if parameters are not properly defined. 
            
            % Check obj.Class
            if ~isequal(obj.Class,1) && ~isequal(obj.Class,2) && ~isequal(obj.Class,3)
                    error('`IECDef.Class` must be equal to 1, 2, or 3');
            end           
            
            % Check obj.fastsim
            if ~strcmp(obj.fastsim,'fast') && ~strcmp(obj.fastsim,'fast simulink') && ~strcmp(obj.fastsim,'adams')
                    error('`IECDef.fastsim` must be equal to "fast", "fast simulink", "adams"');
            end   
            
        end
        
        function setAvgWindSpeed(obj)
            % This method sets the average wind speed based on the
            % obj.Class
            switch obj.Class
                case 1
                    obj.avgws=0.2*50; % m/s, average wind speed of IEC Class I site (Vref=50m/s); IEC Section 6.3.1.1 Eqn (9)
                case 2
                    obj.avgws=0.2*42.5; % m/s, average wind speed of IEC Class II site (Vref=42.5m/s); IEC Section 6.3.1.1 Eqn (9)
                case 3
                    obj.avgws=0.2*37.5; % m/s, average wind speed of IEC Class III site (Vref=37.5m/s); IEC Section 6.3.1.1 Eqn (9)
            end               
        end
        
        function setSimFlag(obj,simFlag)
            %This method sets the simulate flag
            obj.simulate=simFlag;
        end

        function setBladeGageCoordinateRotation(obj,simFlag)
            %This method sets the blade gage coordinate roation
            obj.bladeGageCoordinateRotation=simFlag;
        end        
        
        function setGageLabels(obj,fst,ad)
            % This method sets the blade gage labels
            for ss = 1:length(obj.BldGagNd)
                obj.bladeGageLabels_MLx{ss} = ['Spn' num2str(obj.BldGagNd(ss)) 'MLxb1']; %#ok<*AGROW>
                obj.bladeGageLabels_MLy{ss} = ['Spn' num2str(obj.BldGagNd(ss)) 'MLyb1'];
                obj.bladeGageLabels_MLz{ss} = '';
                obj.bladeGageLabels_FLx{ss} = '';
                obj.bladeGageLabels_FLy{ss} = '';
                obj.bladeGageLabels_FLz{ss} = '';
            end
            obj.bladeGageLabels = [obj.bladeGageLabels_MLx obj.bladeGageLabels_MLy];%...
                %obj.bladeGageLabels_MLz obj.bladeGageLabels_FLx...
                %obj.bladeGageLabels_FLy obj.bladeGageLabels_FLz];       
            % this script must be run twice for the first 9 gages, and then for the additional gages.
            switch obj.gageSetCase
                case 'set1'
                    setGag = 1:9; % *set 1*
                case 'set2'
                    setGag = 10:14; % *set 2*
            end                
            NBlGages = length(setGag);
            totalBladeGageNumber = 9;   % actual number of gages to be used
            BldGagMaxSpanLoc = 0.95;
            BldGagMinSpanLoc = 0;     % minimum span location of a strain gage (no gage here currently; RootMxb1...)
            drSpan = (BldGagMaxSpanLoc - 0) / totalBladeGageNumber;
            gagSpan = BldGagMinSpanLoc + drSpan : drSpan : BldGagMaxSpanLoc;
            bladeCoordinateSpanStations = (ad.RNodes - fst.TurbConf.HubRad) ./ (fst.TurbConf.TipRad-fst.TurbConf.HubRad);
            obj.BldGagNd = 0;
            obj.bladeGageLabels_MLx = {}; obj.bladeGageLabels_MLy = {}; obj.bladeGageLabels_MLz = {};
            obj.bladeGageLabels_FLx = {}; obj.bladeGageLabels_FLy = {}; obj.bladeGageLabels_FLz = {};
            obj.bladeGageLabels = {};
            for bb = 1%:3 %number of blades
                for ss = 1:length(setGag)
                    [~, rGagSpan(ss)] = min(abs(bladeCoordinateSpanStations - gagSpan(setGag(ss))));
                    obj.BldGagNd(ss) = rGagSpan(ss);
                    if bb == 1
                        obj.bladeGageLabels_MLx{ss} = ['Spn' num2str(ss) 'MLxb' num2str(bb)]; %#ok<*AGROW>
                        obj.bladeGageLabels_MLy{ss} = ['Spn' num2str(ss) 'MLyb' num2str(bb)];
                        obj.bladeGageLabels_MLz{ss} = ['Spn' num2str(ss) 'MLzb' num2str(bb)];
                        obj.bladeGageLabels_FLx{ss} = ['Spn' num2str(ss) 'FLxb' num2str(bb)];
                        obj.bladeGageLabels_FLy{ss} = ['Spn' num2str(ss) 'FLyb' num2str(bb)];
                        obj.bladeGageLabels_FLz{ss} = ['Spn' num2str(ss) 'FLzb' num2str(bb)];
                        rowName{ss} = ['gage' num2str(ss)];
                    else
                        obj.bladeGageLabels_MLx{end+1} = ['Spn' num2str(ss) 'MLxb' num2str(bb)]; %#ok<*AGROW>
                        obj.bladeGageLabels_MLy{end+1} = ['Spn' num2str(ss) 'MLyb' num2str(bb)];
                        obj.bladeGageLabels_MLz{end+1} = ['Spn' num2str(ss) 'MLzb' num2str(bb)];
                        obj.bladeGageLabels_FLx{end+1} = ['Spn' num2str(ss) 'FLxb' num2str(bb)];
                        obj.bladeGageLabels_FLy{end+1} = ['Spn' num2str(ss) 'FLyb' num2str(bb)];
                        obj.bladeGageLabels_FLz{end+1} = ['Spn' num2str(ss) 'FLzb' num2str(bb)];
                    end
                end
            end            
            obj.bladeGageLabels = [obj.bladeGageLabels_MLx obj.bladeGageLabels_MLy...
            obj.bladeGageLabels_MLz obj.bladeGageLabels_FLx...
            obj.bladeGageLabels_FLy obj.bladeGageLabels_FLz];
        
            tbl = array2table([gagSpan(setGag)' bladeCoordinateSpanStations(rGagSpan) rGagSpan']);
            tbl.Properties.VariableNames = {'GageSpanLocation_bladeCoordinates' 'ActualSpanLocations_bladeCoordinates' 'ADlocation'};
            tbl.Properties.RowNames = rowName;
            disp(tbl)
            disp('Press F5 to confirm the strain gage locations and proceed...')        
        end
        
                
        function runFullLoads(obj)
            % This methods runs the full loads analysis if specified
            if false(obj.fullLoads)
                % define the necessary FAST outputs
                fst.OutList=[{'WindVxi','WindVyi','WindVzi',...
                    'HorWindV','NacYawErr',...  % add these channels for JCB's yaw controller
                    'GenPwr','HSShftTq','HSSBrTq','GenSpeed','RotSpeed','TSR',...
                    'TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS',...
                    'NcIMUTAxs','NcIMUTAys',...
                    'BldPitch1',...
                    'OoPDefl1','OoPDefl2','OoPDefl3',...
                    'TipClrnc1','TipClrnc2','TipClrnc3',...
                    'TipDxb1','TipDxb2','TipDxb3'}...
                    obj.bladeGageLabels...
                    {'LSShftFxs','LSShftMxs',...
                    ...%'RootMxc1','RootMyc1',...
                    ...%'RootMxc2','RootMyc2',...
                    ...%'RootMxc3','RootMyc3' ...
                    'RootMxb1','RootMyb1','RootMzb1',...
                    'RootMxb2','RootMyb2','RootMzb2',...
                    'RootMxb3','RootMyb3','RootMzb3'}];   
            else
                % set up analysis for complete loads determination to ensure turbine
                % operation within structure loads envelope
                fst.OutList=[{'WindVxi','WindVyi','WindVzi',...
                    'HorWindV','NacYawErr',...  % add these channels for JCB's yaw controller
                    'GenPwr','HSShftTq','HSSBrTq','GenSpeed','RotSpeed','TSR',...
                    'TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS',...
                    'NcIMUTAxs','NcIMUTAys',...
                    'BldPitch1',...
                    'OoPDefl1','OoPDefl2','OoPDefl3',...
                    'TipClrnc1','TipClrnc2','TipClrnc3',...
                    'TipDxb1','TipDxb2','TipDxb3'}...
                    obj.bladeGageLabels...
                    ... SOME NEW CHANNELS TO ADD:
                    {...%'RootFxc1','RootFxc2','RootFxc3',...
                    ...%'RootFyc1','RootFyc2','RootFyc3',...
                    'RootFxb1','RootFxb2','RootFxb3',...
                    'RootFyb1','RootFyb2','RootFyb3',...
                    'RootFzb1','RootFzb2','RootFzb3',...
                    ...%'RootMxc1','RootMxc2','RootMxc3',...
                    ...%'RootMyc1','RootMyc2','RootMyc3',...
                    'RootMxb1','RootMxb2','RootMxb3',...
                    'RootMyb1','RootMyb2','RootMyb3',...
                    'RootMzb1','RootMzb2','RootMzb3',...
                    'YawBrFxn','YawBrFyn','YawBrFzn',...
                    'YawBrFxp','YawBrFyp',...
                    'YawBrMxn','YawBrMyn','YawBrMzn',...
                    'YawBrMxp','YawBrMyp',...
                    'RotThrust','LSShftTq',...
                    ...%'LSShftFya','LSShftFza',...    
                    ...%'LSSTipMya','LSSTipMza',...
                    'LSShftFxs','LSShftFys','LSShftFzs',...
                    'LSShftMxs','LSSTipMys','LSSTipMzs',...
                    'TwrBsFxt','TwrBsFyt','TwrBsFzt',...
                    'TwrBsMxt','TwrBsMyt','TwrBsMzt',...
                    }];
                disp(fst.OutList)
            end            
        
        end
        
        function setRandomSeeds(obj)
            % This method setd random seeds
            if ~exist([obj.parDir 'seeds.mat'],'file')
                seeds=randi(123456,1,obj.numSeeds);
                save([obj.parDir 'seeds.mat'],'seeds')
                obj.seeds = seeds;
            else
                load([obj.parDir 'seeds.mat']);
                % ble: added this check - seeds file can exist but not be correct.
                if length(seeds) ~= obj.numSeeds
                    clear seeds
                    seeds=randi(123456,1,obj.numSeeds);
                    save([obj.parDir 'seeds.mat'],'seeds')
                    obj.seeds = seeds;
                else
                    obj.seeds = seeds;
                end
                % ble: ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            end
            %obj.seeds = seeds; %ble: seeds needs to be saved this way for parallel operation.                                
        end
        
        function saveOuput(obj)
%             consider moving runIEC scripts to save output to this method
        end
        
    end            
end