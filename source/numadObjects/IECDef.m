%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Part of the SNL NuMAD Toolbox                    
%  Developed by Sandia National Laboratories Wind Energy Technologies 
%              See license.txt for disclaimer information             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef IECDef < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ``IECDef``  A class definition for IEC input parameters
%
%   Usage: 
%     params = IECDef();
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    properties (SetAccess = 'public', GetAccess = 'public')        
        fstfn= 'NOT DEFINED';   	% FAST path, Default = "NOT DEFINED"
        numadfn='NOT DEFINED';      % NuMAD path, including extension, Default = "NOT DEFINED"
        fastsim = 'fast';           % Fast simulation, Options: 'fast', 'fast simulink', 'adams', Default: 'fast'
        simulinkModel = 'NOT DEFINED';          % Simulink model file, Default = "NOT DEFINED"
        simulinkModelFolder = 'NOT DEFINED';    % Simulink model directory, Default = "NOT DEFINED"
        operatingPoints=[0 0 0];    % Operating Points: [Cutin RatedSpeed CutOut], Defaults = [0 0 0]
        ws=[3:2:25];                % Range of mean wind speeds for turbulent simulations (m/s)
        wd=[180];                   % Range of "wind direction" bias for look-up table simulations - programmed as yaw position (units)
        yaw = [0];                  % Intentional yaw misalignment, degrees (for DLC 1.1)
        momentMaxRotation = 45;     % Angular discretization for coordinate rotation and maxima moment calculation (used for fatigue and ultimate) (deg)
        ratedSpeed=12;              % Rated speed (rpm)
        lin=10;                     % Range of steady wind speeds for linearizations, Default = 10
        sf_fat=0;                   % Total fatigue safety factor, Default = 0
        sf_uts=0;                   % Total ultimate strength safety factor, Default = 0
        sf_tow=0;                   % Total tower clearance safety factor, Default = 0
        numSeeds=6;                 % Number of seeds - number of 10-minute simulations - for turbulent simulations
        delay=29;                   % Simulation delay, discardedsimulated data at the beginning of each simulation (turbulent and otherwise), Default = 29 (s)
        SimTime=600;                % Total simulation time, Default = 600 (s)
        NumGrid=10              	% Number of grid points in turbsim 4-D wind field, Default = 10
        Class=3;                    % Turbine Class, Options: 1, 2, 3, Default = 3
        TurbClass='C';              % Turbulence Class, Options: 'A', 'B', 'C', Default = 'C'
        designLife = 30;            % Years of life 
        BldGagNd=[1,2,3,4,5,6,7];   % Vector of length 7; blade gage nodes (corresponding to aerodyn nodes) for moment (strain) gages in FAST computations 
        fatigueCriterion = 'Shifted Goodman';   % Fatigue Criteria, Options: 'Shifted Goodman'
    end
    
    properties (SetAccess=private)
        avgws=0                     % Average wind speed, Default = 0 (m/s)
    end
    
    properties (Hidden, SetAccess=private)
    end
    
    properties (Hidden)
    end
    
    methods
        function obj = IECDef()
            % This method initilizes ``IECDef`` and creates a
            % ``params`` object.          
            % Returns
            % ------------
            %     param : obj
            %         IECDef object    
            %             
            if nargin > 0
                
            end

% Kelley: Need to move this so that it's called once params.Class is defined             
            switch obj.Class
                case 1
                    obj.avgws=0.2*50; % m/s, average wind speed of IEC Class I site (Vref=50m/s); IEC Section 6.3.1.1 Eqn (9)
                case 2
                    obj.avgws=0.2*42.5; % m/s, average wind speed of IEC Class II site (Vref=42.5m/s); IEC Section 6.3.1.1 Eqn (9)
                case 3
                    obj.avgws=0.2*37.5; % m/s, average wind speed of IEC Class III site (Vref=37.5m/s); IEC Section 6.3.1.1 Eqn (9)
            end    
            
        end
        
        function checkInputs(obj)
            % This method checks NuMAD user inputs and generates error
            % messages if parameters are not properly defined. 
            
            if ~isequal(obj.Class,1) && ~isequal(obj.Class,2) && ~isequal(obj.Class,3)
                    error('`params.Class` must be equal to 1, 2, or 3');
            end           
            
            if ~strcmp(obj.fastsim,'fast') && ~strcmp(obj.fastsim,'fast simulink') && ~strcmp(obj.fastsim,'adams')
                    error('`params.fastsim` must be equal to "fast", "fast simulink", "adams"');
            end   
            
        end
        
        
        
    end            
end