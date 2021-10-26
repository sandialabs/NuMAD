function writeTurbSim(tsim,output_file)
%WRITETURBSIM  Write a TurbSim input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeTurbSim(tsim,'file_name') 
%      Works with TurbSim v1.50
%
%      tsim = TurbSim data structure; view default structure by reading a 
%           input file using readTurbSim()
%      output_file = filename of file to be created.


if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,tsim.title,'');
wip(fid,'','');
wip(fid,'','---------Runtime Options-----------------------------------');
wip(fid,tsim.RandSeed1,'RandSeed1       - First random seed  (-2147483648 to 2147483647)');
wip(fid,tsim.RandSeed2,'RandSeed2       - Second random seed (-2147483648 to 2147483647) for intrinsic pRNG, or an alternative pRNG: "RanLux" or "RNSNLW"');
wip(fid,tsim.WrBHHTP,  'WrBHHTP         - Output hub-height turbulence parameters in GenPro-binary form?  (Generates RootName.bin)');
wip(fid,tsim.WrFHHTP,  'WrFHHTP         - Output hub-height turbulence parameters in formatted form?  (Generates RootName.dat)');
wip(fid,tsim.WrADHH,   'WrADHH          - Output hub-height time-series data in AeroDyn form?  (Generates RootName.hh)');
wip(fid,tsim.WrADFF,   'WrADFF          - Output full-field time-series data in TurbSim/AeroDyn form? (Generates RootName.bts)');
wip(fid,tsim.WrBLFF,   'WrBLFF          - Output full-field time-series data in BLADED/AeroDyn form?  (Generates RootName.wnd)');
wip(fid,tsim.WrADTWR,  'WrADTWR         - Output tower time-series data? (Generates RootName.twr)');
wip(fid,tsim.WrFMTFF,  'WrFMTFF         - Output full-field time-series data in formatted (readable) form?  (Generates RootName.u, RootName.v, RootName.w)');
wip(fid,tsim.WrACT,    'WrACT           - Output coherent turbulence time steps in AeroDyn form? (Generates RootName.cts)');
wip(fid,tsim.Clockwise,'Clockwise       - Clockwise rotation looking downwind? (used only for full-field binary files - not necessary for AeroDyn)');
switch tsim.ver
    case 'v1.3'
    wip(fid,tsim.ScaleIEC, 'ScaleIEC        - Scale hub-height IEC turbulence to target TI?');
    otherwise
    wip(fid,tsim.ScaleIEC, 'ScaleIEC        - Scale IEC turbulence models to exact target standard deviation? [0=no additional scaling; 1=use hub scale uniformly; 2=use individual scales]');
end
wip(fid,'','');
wip(fid,'','--------Turbine/Model Specifications-----------------------');
wip(fid,tsim.NumGrid_Z,   'NumGrid_Z       - Vertical grid-point matrix dimension');
wip(fid,tsim.NumGrid_Y,   'NumGrid_Y       - Horizontal grid-point matrix dimension');
wip(fid,tsim.TimeStep,    'TimeStep        - Time step [seconds]');
wip(fid,tsim.AnalysisTime,'AnalysisTime    - Length of analysis time series [seconds] (program will add time if necessary: AnalysisTime = MAX(AnalysisTime, UsableTime+GridWidth/MeanHHWS) )');
wip(fid,tsim.UsableTime,  'UsableTime      - Usable length of output time series [seconds] (program will add GridWidth/MeanHHWS seconds)');
wip(fid,tsim.HubHt,       'HubHt           - Hub height [m] (should be > 0.5*GridHeight)');
wip(fid,tsim.GridHeight,  'GridHeight      - Grid height [m]');
wip(fid,tsim.GridWidth,   'GridWidth       - Grid width [m] (should be >= 2*(RotorRadius+ShaftLength))');
wip(fid,tsim.VFlowAng,    'VFlowAng        - Vertical mean flow (uptilt) angle [degrees]');
wip(fid,tsim.HFlowAng,    'HFlowAng        - Horizontal mean flow (skew) angle [degrees]');
wip(fid,'','');
wip(fid,'','--------Meteorological Boundary Conditions-------------------');
wip(fid,tsim.TurbModel,      'TurbModel       - Turbulence model ("IECKAI"=Kaimal, "IECVKM"=von Karman, "GP_LLJ", "NWTCUP", "SMOOTH", "WF_UPW", "WF_07D", "WF_14D", or "NONE")');
wip(fid,tsim.IECstandard,    'IECstandard     - Number of IEC 61400-x standard (x=1,2, or 3 with optional 61400-1 edition number (i.e. "1-Ed2") )');
wip(fid,tsim.IECturbc,       'IECturbc        - IEC turbulence characteristic ("A", "B", "C" or the turbulence intensity in percent) ("KHTEST" option with NWTCUP, not used for other models)');
wip(fid,tsim.IEC_WindType,   'IEC_WindType    - IEC turbulence type ("NTM"=normal, "xETM"=extreme turbulence, "xEWM1"=extreme 1-year wind, "xEWM50"=extreme 50-year wind, where x=wind turbine class 1, 2, or 3)');
wip(fid,tsim.ETMc,           'ETMc            - IEC ETM "c" parameter [m/s] (or "default")');
wip(fid,tsim.WindProfileType,'WindProfileType - Wind profile type ("JET"=Low-level jet,"LOG"=Logarithmic,"PL"=Power law, or "default")');
wip(fid,tsim.RefHt,          'RefHt           - Height of the reference wind speed [m]');
wip(fid,tsim.URef,           'URef            - Mean (total) wind speed at the reference height [m/s]');
wip(fid,tsim.ZJetMax,        'ZJetMax         - Jet height [m] (used only for JET wind profile, valid 70-490 m)');
wip(fid,tsim.PLExp,          'PLExp           - Power law exponent  (or "default")');
wip(fid,tsim.Z0,             'Z0              - Surface roughness length [m] (or "default")');
wip(fid,'','');
wip(fid,'','--------Non-IEC Meteorological Boundary Conditions------------');
wip(fid,tsim.Latitude,'Latitude        - Site latitude [degrees] (or "default")');
wip(fid,tsim.RICH_NO, 'RICH_NO         - Gradient Richardson number');
wip(fid,tsim.UStar,   'UStar           - Friction or shear velocity [m/s] (or "default")');
wip(fid,tsim.ZI,      'ZI              - Mixing layer depth [m] (or "default")');
wip(fid,tsim.PC_UW,   'PC_UW           - Mean hub u''w'' Reynolds stress (or "default" or "none")');
wip(fid,tsim.PC_UV,   'PC_UV           - Mean hub u''v'' Reynolds stress (or "default" or "none")');
wip(fid,tsim.PC_VW,   'PC_VW           - Mean hub v''w'' Reynolds stress (or "default" or "none")');
wip(fid,tsim.IncDec1, 'IncDec1         - u-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")');
wip(fid,tsim.IncDec2, 'IncDec2         - v-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")');
wip(fid,tsim.IncDec3, 'IncDec3         - w-component coherence parameters (e.g. "10.0  0.3e-3" in quotes) (or "default")');
wip(fid,tsim.CohExp,  'CohExp          - Coherence exponent (or "default")');
wip(fid,'','');
wip(fid,'','--------Coherent Turbulence Scaling Parameters-------------------');
wip(fid,tsim.CTEventPath,'CTEventPath     - Name of the path where event data files are located');
wip(fid,tsim.CTEventFile,'CTEventFile     - Type of event files ("random", "les" or "dns")');
wip(fid,tsim.Randomize,  'Randomize       - Randomize disturbance scale and location? (true/false)');
wip(fid,tsim.DistScl,    'DistScl         - Disturbance scale (ratio of dataset height to rotor disk).');
wip(fid,tsim.CTLy,       'CTLy            - Fractional location of tower centerline from right (looking downwind) to left side of the dataset.');
wip(fid,tsim.CTLz,       'CTLz            - Fractional location of hub height from the bottom of the dataset.');
wip(fid,tsim.CTStartTime,'CTStartTime     - Minimum start time for coherent structures in RootName.cts [seconds]');
wip(fid,'','');
wip(fid,'','==================================================');
wip(fid,'','NOTE: Do not add or remove any lines in this file!');
wip(fid,'','==================================================');

if fid~=1, fclose(fid); end
end


function wip(fid,param,descrip)
% write input file parameter
if isempty(param)
    % do nothing if param empty
elseif ischar(param)
    fprintf(fid,'%-16s ',param);  %output string
elseif isfloat(param)
    if numel(param)==1
        fprintf(fid,'%-16g ',param);  %output single number
    else
        str = sprintf(' %g',param);        %create list of numbers
        str = sprintf('"%s"',str(2:end));  %quote the list of numbers
        fprintf(fid,'%-16s ',str);          %output the quoted list
    end
end
fprintf(fid,'%s\n',descrip);
end

