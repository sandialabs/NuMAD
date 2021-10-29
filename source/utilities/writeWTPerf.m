function writeWTPerf(wt,filename,wtperfVersion)
% WRITEWTPERF  Write a WTPerf input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeWTPerf(wt,'filename',[wtperfVersion])
%      Works with WT_Perf v3.10 and v3.05.00a
%
%      wt = WTPerf data structure; view default structure by reading a
%           WTPerf input file using readWTPerf()
%      filename = filename of WTPerf file to be created.
%      wtperfVersion = (optional) select between 'v310' and
%                      'v3.05.00' formats
%

if ~exist('wtperfVersion','var') || isempty(wtperfVersion)
    switch wt.originalVersion
        case {'v310','v300','v3.10'}
            wtperfVersion = 'v3.10';  % flag to switch the output file format
        case {'v3.05','v3.05.00'}
            wtperfVersion = 'v3.05.00';  % flag to switch the output file format
        otherwise
            wtperfVersion = 'v3.05.00';  % flag to switch the output file format
    end
end

switch wtperfVersion
    case {'v310','v300','v3.10'}
        WTPver = 'v3.10';
    case {'v3.05','v3.05.00'}
        WTPver = 'v3.05.00a';
    otherwise
        error('wtperfVersion "%s" not recognized',wtperfVersion);
end

fw = fileWriter(filename);  % create a file writer object
fw.fmt_num = '%-20g';  % set the format strings
fw.fmt_str = '%-20s';
fw.fmt_lst = '%-16g ';

fw.text('-----  WT_Perf Input File  -----------------------------------------------------');
fw.text(wt.title{1});
switch WTPver
    case 'v3.10'
        fw.text('Compatible with WT_Perf v3.10');
    case 'v3.05.00a'
        fw.text('Compatible with WT_Perf v3.05.00a');
end
fw.text('-----  Input Configuration  ----------------------------------------------------');
fw.str( wt.Echo,    'Echo:                      Echo input parameters to "echo.out"?');
fw.str( wt.DimenInp,'DimenInp:                  Turbine parameters are dimensional?');
fw.str( wt.Metric,  'Metric:                    Turbine parameters are Metric (MKS vs FPS)?');
fw.text('-----  Model Configuration  ----------------------------------------------------');
fw.num( wt.NumSect, 'NumSect:                   Number of circumferential sectors.');
switch WTPver
    case 'v3.10'
        if wt.MaxIter < 100
            % MaxIter works differently in v310 and v3.05.00, warn user
            warning('MaxIter may be too low for v310 (%s)',filename);
        end
        fw.num( wt.MaxIter, 'MaxIter:                   Max number of iterations for induction factor.');
    case 'v3.05.00a'
        if wt.MaxIter > 100
            % MaxIter works differently in v310 and v3.05.00, warn user
            warning('MaxIter may be too high for v3.05.00 (%s)',filename);
        end
        fw.num( wt.MaxIter, 'MaxIter:                   Max number of iterations for Newton''s method to find induction factor.');
        fw.num( wt.NSplit,  'NSplit:                    Max number of splits for binary search method.');
end
fw.num( wt.ATol,    'ATol:                      Error tolerance for induction iteration.');
fw.num( wt.SWTol,   'SWTol:                     Error tolerance for skewed-wake iteration.');
fw.text('-----  Algorithm Configuration  ------------------------------------------------');
fw.str( wt.TipLoss, 'TipLoss:                   Use the Prandtl tip-loss model?');
fw.str( wt.HubLoss, 'HubLoss:                   Use the Prandtl hub-loss model?');
fw.str( wt.Swirl,   'Swirl:                     Include Swirl effects?');
fw.str( wt.SkewWake,'SkewWake:                  Apply skewed-wake correction?');
switch WTPver
    case 'v3.10'
        fw.str( wt.AdvBrake,'AdvBrake:                  Use the advanced brake-state model?');
        fw.str( wt.IndProp, 'IndProp:                   Use PROP-PC instead of PROPX induction algorithm?');
    case 'v3.05.00a'
        fw.str( wt.IndType, 'IndType:                   Use BEM induction algorithm?');
end
fw.str( wt.AIDrag,  'AIDrag:                    Use the drag term in the axial induction calculation.');
fw.str( wt.TIDrag,  'TIDrag:                    Use the drag term in the tangential induction calculation.');
switch WTPver
    case 'v3.05.00a'
        fw.str( wt.TISingularity,'TISingularity:             Use the singularity avoidance method in the tangential-induction calculation?');
        fw.str( wt.DAWT,         'DAWT:                      Run Diffuser Augmented Water Turbine Analysis? [feature not implimented yet]');
        fw.str( wt.Cavitation,   'Cavitation:                Run cavitation check? if cavitation, output sevens, check 12 oclock azimuth');
        fw.text('-----  Cavitation Model  -------------------------------------------------------');
        fw.num( wt.PressAtm,  'PressAtm:                  Air Atmospheric Pressure, Pa units, absolute');
        fw.num( wt.PressVapor,'PressVapor:                Vapor Pressure of Water, Pa units, absolute');
        fw.num( wt.CavSF,     'CavSF:                     Cavitation safety factor');
        fw.num( wt.WatDepth,  'WatDepth:                  Depth from water free surface to mudline (tower base)');
end
fw.text('-----  Turbine Data  -----------------------------------------------------------');
fw.num( wt.NumBlade,'NumBlade:                  Number of blades.');
fw.num( wt.RotorRad,'RotorRad:                  Rotor radius [length].');
fw.num( wt.HubRad,  'HubRad:                    Hub radius [length or div by radius].');
fw.num( wt.PreCone, 'PreCone:                   Precone angle, positive downwind [deg].');
fw.num( wt.Tilt,    'Tilt:                      Shaft tilt [deg].');
fw.num( wt.Yaw,     'Yaw:                       Yaw error [deg].');
fw.num( wt.HubHt,   'HubHt:                     Hub height [length or div by radius].');
fw.num( wt.NumSeg,  'NumSeg:                    Number of blade segments (entire rotor radius).');
fw.text('   RElm   Twist     Chord   AFfile  PrntElem');

for k=1:wt.NumSeg
    fprintf(fw.fileID,' %6g' ,wt.seg.RElm(k));
    fprintf(fw.fileID,' %7g' ,wt.seg.Twist(k));
    fprintf(fw.fileID,' %10g',wt.seg.Chord(k));
    fprintf(fw.fileID,' %5d  '  ,wt.seg.AFfile(k));
    fprintf(fw.fileID,'    %s',wt.seg.PrntElem{k});
    fprintf(fw.fileID,'\n');
end

fw.text('-----  Aerodynamic Data  -------------------------------------------------------');
fw.num( wt.Rho,     'Rho:                 Air density [mass/volume].');
fw.num( wt.KinVisc, 'KinVisc:             Kinematic air viscosity');
fw.num( wt.ShearExp,'ShearExp:            Wind shear exponent (1/7 law = 0.143).');
fw.str( wt.UseCm,   'UseCm                Are Cm data included in the airfoil tables?');
switch WTPver
    case 'v3.05.00a'
        fw.str( wt.UseCpmin,'UseCpmin:            Are Cp,min data included in the airfoil tables?');
end
fw.num( wt.NumAF,   'NumAF:               Number of airfoil files.');
fw.str( wt.AF_File{1},'AF_File:             List of NumAF airfoil files.');

for i=2:wt.NumAF
    fprintf(fw.fileID,'%s\n',wt.AF_File{i});
end

fw.text('-----  Output Configuration  ---------------------------------------------------');
switch WTPver
    case 'v3.05.00a'
        fw.str( wt.UnfPower, 'UnfPower:                  Write parametric power to an unformatted file?');
end
fw.str( wt.TabDel,   'TabDel:                    Make output tab-delimited (fixed-width otherwise).');
switch WTPver
    case 'v3.05.00a'
        fw.num( wt.ConvFlag, 'ConvFlag:                  For non-converging cases, 0 to output the result, 1 to output nines, 2 to output NaN (safest).');
        fw.str( wt.Beep,     'Beep:                      Beep if errors occur.');
end
fw.str( wt.KFact,    'KFact:                     Output dimensional parameters in K (e.g., kN instead on N)');
fw.str( wt.WriteBED, 'WriteBED:                  Write out blade element data to "<rootname>.bed"?');
fw.str( wt.InputTSR, 'InputTSR:                  Input speeds as TSRs?');
switch WTPver
    case 'v3.05.00a'
        fw.str( wt.OutMaxCp,'OutMaxCp:                  Output conditions for the maximum Cp?');
end
fw.str( wt.SpdUnits,'SpdUnits:                  Wind-speed units (mps, fps, mph).');

fw.text('-----  Combined-Case Analysis  -------------------------------------------------');
fw.num( wt.NumCases,'NumCases:                  Number of cases to run.  Enter zero for parametric analysis.');
fw.text('WS or TSR   RotSpd   Pitch                      Remove following block of lines if NumCases is zero.');

if wt.NumCases>0
    for k=1:wt.NumCases
        fprintf(fw.fileID,' %8g',wt.case{k}.WSorTSR);
        fprintf(fw.fileID,' %8g',wt.case{k}.RotSpd);
        fprintf(fw.fileID,' %8g',wt.case{k}.Pitch);
        fprintf(fw.fileID,'\n');
    end
end

fw.text('-----  Parametric Analysis (Ignored if NumCases > 0 )  -------------------------');
if wt.NumCases==0
    fw.num( wt.par.ParRow,'ParRow:                    Row parameter    (1-rpm, 2-pitch, 3-tsr/speed).');
    fw.num( wt.par.ParCol,'ParCol:                    Column parameter (1-rpm, 2-pitch, 3-tsr/speed).');
    fw.num( wt.par.ParTab,'ParTab:                    Table parameter  (1-rpm, 2-pitch, 3-tsr/speed).');
    fw.str( wt.par.OutPwr,'OutPwr:                    Request output of rotor power?');
    fw.str( wt.par.OutCp ,'OutCp:                     Request output of Cp?');
    fw.str( wt.par.OutTrq,'OutTrq:                    Request output of shaft torque?');
    fw.str( wt.par.OutFlp,'OutFlp:                    Request output of flap bending moment?');
    fw.str( wt.par.OutThr,'OutThr:                    Request output of rotor thrust?');
    fw.csv( wt.par.Pit   ,'PitSt, PitEnd, PitDel:     First, last, delta blade pitch (deg).');
    fw.csv( wt.par.Omg   ,'OmgSt, OmgEnd, OmgDel:     First, last, delta rotor speed (rpm).');
    fw.csv( wt.par.Spd   ,'SpdSt, SpdEnd, SpdDel:     First, last, delta speeds.');
end

end

