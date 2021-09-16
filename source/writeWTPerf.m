function writeWTPerf(wt,output_file)
% WRITEWTPERF  Write a WTPerf input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeWTPerf(wt,'file_name')
%      Works with WT_Perf v3.10
%
%      wt = WTPerf data structure; view default structure by reading a 
%           WTPerf input file using readWTPerf()
%      output_file = filename of WTPerf file to be created.
%

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],'-----  WT_Perf Input File  -----------------------------------------------------');
wip(fid,wt.title{1},[]);
wip(fid,[],'Input file for WT Perf v3.10');
wip(fid,[],'-----  Input Configuration  ----------------------------------------------------');
wip(fid,wt.Echo,'Echo:                      Echo input parameters to "echo.out"?');
wip(fid,wt.DimenInp,'DimenInp:                  Turbine parameters are dimensional?');
wip(fid,wt.Metric,'Metric:                    Turbine parameters are Metric (MKS vs FPS)?');
wip(fid,[],'-----  Model Configuration  ----------------------------------------------------');
wip(fid,wt.NumSect,'NumSect:                   Number of circumferential sectors.');
wip(fid,wt.MaxIter,'MaxIter:                   Max number of iterations for induction factor.');
wip(fid,wt.ATol,'ATol:                      Error tolerance for induction iteration.');
wip(fid,wt.SWTol,'SWTol:                     Error tolerance for skewed-wake iteration.');
wip(fid,[],'-----  Algorithm Configuration  ------------------------------------------------');
wip(fid,wt.TipLoss,'TipLoss:                   Use the Prandtl tip-loss model?');
wip(fid,wt.HubLoss,'HubLoss:                   Use the Prandtl hub-loss model?');
wip(fid,wt.Swirl,'Swirl:                     Include Swirl effects?');
wip(fid,wt.SkewWake,'SkewWake:                  Apply skewed-wake correction?');
wip(fid,wt.AdvBrake,'AdvBrake:                  Use the advanced brake-state model?');
wip(fid,wt.IndProp,'IndProp:                   Use PROP-PC instead of PROPX induction algorithm?');
wip(fid,wt.AIDrag,'AIDrag:                    Use the drag term in the axial induction calculation.');
wip(fid,wt.TIDrag,'TIDrag:                    Use the drag term in the tangential induction calculation.');
wip(fid,[],'-----  Turbine Data  -----------------------------------------------------------');
wip(fid,wt.NumBlade,'NumBlade:                  Number of blades.');
wip(fid,wt.RotorRad,'RotorRad:                  Rotor radius [length].');
wip(fid,wt.HubRad,'HubRad:                    Hub radius [length or div by radius].');
wip(fid,wt.PreCone,'PreCone:                   Precone angle, positive downwind [deg].');
wip(fid,wt.Tilt,'Tilt:                      Shaft tilt [deg].');
wip(fid,wt.Yaw,'Yaw:                       Yaw error [deg].');
wip(fid,wt.HubHt,'HubHt:                     Hub height [length or div by radius].');
wip(fid,wt.NumSeg,'NumSeg:                    Number of blade segments (entire rotor radius).');
wip(fid,[],'   RElm   Twist     Chord   AFfile  PrntElem');

ElmTable = {wt.seg.RElm, wt.seg.Twist, wt.seg.Chord, wt.seg.AFfile, wt.seg.PrntElem};
wiptbl(fid,{'%8.5f','%7.2f','%7.3f','%4d','     %s'},ElmTable,wt.NumSeg);

wip(fid,[],'-----  Aerodynamic Data  -------------------------------------------------------');
wip(fid,wt.Rho,'Rho:                 Air density [mass/volume].');
wip(fid,wt.KinVisc,'KinVisc:             Kinematic air viscosity');
wip(fid,wt.ShearExp,'ShearExp:            Wind shear exponent (1/7 law = 0.143).');
wip(fid,wt.UseCm,'UseCm                Are Cm data included in the airfoil tables?');
% wip(fid,wt.UseCpmin,'UseCpmin:            Are Cp,min data included in the airfoil tables?');
wip(fid,wt.NumAF,'NumAF:               Number of airfoil files.');
wip(fid,wt.AF_File{1},'AF_File:             List of NumAF airfoil files.');

for i=2:wt.NumAF
    wip(fid,wt.AF_File{i},[]);
end

wip(fid,[],'-----  Output Configuration  ---------------------------------------------------');
% wip(fid,[],'False                UnfPower:                  Write parametric power to an unformatted file?');
wip(fid,wt.TabDel,'TabDel:                    Make output tab-delimited (fixed-width otherwise).');
% wip(fid,[],'True                 OutNines:                  Output nines if the solution doesn''t fully converge to the specified tolerences.');
% wip(fid,[],'True                 Beep:                      Beep if errors occur.');
wip(fid,wt.KFact,'KFact:                     Output dimensional parameters in K (e.g., kN instead on N)');
wip(fid,wt.WriteBED,'WriteBED:                  Write out blade element data to "<rootname>.bed"?');
wip(fid,wt.InputTSR,'InputTSR:                  Input speeds as TSRs?');
% wip(fid,wt.OutMaxCp,'OutMaxCp:                  Output conditions for the maximum Cp?');
wip(fid,wt.SpdUnits,'SpdUnits:                  Wind-speed units (mps, fps, mph).');
wip(fid,[],'-----  Combined-Case Analysis  -------------------------------------------------');
wip(fid,wt.NumCases,'NumCases:                  Number of cases to run.  Enter zero for parametric analysis.');
wip(fid,[],'WS or TSR   RotSpd   Pitch                      Remove following block of lines if NumCases is zero.');
if wt.NumCases==0
    wip(fid,[],'-----  Parametric Analysis (Ignored if NumCases > 0 )  -------------------------');
    wip(fid,wt.par.ParRow,'ParRow:                    Row parameter    (1-rpm, 2-pitch, 3-tsr/speed).');
    wip(fid,wt.par.ParCol,'ParCol:                    Column parameter (1-rpm, 2-pitch, 3-tsr/speed).');
    wip(fid,wt.par.ParTab,'ParTab:                    Table parameter  (1-rpm, 2-pitch, 3-tsr/speed).');
    wip(fid,wt.par.OutPwr,'OutPwr:                    Request output of rotor power?');
    wip(fid,wt.par.OutCp,'OutCp:                     Request output of Cp?');
    wip(fid,wt.par.OutTrq,'OutTrq:                    Request output of shaft torque?');
    wip(fid,wt.par.OutFlp,'OutFlp:                    Request output of flap bending moment?');
    wip(fid,wt.par.OutThr,'OutThr:                    Request output of rotor thrust?');
    
    wipcsv(fid,wt.par.Pit,'PitSt, PitEnd, PitDel:     First, last, delta blade pitch (deg).');
    wipcsv(fid,wt.par.Omg,'OmgSt, OmgEnd, OmgDel:     First, last, delta rotor speed (rpm).');
    wipcsv(fid,wt.par.Spd,'SpdSt, SpdEnd, SpdDel:     First, last, delta speeds.');
end

if fid~=1, fclose(fid); end
end

%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function wip(fid,param,descrip)
% write input file parameter
if ~any(size(param)) && ~ischar(param)
    % do nothing if param = []
    % note: used ~any(size(param)) rather than isempty(param)
    %       so that unset parameters (size = [1 0]) will still
    %       get through to the following elseif statements
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

function wiplst(fid,param,descrip)
% write input file parameter list (one per line)
for k = 1:length(param)
    fprintf(fid,'%-16s ',param{k});  %output string
    if k==1
        fprintf(fid,'%s',descrip);
    end
    fprintf(fid,'\n');
end
end

function wiptbl(fid,frmt,table,nrows)
for r=1:nrows
    for c=1:length(table)
        col = table{c};
        if iscell(col)
            param = col{r};
        else
            param = col(r);
        end
        fprintf(fid,frmt{c},param);
    end
    fprintf(fid,'\n');
end
end

function wipcsv(fid,param,descrip)
% write input file parameter list as comma separated values
str = '';
for k = 1:length(param)
    str = strcat(str,sprintf('%d, ',param(k)));
end
if length(str) > 1
    str(end) = [];
end
fprintf(fid,'%-16s %s\n',str,descrip);
end