function writeFastBlade(blade,output_file)
%WRITEFASTBLADE  Write a FAST blade input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeFastBlade(blade,'file_name') 
%      Created for Fast v7.00
%
%      blade = FAST blade data structure; view default structure 
%          by utilizing readFastBlade()
%      output_file = name of file to be created.

if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],'--------------------------------------------------------------------------------')
wip(fid,[],'---------------------- FAST INDIVIDUAL BLADE FILE ------------------------------')
wip(fid,[],blade.title{1});
wip(fid,[],'---------------------- BLADE PARAMETERS ----------------------------------------')
wip(fid,blade.NBlInpSt,    'NBlInpSt    - Number of blade input stations (-)')
wip(fid,blade.CalcBMode,   'CalcBMode   - Calculate blade mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)')
wip(fid,blade.BldFlDmp(1), 'BldFlDmp(1) - Blade flap mode #1 structural damping in percent of critical (%%)')
wip(fid,blade.BldFlDmp(2), 'BldFlDmp(2) - Blade flap mode #2 structural damping in percent of critical (%%)')
wip(fid,blade.BldEdDmp(1), 'BldEdDmp(1) - Blade edge mode #1 structural damping in percent of critical (%%)')
wip(fid,[],'---------------------- BLADE ADJUSTMENT FACTORS --------------------------------')
wip(fid,blade.FlStTunr(1), 'FlStTunr(1) - Blade flapwise modal stiffness tuner, 1st mode (-)')
wip(fid,blade.FlStTunr(2), 'FlStTunr(2) - Blade flapwise modal stiffness tuner, 2nd mode (-)')
wip(fid,blade.AdjBlMs,     'AdjBlMs     - Factor to adjust blade mass density (-)')
wip(fid,blade.AdjFlSt,     'AdjFlSt     - Factor to adjust blade flap stiffness (-)')
wip(fid,blade.AdjEdSt,     'AdjEdSt     - Factor to adjust blade edge stiffness (-)')
wip(fid,[],'---------------------- DISTRIBUTED BLADE PROPERTIES ----------------------------')
wip(fid,[],'BlFract  AeroCent  StrcTwst  BMassDen  FlpStff       EdgStff       GJStff        EAStff        Alpha  FlpIner  EdgIner  PrecrvRef  PreswpRef  FlpcgOf  EdgcgOf  FlpEAOf  EdgEAOf')
wip(fid,[],'(-)      (-)       (deg)     (kg/m)    (Nm^2)        (Nm^2)        (Nm^2)        (N)           (-)    (kg m)   (kg m)   (m)        (m)        (m)      (m)      (m)      (m)')
PropTable = {blade.prop.BlFract, blade.prop.AeroCent, blade.prop.StrcTwst, ...
             blade.prop.BMassDen, blade.prop.FlpStff, blade.prop.EdgStff, ...
             blade.prop.GJStff, blade.prop.EAStff, blade.prop.Alpha, ...
             blade.prop.FlpIner, blade.prop.EdgIner, blade.prop.PrecrvRef, ...
             blade.prop.PreswpRef, blade.prop.FlpcgOf, blade.prop.EdgcgOf, ...
             blade.prop.FlpEAOf, blade.prop.EdgEAOf};
Frmt = {'%7.5f  ', '%5.3f     ', '%6.3f    ', ...
        '%8.3f  ', '%12.5e  ', '%12.5e  ', ...
        '%12.5e  ', '%12.5e  ', '%5.3f  ', ...
        '%7.3f  ', '%7.3f ', '%6.3f     ',...
        '%6.3f     ', '%6.3f   ', '%6.3f   ',...
        '%6.3f   ', '%6.3f'};
wiptbl(fid,Frmt,PropTable,blade.NBlInpSt);
wip(fid,[],'---------------------- BLADE MODE SHAPES ---------------------------------------')
wip(fid,blade.BldFl1Sh(2), 'BldFl1Sh(2) - Flap mode 1, coeff of x^2')
wip(fid,blade.BldFl1Sh(3), 'BldFl1Sh(3) -            , coeff of x^3')
wip(fid,blade.BldFl1Sh(4), 'BldFl1Sh(4) -            , coeff of x^4')
wip(fid,blade.BldFl1Sh(5), 'BldFl1Sh(5) -            , coeff of x^5')
wip(fid,blade.BldFl1Sh(6), 'BldFl1Sh(6) -            , coeff of x^6')
wip(fid,blade.BldFl2Sh(2), 'BldFl2Sh(2) - Flap mode 2, coeff of x^2')
wip(fid,blade.BldFl2Sh(3), 'BldFl2Sh(3) -            , coeff of x^3')
wip(fid,blade.BldFl2Sh(4), 'BldFl2Sh(4) -            , coeff of x^4')
wip(fid,blade.BldFl2Sh(5), 'BldFl2Sh(5) -            , coeff of x^5')
wip(fid,blade.BldFl2Sh(6), 'BldFl2Sh(6) -            , coeff of x^6')
wip(fid,blade.BldEdgSh(2), 'BldEdgSh(2) - Edge mode 1, coeff of x^2')
wip(fid,blade.BldEdgSh(3), 'BldEdgSh(3) -            , coeff of x^3')
wip(fid,blade.BldEdgSh(4), 'BldEdgSh(4) -            , coeff of x^4')
wip(fid,blade.BldEdgSh(5), 'BldEdgSh(5) -            , coeff of x^5')
wip(fid,blade.BldEdgSh(6), 'BldEdgSh(6) -            , coeff of x^6')

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