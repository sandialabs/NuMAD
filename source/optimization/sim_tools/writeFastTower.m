function writeFastTower(twr,output_file)
%WRITEFASTTOWER  Write a FAST tower input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeFastTower(twr,'file_name') 
%      Created for Fast v7.00
%
%      twr = FAST Tower data structure; view default structure 
%          by utilizing readFastTower()
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
wip(fid,[],'---------------------- FAST TOWER FILE -----------------------------------------')
wip(fid,[],twr.title{1})
wip(fid,[],'---------------------- TOWER PARAMETERS ----------------------------------------')
wip(fid,twr.NTwInpSt,    'NTwInpSt    - Number of input stations to specify tower geometry')
wip(fid,twr.CalcTMode,   'CalcTMode   - Calculate tower mode shapes internally {T: ignore mode shapes from below, F: use mode shapes from below} [CURRENTLY IGNORED] (flag)')
wip(fid,twr.TwrFADmp(1), 'TwrFADmp(1) - Tower 1st fore-aft mode structural damping ratio (%%)')
wip(fid,twr.TwrFADmp(2), 'TwrFADmp(2) - Tower 2nd fore-aft mode structural damping ratio (%%)')
wip(fid,twr.TwrSSDmp(1), 'TwrSSDmp(1) - Tower 1st side-to-side mode structural damping ratio (%%)')
wip(fid,twr.TwrSSDmp(2), 'TwrSSDmp(2) - Tower 2nd side-to-side mode structural damping ratio (%%)')
wip(fid,[],'---------------------- TOWER ADJUSTMUNT FACTORS --------------------------------')
wip(fid,twr.FAStTunr(1), 'FAStTunr(1) - Tower fore-aft modal stiffness tuner, 1st mode (-)')
wip(fid,twr.FAStTunr(2), 'FAStTunr(2) - Tower fore-aft modal stiffness tuner, 2nd mode (-)')
wip(fid,twr.SSStTunr(1), 'SSStTunr(1) - Tower side-to-side stiffness tuner, 1st mode (-)')
wip(fid,twr.SSStTunr(2), 'SSStTunr(2) - Tower side-to-side stiffness tuner, 2nd mode (-)')
wip(fid,twr.AdjTwMa,     'AdjTwMa     - Factor to adjust tower mass density (-)')
wip(fid,twr.AdjFASt,     'AdjFASt     - Factor to adjust tower fore-aft stiffness (-)')
wip(fid,twr.AdjSSSt,     'AdjSSSt     - Factor to adjust tower side-to-side stiffness (-)')
wip(fid,[],'---------------------- DISTRIBUTED TOWER PROPERTIES ----------------------------')
wip(fid,[],'HtFract  TMassDen  TwFAStif      TwSSStif      TwGJStif      TwEAStif      TwFAIner  TwSSIner  TwFAcgOf  TwSScgOf')
wip(fid,[],'(-)      (kg/m)    (Nm^2)        (Nm^2)        (Nm^2)        (N)           (kg m)    (kg m)    (m)       (m)')
PropTable = {twr.prop.HtFract, twr.prop.TMassDen, twr.prop.TwFAStif, ...
             twr.prop.TwSSStif, twr.prop.TwGJStif, twr.prop.TwEAStif, ...
             twr.prop.TwFAIner, twr.prop.TwSSIner, twr.prop.TwFAcgOf, ...
             twr.prop.TwSScgOf};
Frmt = {'%7.5f  ',  '%8.3f  ',  '%12.5e  ', ...
        '%12.5e  ', '%12.5e  ', '%12.5e  ', ...
        '%7.2f   ',  '%7.2f  ',  '%6.3f    ', ...
        '%6.3f'};
wiptbl(fid,Frmt,PropTable,twr.NTwInpSt);
wip(fid,[],'---------------------- TOWER FORE-AFT MODE SHAPES ------------------------------')
wip(fid,twr.TwFAM1Sh(2), 'TwFAM1Sh(2) - Mode 1, coefficient of x^2 term')
wip(fid,twr.TwFAM1Sh(3), 'TwFAM1Sh(3) -       , coefficient of x^3 term')
wip(fid,twr.TwFAM1Sh(4), 'TwFAM1Sh(4) -       , coefficient of x^4 term')
wip(fid,twr.TwFAM1Sh(5), 'TwFAM1Sh(5) -       , coefficient of x^5 term')
wip(fid,twr.TwFAM1Sh(6), 'TwFAM1Sh(6) -       , coefficient of x^6 term')
wip(fid,twr.TwFAM2Sh(2), 'TwFAM2Sh(2) - Mode 2, coefficient of x^2 term')
wip(fid,twr.TwFAM2Sh(3), 'TwFAM2Sh(3) -       , coefficient of x^3 term')
wip(fid,twr.TwFAM2Sh(4), 'TwFAM2Sh(4) -       , coefficient of x^4 term')
wip(fid,twr.TwFAM2Sh(5), 'TwFAM2Sh(5) -       , coefficient of x^5 term')
wip(fid,twr.TwFAM2Sh(6), 'TwFAM2Sh(6) -       , coefficient of x^6 term')
wip(fid,[],'---------------------- TOWER SIDE-TO-SIDE MODE SHAPES --------------------------')
wip(fid,twr.TwSSM1Sh(2), 'TwSSM1Sh(2) - Mode 1, coefficient of x^2 term')
wip(fid,twr.TwSSM1Sh(3), 'TwSSM1Sh(3) -       , coefficient of x^3 term')
wip(fid,twr.TwSSM1Sh(4), 'TwSSM1Sh(4) -       , coefficient of x^4 term')
wip(fid,twr.TwSSM1Sh(5), 'TwSSM1Sh(5) -       , coefficient of x^5 term')
wip(fid,twr.TwSSM1Sh(6), 'TwSSM1Sh(6) -       , coefficient of x^6 term')
wip(fid,twr.TwSSM2Sh(2), 'TwSSM2Sh(2) - Mode 2, coefficient of x^2 term')
wip(fid,twr.TwSSM2Sh(3), 'TwSSM2Sh(3) -       , coefficient of x^3 term')
wip(fid,twr.TwSSM2Sh(4), 'TwSSM2Sh(4) -       , coefficient of x^4 term')
wip(fid,twr.TwSSM2Sh(5), 'TwSSM2Sh(5) -       , coefficient of x^5 term')
wip(fid,twr.TwSSM2Sh(6), 'TwSSM2Sh(6) -       , coefficient of x^6 term')

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