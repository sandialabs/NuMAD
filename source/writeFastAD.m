function writeFastAD(ad,output_file)
%WRITEFASTAD  Write a FAST Aerodyn input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeFastAD(ad,file_name) 
%      For use with AeroDyn v13.00.00
% 
%      ad = AeroDyn data structure; view default structure by utilizing 
%           readFastAD()
%      file_name = name of file to be created.


if ~exist('output_file','var')  %print to command window if no output_file given
    fid=1;
else
    fid=fopen(output_file,'wt');   %try to open output_file for Writing in Text mode
    if (fid == -1)
        error('Could not open file "%s"\n',output_file);
        return
    end
end

wip(fid,[],ad.title{1});
wip(fid,ad.SysUnits,      'SysUnits - System of units for used for input and output [must be SI for FAST] (unquoted string)')
wip(fid,ad.StallMod,      'StallMod - Dynamic stall included [BEDDOES or STEADY] (unquoted string)')
wip(fid,ad.UseCm,         'UseCm    - Use aerodynamic pitching moment model? [USE_CM or NO_CM] (unquoted string)')
wip(fid,ad.InfModel,      'InfModel - Inflow model [DYNIN or EQUIL] (unquoted string)')
wip(fid,ad.IndModel,      'IndModel - Induction-factor model [NONE or WAKE or SWIRL] (unquoted string)')
wip(fid,ad.AToler,        'AToler   - Induction-factor tolerance (convergence criteria) (-)')
wip(fid,ad.TLModel,       'TLModel  - Tip-loss model (EQUIL only) [PRANDtl, GTECH, or NONE] (unquoted string)')
wip(fid,ad.HLModel,       'HLModel  - Hub-loss model (EQUIL only) [PRANdtl or NONE] (unquoted string)')
wip(fid,ad.WindFile,      'WindFile - Name of file containing wind data (quoted string)')
wip(fid,ad.HH,            'HH       - Wind reference (hub) height [TowerHt+Twr2Shft+OverHang*SIN(ShftTilt)] (m)')
wip(fid,ad.TwrShad,       'TwrShad  - Tower-shadow velocity deficit (-)')
wip(fid,ad.ShadHWid,      'ShadHWid - Tower-shadow half width (m)')
wip(fid,ad.T_Shad_Refpt,  'T_Shad_Refpt - Tower-shadow reference point (m)')
wip(fid,ad.Rho,           'Rho      - Air density (kg/m^3)')
wip(fid,ad.KinVisc,       'KinVisc  - Kinematic air viscosity [CURRENTLY IGNORED] (m^2/sec)')
wip(fid,ad.DTAero,        'DTAero   - Time interval for aerodynamic calculations (sec)')
wip(fid,ad.NumFoil,       'NumFoil  - Number of airfoil files (-)')
wiplst(fid,ad.FoilNm,     'FoilNm   - Names of the airfoil files [NumFoil lines] (quoted strings)')
wip(fid,ad.BldNodes,      'BldNodes - Number of blade nodes used for analysis (-)')
wip(fid,[],               'RNodes    AeroTwst  DRNodes  Chord  NFoil  PrnElm')
ElmTable = {ad.RNodes, ad.AeroTwst, ad.DRNodes, ad.Chord, ad.NFoil, ad.PrnElm};
wiptbl(fid,{'%8.5f','%7.2f','%12.5f','%7.3f','%4d','     %s'},ElmTable,ad.BldNodes);
wip(fid,ad.MultTab,       '');

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