function writeBModesSecProps(bmodes,output_file)
% WRITEBMODESSECPROP  Write a BModes section properties input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeBModesSecProps(bmodes,output_file)
%      bmodes = BModes data structure; view default structure by examining
%               this script
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

n_secs=length(bmodes.props.sec_loc);
wip(fid,[],'Blade section properties')
wip(fid,n_secs,'n_secs:     number of blade sections at which properties are specified (-)')
wip(fid,[],' ')
wip(fid,[],'sec_loc     str_tw     tw_iner     mass_den   flp_iner    edge_iner    flp_stff      edge_stff     tor_stff    axial_stff   cg_offst   sc_offst  tc_offst')
wip(fid,[],'  (-)       (deg)       (deg)       (kg/m)    (kg-m)      (kg-m)        (Nm^2)        (Nm^2)        (Nm^2)        (N)         (m)        (m)       (m)')

PropTable = {bmodes.props.sec_loc, bmodes.props.str_tw, bmodes.props.tw_iner,...
    bmodes.props.mass_den, bmodes.props.flp_iner, bmodes.props.edge_iner, ...
    bmodes.props.flp_stff, bmodes.props.edge_stff, bmodes.props.tor_stff, ...
    bmodes.props.axial_stff, bmodes.props.cg_offst, bmodes.props.sc_offst, ...
    bmodes.props.tc_offst};

Frmt = {'  %7.5f\t', '% 7.5f\t', '% 7.5f\t', ...
    '% 7.2f\t', ' % 7.5f \t', ' % 7.5f \t', ...
    ' % 5.3e\t', '% 5.3e\t', ' % 5.3e\t', ...
    ' % 5.3e\t', '% 7.5f\t', '% 7.5f\t', ...
    '% 7.5f'};

wiptbl(fid,Frmt,PropTable,n_secs);

wip(fid,[],' ')
wip(fid,[],' ')
wip(fid,[],'**Note: If the above data represents TOWER properties, the following are overwritten:')
wip(fid,[],'  str_tw is set to zero')
wip(fid,[],'  tw_iner is set to zero')
wip(fid,[],'  cg_offst is set to zero')
wip(fid,[],'  sc_offst is set to zero')
wip(fid,[],'  tc_offst is set to zero')
wip(fid,[],'  edge_iner is set equal to flp_iner')
wip(fid,[],'  edge_stff is set equal to flp_stff')

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