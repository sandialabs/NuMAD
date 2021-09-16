function writeBModesMain(bmodes,output_file)
%WRITEBMODESMAIN  Write a BModes primary input file.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   writeBModesMain(bmodes,file_name)
% 
%      bmodes = BModes input structure; view default structure by 
%               examining this script
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

wip(fid,[],'======================   BModes v3.00 Main Input File  ==================')
wip(fid,[],bmodes.title{1});
wip(fid,[],bmodes.title{2});
wip(fid,[],'--------- General parameters ---------------------------------------------------------------------');
wip(fid,bmodes.Echo,'Echo        Echo input file contents to *.echo file if true.');
wip(fid,bmodes.beam_type,'beam_type   1: blade, 2: tower (-)');
wip(fid,bmodes.rot_rpm,'rot_rpm:    rotor speed, automatically set to zero for tower modal analysis (rpm)');
wip(fid,bmodes.rpm_mult,'rpm_mult:   rotor speed muliplicative factor (-)');
wip(fid,bmodes.radius,'radius:     rotor tip radius measured along coned blade axis OR tower height (m)');
wip(fid,bmodes.hub_rad,'hub_rad:    hub radius measured along coned blade axis OR tower rigid-base height (m)');
wip(fid,bmodes.precone,'precone:    built-in precone angle, automatically set to zero for a tower (deg)');
wip(fid,bmodes.bl_thp,'bl_thp:     blade pitch setting, automatically set to zero for a tower (deg)');
wip(fid,bmodes.hub_conn,'hub_conn:   hub-to-blade or tower-base boundary condition [1: cantilevered; 2: free-free; 3: only axial and torsion constraints] (-)');
wip(fid,bmodes.modepr,'modepr:     number of modes to be printed (-)');
wip(fid,bmodes.TabDelim,'TabDelim    (true: tab-delimited output tables; false: space-delimited tables)');
wip(fid,bmodes.mid_node_tw,'mid_node_tw  (true: output twist at mid-node of elements; false: no mid-node outputs)');
wip(fid,[],'');
wip(fid,[],'--------- Blade-tip or tower-top mass properties --------------------------------------------');
wip(fid,bmodes.tip_mass,'tip_mass    blade-tip or tower-top mass (kg)');
wip(fid,bmodes.cm_loc,'cm_loc      tip-mass c.m. offset from the blade axis measured along the tip section y reference axis (m)');
wip(fid,bmodes.cm_axial,'cm_axial    tip-mass c.m. offset tower tip measures axially along the z axis (m)');
wip(fid,bmodes.ixx_tip,'ixx_tip     blade lag mass moment of inertia about the tip-section x reference axis (kg-m^2)');
wip(fid,bmodes.iyy_tip,'iyy_tip     blade flap mass moment of inertia about the tip-section y reference axis (kg-m^2)');
wip(fid,bmodes.izz_tip,'izz_tip     torsion mass moment of inertia about the tip-section z reference axis (kg-m^2)');
wip(fid,bmodes.ixy_tip,'ixy_tip     cross product of inertia about x and y reference axes(kg-m^2)');
wip(fid,bmodes.izx_tip,'izx_tip     cross product of inertia about z and x reference axes(kg-m^2)');
wip(fid,bmodes.iyz_tip,'iyz_tip     cross product of inertia about y and z reference axes(kg-m^2)');
wip(fid,[],'');
wip(fid,[],'--------- Distributed-property identifiers --------------------------------------------------------');
wip(fid,bmodes.id_mat,'id_mat:     material_type [1: isotropic; non-isotropic composites option not yet available]');
wip(fid,bmodes.sec_props_file,'sec_props_file   name of beam section properties file (-)');
wip(fid,[],'');
wip(fid,[],'Property scaling factors..............................');
wip(fid,bmodes.sec_mass_mult,'sec_mass_mult:   mass density multiplier (-)');
wip(fid,bmodes.flp_iner_mult,'flp_iner_mult:   blade flap or tower f-a inertia multiplier (-)');
wip(fid,bmodes.lag_iner_mult,'lag_iner_mult:   blade lag or tower s-s inertia multiplier (-)');
wip(fid,bmodes.flp_stff_mult,'flp_stff_mult:   blade flap or tower f-a bending stiffness multiplier (-)');
wip(fid,bmodes.edge_stff_mult,'edge_stff_mult:  blade lag or tower s-s bending stiffness multiplier (-)');
wip(fid,bmodes.tor_stff_mult,'tor_stff_mult:   torsion stiffness multiplier (-)');
wip(fid,bmodes.axial_stff_mult,'axial_stff_mult: axial stiffness multiplier (-)');
wip(fid,bmodes.cg_offst_mult,'cg_offst_mult:   cg offset multiplier (-)');
wip(fid,bmodes.sc_offst_mult,'sc_offst_mult:   shear center multiplier (-)');
wip(fid,bmodes.tc_offst_mult,'tc_offst_mult:   tension center multiplier (-)');
wip(fid,[],'');
wip(fid,[],'--------- Finite element discretization --------------------------------------------------');
wip(fid,length(bmodes.el_loc)-1,'nselt:     no of blade or tower elements (-)');
wip(fid,[],'Distance of element boundary nodes from blade or flexible-tower root (normalized wrt blade or tower length), el_loc()');
wiphlst(fid,'%5.3f  ',bmodes.el_loc);
wip(fid,[],'');
wip(fid,[],'END of Main Input File Data *********************************************************************');
wip(fid,[],'*************************************************************************************************');

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

function wiphlst(fid,Frmt,param)
% write input file horizontal list of numbers (all on one line, separated by spaces)
for k = 1:length(param)
    fprintf(fid,Frmt,param(k));
end
fprintf(fid,'\n');
end
