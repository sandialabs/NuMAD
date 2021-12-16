function varargout = FASTBlade2BModes(bld,blLen,hubRad,rpmOper,bmodes_path)

% bld = data structure as read in by readFASTBlade()
% blLen [m]
% hubRad [m]
% rpmOper [rpm]

BlLength=blLen;

bmodes.title{1}='BModes analysis to based on FAST Blade file data';
bmodes.title{2}='';
bmodes.hub_rad=hubRad;  % hub radius measured along coned blade axis OR tower rigid-base height (m)
bmodes.radius=bmodes.hub_rad+BlLength; % rotor tip radius measured along coned blade axis OR tower height (m)
bmodes.precone=0; % built-in precone angle, automatically set to zero for a tower (deg)
bmodes.bl_thp=0; % blade pitch setting, automatically set to zero for a tower (deg)
bmodes.modepr=6; % number of modes to be printed (-)
bmodes.el_loc=[0.  0.08  0.16  0.24  0.32  0.40  0.48  0.56  0.64  0.72  0.80  0.90   1.0];
bmodes.el_loc=linspace(0,1,30);

bmodes.Echo='T';
bmodes.beam_type=1;
bmodes.rpm_mult=1;
bmodes.hub_conn=1;
bmodes.TabDelim='T';
bmodes.mid_node_tw='F';
bmodes.tip_mass=0;
bmodes.cm_loc=0;
bmodes.cm_axial=0;
bmodes.ixx_tip=0;
bmodes.iyy_tip=0;
bmodes.izz_tip=0;
bmodes.ixy_tip=0;
bmodes.izx_tip=0;
bmodes.iyz_tip=0;
bmodes.id_mat=1;
bmodes.sec_props_file='"blade_sec_props.dat"';
bmodes.sec_mass_mult=1;
bmodes.flp_iner_mult=1;
bmodes.lag_iner_mult=1;
bmodes.flp_stff_mult=1;
bmodes.edge_stff_mult=1;
bmodes.tor_stff_mult=1;
bmodes.axial_stff_mult=1;
bmodes.cg_offst_mult=1;
bmodes.sc_offst_mult=1;
bmodes.tc_offst_mult=1;

n=length(bld.prop.BlFract);

bmodes.props.sec_loc=bld.prop.BlFract';
bmodes.props.str_tw=bld.prop.StrcTwst';
bmodes.props.tw_iner=ones(n,1)*1e12;
bmodes.props.mass_den=bld.prop.BMassDen';
bmodes.props.flp_iner=bld.prop.FlpIner';
bmodes.props.edge_iner=bld.prop.EdgIner';
bmodes.props.flp_stff=bld.prop.FlpStff';
bmodes.props.edge_stff=bld.prop.EdgStff';
bmodes.props.tor_stff=bld.prop.GJStff';
bmodes.props.axial_stff=bld.prop.EAStff';
bmodes.props.cg_offst=bld.prop.FlpcgOf';
bmodes.props.sc_offst=bld.prop.FlpEAOf';
bmodes.props.tc_offst=bld.prop.FlpEAOf';  % assuming that tension center and shear center are co-located in this case

bmodes.rot_rpm=rpmOper; % rotor speed, automatically set to zero for tower modal analysis (rpm)
writeBModesMain(bmodes,'bmodes.bmi');
writeBModesSecProps(bmodes,'blade_sec_props.dat');

% Run BModes
[~, w] = dos(sprintf('"%s" bmodes.bmi',bmodes_path));
disp(w)

end
