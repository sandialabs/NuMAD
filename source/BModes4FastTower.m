function BModes4FastTower(twr,h,bmodes_path)
% BMODES4FASTTOWER description
%                         Under construction
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%

bmodes.title{1}='BModes analysis to provide mode shapes for FAST Tower';
bmodes.title{2}='';
bmodes.rot_rpm=0; % rotor speed, automatically set to zero for tower modal analysis (rpm)
bmodes.hub_rad=0;  % hub radius measured along coned blade axis OR tower rigid-base height (m)
bmodes.radius=h; % rotor tip radius measured along coned blade axis OR tower height (m)
bmodes.precone=0; % built-in precone angle, automatically set to zero for a tower (deg)
bmodes.bl_thp=0; % blade pitch setting, automatically set to zero for a tower (deg)
bmodes.modepr=6; % number of modes to be printed (-)
bmodes.el_loc=[0.  0.08  0.16  0.24  0.32  0.40  0.48  0.56  0.64  0.72  0.80  0.90   1.0];

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
bmodes.sec_props_file='"tower_sec_props.dat"';
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

z=twr.prop.HtFract*0; %zeros

bmodes.props.sec_loc=twr.prop.HtFract;
bmodes.props.str_tw=z;
bmodes.props.tw_iner=z;
bmodes.props.mass_den=twr.prop.TMassDen;
bmodes.props.flp_iner=twr.prop.TwFAIner;
bmodes.props.edge_iner=bmodes.props.flp_iner;
bmodes.props.flp_stff=twr.prop.TwFAStif;
bmodes.props.edge_stff=bmodes.props.flp_stff;
bmodes.props.tor_stff=twr.prop.TwGJStif;
bmodes.props.axial_stff=twr.prop.TwEAStif;
bmodes.props.cg_offst=z;
bmodes.props.sc_offst=z;
bmodes.props.tc_offst=z;

writeBModesMain(bmodes,'bmodes.bmi');
writeBModesSecProps(bmodes,'tower_sec_props.dat');

% Run BModes
dos([bmodes_path ' bmodes.bmi']);

end
