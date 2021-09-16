function BModes4NuMAD2BPE2FASTBlade(bmodes_path)
%BModes4NuMAD2BPE2FASTBlade  Generate BModes analysis to support BPE2FASTBlade
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   BModes4NuMAD2BPE2FASTBlade(bmodes_path)
%
%   Generate BModes analysis to support BPE2FASTBlade
%
%   Requires availability of post-processed BPE results.  Use BPEPost.m to
%   generate "BPE_SectionData.mat".
%

%===== CREDITS & CHANGELOG ================================================
% yyyy.mm.dd  initials: description
% yyyy.mm.dd  initials: description

load BPE_SectionData
bmodes_input='bmodes.bmi';

bmodes.title{1}='Bmodes title 1';
bmodes.title{2}='Bmodes title 2';
bmodes.Echo='True';
bmodes.beam_type=1;
bmodes.rot_rpm=15;
bmodes.rpm_mult=1;
bmodes.radius=zbeamnode_mid(end);
bmodes.hub_rad=0;
bmodes.precone=0;
bmodes.bl_thp=0;
bmodes.hub_conn=1;
bmodes.modepr=6;
bmodes.TabDelim='True';
bmodes.mid_node_tw='False';
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
bmodes.sec_props_file='bmodes_sec_props.dat';
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
% bmodes.el_loc=[0.  0.08  0.16  0.24  0.32  0.40  0.48  0.56  0.64  0.72  0.80  0.90  1.0];
bmodes.el_loc=linspace(0,1,11);

writeBModesMain(bmodes,bmodes_input);

bmodes.props.sec_loc=zbeamnode_mid/zbeamnode_mid(end);
bmodes.props.str_tw=rot_aero*-57.3*rotdir;
bmodes.props.tw_iner=rot_aero*-57.3*rotdir;
bmodes.props.mass_den=massperlen;
bmodes.props.flp_iner=Iyyperlen_cg;
bmodes.props.edge_iner=Ixxperlen_cg;
bmodes.props.flp_stff=EI_flap;
bmodes.props.edge_stff=EI_edge;
bmodes.props.tor_stff=GJ;
bmodes.props.axial_stff=EA;
bmodes.props.cg_offst=ycg;
% The following is more correct:
% bmodes.props.sc_offst=y_sh_off_sect;
% bmodes.props.tc_offst=y_el_off_sect;
% but the following is more robust, and is close to correct for use in getting basic mode shapes:
bmodes.props.sc_offst=zeros(length(zbeamnode_mid),1);
bmodes.props.tc_offst=zeros(length(zbeamnode_mid),1);

writeBModesSecProps(bmodes,'bmodes_sec_props.dat');

dos(sprintf('"%s" %s',bmodes_path,bmodes_input));

end



