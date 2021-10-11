function varargout = BModes4Blade2PreComp2FASTBlade(blade,bmodes_path,PreComp_SectionData)

BlLength=blade.ispan(end);

bmodes.title{1}='BModes analysis to support NuMAD to FAST Blade using PreComp';
bmodes.title{2}='';
bmodes.hub_rad=1.5;  % hub radius measured along coned blade axis OR tower rigid-base height (m)
bmodes.radius=bmodes.hub_rad+BlLength; % rotor tip radius measured along coned blade axis OR tower height (m)
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

bmodes.props.sec_loc=PreComp_SectionData(:,1)';
bmodes.props.str_tw=PreComp_SectionData(:,3)';
bmodes.props.tw_iner=PreComp_SectionData(:,21)';
bmodes.props.mass_den=PreComp_SectionData(:,18)';
bmodes.props.flp_iner=PreComp_SectionData(:,19)';
bmodes.props.edge_iner=PreComp_SectionData(:,20)';
bmodes.props.flp_stff=PreComp_SectionData(:,4)';
bmodes.props.edge_stff=PreComp_SectionData(:,5)';
bmodes.props.tor_stff=PreComp_SectionData(:,6)';
bmodes.props.axial_stff=PreComp_SectionData(:,7)';
bmodes.props.cg_offst=PreComp_SectionData(:,23)';
bmodes.props.sc_offst=PreComp_SectionData(:,15)';
bmodes.props.tc_offst=PreComp_SectionData(:,17)';

%for rpm=0:5:60
for rpm=20
    bmodes.rot_rpm=rpm; % rotor speed, automatically set to zero for tower modal analysis (rpm)
    writeBModesMain(bmodes,sprintf('bmodes_%drpm.bmi',rpm));
    writeBModesSecProps(bmodes,'blade_sec_props.dat');
    % Run BModes
    [~, w] = dos(sprintf('"%s" bmodes_%drpm.bmi',bmodes_path,rpm));
    disp(w)
    if (nargout > 0 && rpm==20)
        % use regular expression to parse 'mode   1  frequency (hz) =    0.615'
        t = regexp(w,'mode\s*(?<mode>\d+)\s*frequency \(hz\) =\s*(?<freq>[0-9\-\.]+)','names');
        modenumber = str2double({t.mode});
        if modenumber(:) ~= transpose(1:numel(modenumber))
            disp('BModes4NuMAD2PreComp2FASTBlade: Warning - BModes results not in expected order.')
        end
        freq = str2double({t.freq});
        varargout{1} = freq(:);
    end
end

end
