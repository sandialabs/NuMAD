########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import yaml
import numpy as np
from scipy.stats import mode

from pynumad.utils.misc_utils import LARCetaT, LARCetaL, _parse_data
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.objects.Component import Component
from pynumad.objects.Airfoil import Airfoil
from pynumad.objects.Material import Material


def yaml_to_blade(blade, filename: str, write_airfoils: bool = False):
    """
    This method writes blade information from a .yaml file to a Blade object.
    The yaml file is expected to be formatted according to the WindIO ontology.
    See https://windio.readthedocs.io/en/stable/source/turbine.html.

    Parameters
    ----------
    blade : Blade
    filename : string 
        path to .yaml file
    write_airfoils : bool
        Set true to write airfoil files while reading in data. Defaults to false.

    Returns
    -------
    blade : Blade
        input blade object populated with yaml data
    """

    # Read in yaml file as a nested dictionary
    with open(filename) as blade_yaml:
        # data = yaml.load(blade_yaml,Loader=yaml.FullLoader)
        data = yaml.load(blade_yaml,Loader=yaml.Loader)

    # Name some key subdata
    blade_outer_shape_bem = data['components']['blade']['outer_shape_bem']
    
    # older versions of wind ontology do not have 'outer_shape_bem' subsection for hub data
    try:
        hub_outer_shape_bem = data['components']['hub']['outer_shape_bem']
    except KeyError:
        hub_outer_shape_bem = data['components']['hub']
    
    blade_internal_structure = data['components']['blade']['internal_structure_2d_fem']
    af_data = data['airfoils']
    mat_data = data['materials']

    ### STATIONS / AIRFOILS
    _add_stations(blade, blade_outer_shape_bem, hub_outer_shape_bem, 
                    af_data, filename, write_airfoils)
    
    ### MATERIALS
    _add_materials(blade, mat_data)

    ## Blade Components
    N_layer_comp = len(blade_internal_structure['layers'])
    
    # Spar Cap Width and Offset
    # Obtain component name and index for hp and lp sides of sparcap
    for i in range(N_layer_comp):
        if 'spar' in blade_internal_structure['layers'][i]['name'].lower():
            name = blade_internal_structure['layers'][i]['name']
            if 'suc' in blade_internal_structure['layers'][i]['side'].lower():
                spar_lp_index = i
                spar_lp_name = name
            if 'pres' in blade_internal_structure['layers'][i]['side'].lower():
                spar_hp_index = i
                spar_hp_name = name

    # Because spar cap width must be constant, average the yaml file on
    # pressure and suction surfaces across span    
    # blade.sparcapwidth = np.zeros((2))
    blade.sparcapoffset = np.zeros((2))
    blade.sparcapwidth_hp = np.array(blade_internal_structure['layers'][spar_hp_index]['width']['values'])*1000
    blade.sparcapwidth_lp  = np.array(blade_internal_structure['layers'][spar_lp_index]['width']['values'])*1000

    blade.sparcapoffset_hp = np.array(blade_internal_structure['layers'][spar_hp_index]['offset_y_pa']['values'])*1000
    blade.sparcapoffset_lp = np.array(blade_internal_structure['layers'][spar_lp_index]['offset_y_pa']['values'])*1000
    
    # TE and LE Bands
    for i in range(N_layer_comp):
        if 'reinf' in blade_internal_structure['layers'][i]['name'].lower():
            if 'le' in blade_internal_structure['layers'][i]['name'].lower():
                I_LE = i
            else:
                I_TE = i
    
    # Leading and Trailing Edge bands are constants in millimeters

    blade.leband = np.array(blade_internal_structure['layers'][I_LE]['width']['values'])*1000 / 2
    blade.teband = np.array(blade_internal_structure['layers'][I_TE]['width']['values'])*1000 / 2
    ### COMPONENTS
    _add_components(blade, blade_internal_structure, spar_hp_index, spar_lp_index)
    
    blade.updateBlade()
    # save(blade_name)
    # BladeDef_to_NuMADfile(obj,numad_name,matdb_name,numad_af_folder)
    return blade


def _add_stations(blade,blade_outer_shape_bem, hub_outer_shape_bem,
                    af_data, file: str, write_airfoils):

    # Obtaining some parameters not explicitly given in YAML file
    L = np.ceil(blade_outer_shape_bem['reference_axis']['z']['values'][-1])
    R = L + hub_outer_shape_bem['diameter'] / 2
    L = R - hub_outer_shape_bem['diameter'] / 2
    blade.ispan = np.multiply(np.transpose(blade_outer_shape_bem['chord']['grid']),L)
    
    
    #Aerodynamic properties
    # using interp because yaml can have different r/R for twist and chord
    temp_x = np.transpose(blade_outer_shape_bem['twist']['grid'])
    temp_y = blade_outer_shape_bem['twist']['values']
    blade.degreestwist = interpolator_wrap(np.multiply(temp_x,L),np.transpose(temp_y),blade.ispan) * 180.0 / np.pi
    blade.chord = interpolator_wrap(
        np.multiply(np.transpose(blade_outer_shape_bem['chord']['grid']),L),
        np.transpose(blade_outer_shape_bem['chord']['values']),blade.ispan)
    af_dir_names = []
    for i in range(len(af_data)):
        af_dir_names.append(af_data[i]['name'])
    numstations = len(blade_outer_shape_bem['airfoil_position']['labels'])
    tc = [None]*numstations
    aero_cent = [None]*numstations

    for i in range(numstations):
        _,_,iaf_temp = np.intersect1d(blade_outer_shape_bem['airfoil_position']['labels'][i],af_dir_names,'stable',return_indices=True)
        IAF = iaf_temp[0] # Expect only one index of intersection
        tc[i] = af_data[IAF]['relative_thickness']
        tc_xL = blade_outer_shape_bem['airfoil_position']['grid'][i]
        aero_cent[i] = af_data[IAF]['aerodynamic_center']
        xf_coords = np.stack((af_data[IAF]['coordinates']['x'],
            af_data[IAF]['coordinates']['y']),1)

        # find coordinate direction (clockwise or counter-clockwise) Winding
        # Number. clockwise starting at (1,0) is correct
        with np.errstate(divide='ignore', invalid='ignore'):
            if np.nanmean(np.gradient(np.arctan(xf_coords[:,1] / xf_coords[:,0]))) > 0:
                xf_coords = np.flipud(xf_coords)

        if write_airfoils:
            import os
            out_folder = 'yaml2BladeDef_' + file.replace('.yaml','')
            # blade_name = out_folder + '/' + file.replace('.yaml','') + '_blade.mat'
            # matdb_name =...
            # numade_name =...

            # Creating folders
            os.makedirs(out_folder+'/af_coords/', exist_ok = True)
            # os.makedirs(out_folder+'/af_polars/', exist_ok = True)
            os.makedirs(out_folder+'/airfoil/', exist_ok = True)
            writeNuMADAirfoil(xf_coords, 
                blade_outer_shape_bem['airfoil_position']['labels'][i],
                out_folder + '/af_coords/' +
                blade_outer_shape_bem['airfoil_position']['labels'][i]+'.txt')

        ref = blade_outer_shape_bem['airfoil_position']['labels'][i]
        af = Airfoil(coords = xf_coords, ref = ref)
        af.resample(spacing='half-cosine')
        blade.addStation(af,tc_xL*L)
    # Obtain some key blade attributes
    blade.span = blade.ispan
    blade.percentthick = np.multiply(interpolator_wrap(np.multiply(blade_outer_shape_bem['airfoil_position']['grid'],L),tc,blade.ispan),100)
    blade.aerocenter = interpolator_wrap(np.multiply(blade_outer_shape_bem['airfoil_position']['grid'],L),aero_cent,blade.span)
    blade.chordoffset = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['pitch_axis']['grid']),L),
        np.transpose(blade_outer_shape_bem['pitch_axis']['values']),blade.span)
    blade.naturaloffset = 0
    blade.prebend = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['reference_axis']['x']['grid']),L),
        np.transpose(blade_outer_shape_bem['reference_axis']['x']['values']),blade.span)
    blade.sweep = interpolator_wrap(np.multiply(np.transpose(blade_outer_shape_bem['reference_axis']['y']['grid']),L),
        np.transpose(blade_outer_shape_bem['reference_axis']['y']['values']),blade.span)

    # for i in range(len(tc)):
    #     afc = AirfoilDef(out_folder + 
    #     '/af_coords/' + 
    #     blade_outer_shape_bem['airfoil_position']['labels'][i] +
    #     '.txt')
    #     blade.addStation(afc,np.multiply(tc_xL[i],L))

    #NOTE nothing happens to afc? Tentatively ignoring...
    # If i return to this make sure to listify the afcs
    ### AIRFOILS
    # for i in range(len(tc)):
    #     afc = AirfoilDef(out_folder + '/af_coords/' + 
    #         blade_outer_shape_bem['airfoil_position']['labels'][i] + 
    #         '.txt')
    #     blade.addStation(afc,np.multiply(tc_xL[i],L))
    # afc.resample #NOTE afc isn't used after this... why resample?
    return


def _add_materials(blade, material_data):
    materials_dict =dict()
    for i in range(len(material_data)):
        cur_mat = Material()
        cur_mat.name = material_data[i]['name']
        if material_data[i]['orth'] == 1:
            cur_mat.type = 'orthotropic'
        else:
            cur_mat.type = 'isotropic'
        # Add ply thickness option if ply thickness exists in yaml
        try:
            cur_mat.layerthickness = material_data[i]['ply_t'] * 1000
        except KeyError:
                print('Warning! material ply thickness ' + 
                        material_data[i]['name'] + 
                        ' not defined, assuming 1 mm thickness')
                cur_mat.layerthickness = 1
            
        finally:
            pass

        # first
        cur_mat.uts = _parse_data(material_data[i]['Xt'])
        cur_mat.ucs = -_parse_data(material_data[i]['Xc'])
        cur_mat.uss = _parse_data(material_data[i]['S'])
        cur_mat.xzit = 0.3
        cur_mat.xzic = 0.25
        cur_mat.yzit = 0.3
        cur_mat.yzic = 0.25
        try: 
            cur_mat.g1g2 = material_data[i].get('GIc',0) / material_data[i].get('GIIc',0)
        except ZeroDivisionError:
            cur_mat.g1g2 = np.nan
        if 'alp0' in material_data[i]:
            cur_mat.alp0 = _parse_data(material_data[i]['alp0'])
            cur_mat.etat = LARCetaT(cur_mat.alp0)
        else: 
            cur_mat.alp0 = None
            cur_mat.etat = None
        try:
            #test if property is a list
            material_data[i]['E']+[]
        except TypeError:
            cur_mat.ex = _parse_data(material_data[i]['E'])
            cur_mat.ey = _parse_data(material_data[i]['E'])
            cur_mat.ez = _parse_data(material_data[i]['E'])
            cur_mat.gxy = _parse_data(material_data[i]['G'])
            cur_mat.gxz = _parse_data(material_data[i]['G'])
            cur_mat.gyz = _parse_data(material_data[i]['G'])
            cur_mat.prxy = _parse_data(material_data[i]['nu'])
            cur_mat.prxz = _parse_data(material_data[i]['nu'])
            cur_mat.pryz = _parse_data(material_data[i]['nu'])
            cur_mat.etal = LARCetaL(cur_mat.uss,cur_mat.ucs,cur_mat.alp0)
        else:
            cur_mat.ex = _parse_data(material_data[i]['E'][0])
            cur_mat.ey = _parse_data(material_data[i]['E'][1])
            cur_mat.ez = _parse_data(material_data[i]['E'][2])
            cur_mat.gxy = _parse_data(material_data[i]['G'][0])
            cur_mat.gxz = _parse_data(material_data[i]['G'][1])
            cur_mat.gyz = _parse_data(material_data[i]['G'][2])
            cur_mat.prxy = _parse_data(material_data[i]['nu'][0])
            cur_mat.prxz = _parse_data(material_data[i]['nu'][1])
            cur_mat.pryz = _parse_data(material_data[i]['nu'][2])
            cur_mat.etal = LARCetaL(cur_mat.uss[0],cur_mat.ucs[1],cur_mat.alp0)
        try:
            cur_mat.m = material_data[i]['m']
        except KeyError:
            print(f"No fatigue exponent found for material: {material_data[i]['name']}")
        cur_mat.density = material_data[i]['rho']
        # cur_mat.dens = mat_data[i]['rho']
        cur_mat.drydensity = material_data[i]['rho']
        if 'description' in material_data[i].keys() and 'source' in material_data[i].keys():
            desc_sourc = [material_data[i]['description'],', ',material_data[i]['source']]
            cur_mat.reference = ''.join(desc_sourc)
        else:
            cur_mat.reference = []

        materials_dict[cur_mat.name] = cur_mat
    blade.materials = materials_dict
    return


def _add_components(blade, blade_internal_structure, spar_hp, spar_lp):
    N_layer_comp = len(blade_internal_structure['layers'])
    component_list = list()
    for i in range(N_layer_comp):
        i_component_data = blade_internal_structure['layers'][i]
        cur_comp = Component()
        cur_comp.group = 0
        cur_comp.name = i_component_data['name']
        #   comp['material'] = blade_internal_structure['layers']{i}['material'];
        # mat_names = [mat.name for mat in blade.materials]
        # C,IA,IB = np.intersect1d(mat_names,i_component_data['material'],return_indices=True)
        cur_comp.materialid = i_component_data['material']
        try:
            cur_comp.fabricangle = np.mean(i_component_data['fiber_orientation']['values'])
        finally:
            pass
        if 'spar' in i_component_data['name'].lower():
            cur_comp.imethod = 'pchip'
        else:
            cur_comp.imethod = 'linear'
        # cur_comp.cp[:,0] = np.transpose(i_component_data['thickness']['grid'])
        cptemp1 = np.transpose(i_component_data['thickness']['grid'])
        temp_n_layer = np.multiply(np.transpose(i_component_data['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness
        I_round_up = np.flatnonzero((temp_n_layer > 0.05) & (temp_n_layer < 0.5))
        cptemp2 = np.round(np.multiply(np.transpose(i_component_data['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness)
        cur_comp.cp = np.stack((cptemp1,cptemp2),axis=1)
        # if I_round_up.size > 0:
        #     cur_comp.cp[I_round_up,1] = 1 # increase n_layers from 0 to 1 for 0.05<n_layers<0.5
        #     comp['cp'](:,2) = cell2mat(blade_internal_structure['layers']{i}['thickness']['values'])'.*1000;  # use when each material ply is 1 mm
        cur_comp.pinnedends = 0
        component_list.append(cur_comp)



    # Spar Caps (pressure and suction)
    component_list[spar_hp].hpextents = ['b','c']
    component_list[spar_lp].lpextents = ['b','c']
    
    for comp in range(len(component_list)):

        # uv coating
        if 'uv' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['le','te']
            component_list[comp].lpextents = ['le','te']
            component_list[comp].cp[:,1] = component_list[comp].cp[:,1]

        # Shell skin1
        if 'shell_skin_outer' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['le','te']
            component_list[comp].lpextents = ['le','te']
            # CK Change me when yaml is fixed!!!!
            component_list[comp].cp[:,1] = component_list[comp].cp[:,1]
        
        # LE Band
        if 'le_reinf' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['le','a']
            component_list[comp].lpextents = ['le','a']  
        
        # TE Band
        if 'te_reinf' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['d','te']
            component_list[comp].lpextents = ['d','te']
        
        # Trailing edge suction surface panel
        if 'te_ss' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].lpextents = ['c','d']
    
        # Leading edge suction surface panel
        if 'le_ss' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].lpextents = ['a','b']

        # Leading edge pressure surface panel)
        if 'le_ps' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['a','b']

        # Trailing edge pressure surface panel
        if 'te_ps' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = ['c','d']
    
        # Shell skin2
        if 'shell_skin_inner' in blade_internal_structure['layers'][comp]['name'].lower():
            component_list[comp].hpextents = np.array(['le','te'])
            component_list[comp].lpextents = np.array(['le','te'])
            # CK Change me when yaml is fixed!!!!
            component_list[comp].cp[:,1] = component_list[comp].cp[:,1]

        # Forward Shear
        if 'web' in blade_internal_structure['layers'][comp].keys():
            if 'fore' in blade_internal_structure['layers'][comp]['web'].lower():
                # Web Skin1
                if 'skin_le' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z+' sw_offset]};
                    # comp[comp].lpextents = {['z+' sw_offset]};
                    component_list[comp].hpextents = ['b']
                    component_list[comp].lpextents = ['b']
                    component_list[comp].group = 1
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']
                    # CK Change me when yaml is fixed!!!!
                    component_list[comp].cp[:,1] = component_list[comp].cp[:,1]
                
                # Web Filler
                if 'filler' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z+' sw_offset]};
                    # comp[comp].lpextents = {['z+' sw_offset]};
                    component_list[comp].hpextents = ['b']
                    component_list[comp].lpextents = ['b']
                    component_list[comp].group = 1
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']

                # Web Skin2
                if 'skin_te' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z+' sw_offset]};
                    # comp[comp].lpextents = {['z+' sw_offset]};
                    component_list[comp].hpextents = ['b']
                    component_list[comp].lpextents = ['b']
                    component_list[comp].group = 1
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']
                    # CK Change me when yaml is fixed!!!!
                    component_list[comp].cp[:,1] = component_list[comp].cp[:,1]

        # Rear Shear
        if 'web' in blade_internal_structure['layers'][comp].keys():
            if 'rear' in blade_internal_structure['layers'][comp]['web'].lower():
                # Web Skin1
                if 'skin_le' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z-' sw_offset]};
                    # comp[comp].lpextents = {['z-' sw_offset]};
                    component_list[comp].hpextents = ['c']
                    component_list[comp].lpextents = ['c']
                    component_list[comp].group = 2
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']
                    # CK Change me when yaml is fixed!!!!
                    component_list[comp].cp[:,1] = component_list[comp].cp[:,1]
                # Web Filler
                if 'filler' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z-' sw_offset]};
                    # comp[comp].lpextents = {['z-' sw_offset]};
                    component_list[comp].hpextents = ['c']
                    component_list[comp].lpextents = ['c']
                    component_list[comp].group = 2
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']
    
                # Web Skin2
                if 'skin_te' in blade_internal_structure['layers'][comp]['name'].lower():
                    # comp[comp].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[comp].hpextents = {['z-' sw_offset]};
                    # comp[comp].lpextents = {['z-' sw_offset]};
                    component_list[comp].hpextents = ['c']
                    component_list[comp].lpextents = ['c']
                    component_list[comp].group = 2
                    component_list[comp].name = blade_internal_structure['layers'][comp]['name']
                    # CK Change me when yaml is fixed!!!!
                    component_list[comp].cp[:,1] = component_list[comp].cp[:,1]

    ### add components to blade
    component_dict = dict()
    for comp in component_list:
        component_dict[comp.name] = comp
    blade.components = component_dict
    return


def writeNuMADAirfoil(coords, reftext, fname): 
    """
    WriteNuMADAirfoil  Write NuMAD airfoil files
    **********************************************************************
    *                   Part of the SNL NuMAD Toolbox                    *
    * Developed by Sandia National Laboratories Wind Energy Technologies *
    *             See license.txt for disclaimer information             *
    **********************************************************************
      WriteNuMADAirfoil(coords,reftext,fname)
        
            fname - full filename, incl extension, of NuMAD airfoil file to write
        coords - Nx2 array of airfoil coordinate data.  First column contains
        x-values, second column contains y-values.  Airfoil coordinates are in
        order as specified by NuMAD (i.e. trailing edge = (1,0) and leading
        edge = (0,0)
        reftext = string representing reference text
    """
    with open(fname,'wt') as fid:
        fid.write('<reference>\n%s</reference>\n' % (reftext))
        fid.write('<coords>\n' % ())
        for i in range(coords.shape[0]):
            fid.write('%8.12f\t%8.12f\n' % tuple(coords[i,:]))
        
        fid.write('</coords>' % ())