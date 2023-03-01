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


def yaml_to_blade(blade, filename: str):
    """
    TODO docstring
    """
    #NOTE this ignores warnings from numpy - comment for debugging
    # np.seterr(invalid='ignore')

    # Read in NREL yaml file
    with open(filename) as blade_yaml:
        data = yaml.load(blade_yaml,Loader=yaml.FullLoader)

    # Name some key subdata
    blade_outer_shape_bem = data['components']['blade']['outer_shape_bem']
    hub_outer_shape_bem = data['components']['hub']['outer_shape_bem']
    blade_internal_structure = data['components']['blade']['internal_structure_2d_fem']
    af_data = data['airfoils']
    mat_data = data['materials']

    ### STATIONS / AIRFOILS
    add_stations(blade, blade_outer_shape_bem, hub_outer_shape_bem, 
                    af_data, filename)
    
    ### MATERIALS
    add_materials(blade, mat_data)

    ## Blade Components
    N_layer_comp = len(blade_internal_structure['layers'])
    
    # Spar Cap Width and Offset
    I_spar_lp = []
    I_spar_hp = []
    for i in range(N_layer_comp):
        if 'spar' in blade_internal_structure['layers'][i]['name'].lower():
            if 'suc' in blade_internal_structure['layers'][i]['side'].lower():
                I_spar_lp.append(i)
            if 'pres' in blade_internal_structure['layers'][i]['side'].lower():
                I_spar_hp.append(i)

    # Syntax from NuMAD 3.0 suggests that I_spar_hp/lp should always be
    # length 1. -kb
    I_spar_hp = I_spar_hp[0]
    I_spar_lp = I_spar_lp[0]
    # Because spar cap width must be constant, average the yaml file on
    # pressure and suction surfaces across span
    blade.sparcapwidth = np.zeros((2))
    blade.sparcapoffset = np.zeros((2))
    blade.sparcapwidth[0] = np.multiply(mode(blade_internal_structure['layers'][I_spar_hp]['width']['values'], keepdims = True).mode[0],1000)
    blade.sparcapwidth[1] = np.multiply(mode(blade_internal_structure['layers'][I_spar_lp]['width']['values'], keepdims = True).mode[0],1000)
    blade.sparcapoffset[0] = np.multiply(np.mean(blade_internal_structure['layers'][I_spar_hp]['offset_y_pa']['values']),1000)
    blade.sparcapoffset[1] = np.multiply(np.mean(blade_internal_structure['layers'][I_spar_lp]['offset_y_pa']['values']),1000)
    
    # TE and LE Bands
    I_reinf = []
    I_LE = []
    I_TE = []
    for i in range(N_layer_comp):
        if 'reinf' in blade_internal_structure['layers'][i]['name'].lower():
            I_reinf.append(i) # I_reinf not used for anything currently
            if 'le' in blade_internal_structure['layers'][i]['name'].lower():
                I_LE.append(i)
            else:
                I_TE.append(i)
    
    # Syntax from NuMAD 3.0 suggests that I_LE/TE should always be
    # length 1. -kb
    I_LE = I_LE[0]
    I_TE = I_TE[0]
    # Leading and Trailing Edge bands are constants in millimeters
    blade.leband = np.multiply(np.mean(blade_internal_structure['layers'][I_LE]['width']['values']),1000) / 2
    blade.teband = np.multiply(np.mean(blade_internal_structure['layers'][I_TE]['width']['values']),1000) / 2
    
    ### COMPONENTS
    add_components(blade, blade_internal_structure, I_spar_hp, I_spar_lp)
    
    blade.updateBlade()
    # save(blade_name)
    # BladeDef_to_NuMADfile(obj,numad_name,matdb_name,numad_af_folder)


def add_stations(blade,blade_outer_shape_bem, hub_outer_shape_bem,
                    af_data, file: str):

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

        if blade.write_airfoils:
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
        af.resample()
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
    pass


def add_materials(blade, mat_data):
    for i in range(len(mat_data)):
        cur_mat = Material()
        cur_mat.name = mat_data[i]['name']
        if mat_data[i]['orth'] == 1:
            cur_mat.type = 'orthotropic'
        else:
            cur_mat.type = 'isotropic'
        # Add ply thickness option if ply thickness exists in yaml
        try:
            cur_mat.layerthickness = mat_data[i]['ply_t'] * 1000
        except KeyError:
                print('Warning! material ply thickness ' + 
                        mat_data[i]['name'] + 
                        ' not defined, assuming 1 mm thickness')
                cur_mat.layerthickness = 1
            
        finally:
            pass

        # first
        cur_mat.uts = _parse_data(mat_data[i]['Xt'])
        cur_mat.ucs = -_parse_data(mat_data[i]['Xc'])
        cur_mat.uss = _parse_data(mat_data[i]['S'])
        cur_mat.xzit = 0.3
        cur_mat.xzic = 0.25
        cur_mat.yzit = 0.3
        cur_mat.yzic = 0.25
        with np.errstate(divide='ignore', invalid='ignore'):
            cur_mat.g1g2 = np.divide(mat_data[i]['GIc'], mat_data[i]['GIIc'])
        cur_mat.alp0 = _parse_data(mat_data[i]['alp0'])
        cur_mat.etat = LARCetaT(cur_mat.alp0)
        try:
            #test if it is a list
            mat_data[i]['E']+[]
        except TypeError:
            cur_mat.ex = _parse_data(mat_data[i]['E'])
            cur_mat.ey = _parse_data(mat_data[i]['E'])
            cur_mat.ez = _parse_data(mat_data[i]['E'])
            cur_mat.gxy = _parse_data(mat_data[i]['G'])
            cur_mat.gxz = _parse_data(mat_data[i]['G'])
            cur_mat.gyz = _parse_data(mat_data[i]['G'])
            cur_mat.prxy = _parse_data(mat_data[i]['nu'])
            cur_mat.prxz = _parse_data(mat_data[i]['nu'])
            cur_mat.pryz = _parse_data(mat_data[i]['nu'])
            cur_mat.etal = LARCetaL(cur_mat.uss,cur_mat.ucs,cur_mat.alp0)
        else:
            cur_mat.ex = _parse_data(mat_data[i]['E'][0])
            cur_mat.ey = _parse_data(mat_data[i]['E'][1])
            cur_mat.ez = _parse_data(mat_data[i]['E'][2])
            cur_mat.gxy = _parse_data(mat_data[i]['G'][0])
            cur_mat.gxz = _parse_data(mat_data[i]['G'][1])
            cur_mat.gyz = _parse_data(mat_data[i]['G'][2])
            cur_mat.prxy = _parse_data(mat_data[i]['nu'][0])
            cur_mat.prxz = _parse_data(mat_data[i]['nu'][1])
            cur_mat.pryz = _parse_data(mat_data[i]['nu'][2])
            cur_mat.etal = LARCetaL(cur_mat.uss[0],cur_mat.ucs[1],cur_mat.alp0)
        try:
            cur_mat.m = mat_data[i]['m']
        except KeyError:
            print(f"No fatigue exponent found for material: {mat_data[i]['name']}")
        cur_mat.density = mat_data[i]['rho']
        # cur_mat.dens = mat_data[i]['rho']
        cur_mat.drydensity = mat_data[i]['rho']
        if 'description' in mat_data[i].keys() and 'source' in mat_data[i].keys():
            desc_sourc = [mat_data[i]['description'],', ',mat_data[i]['source']]
            cur_mat.reference = ''.join(desc_sourc)
        else:
            cur_mat.reference = []
        if blade.materials:
            blade.materials.append(cur_mat)
        else:
            blade.materials = []
            blade.materials.append(cur_mat)


def add_components(blade, blade_internal_structure, I_spar_hp, I_spar_lp):
    N_layer_comp = len(blade_internal_structure['layers'])
    for i in range(N_layer_comp):
        cur_comp = Component()
        cur_comp.group = 0
        cur_comp.name = blade_internal_structure['layers'][i]['name']
        #   comp['material'] = blade_internal_structure['layers']{i}['material'];
        mat_names = [mat.name for mat in blade.materials]
        C,IA,IB = np.intersect1d(mat_names,blade_internal_structure['layers'][i]['material'],return_indices=True)
        cur_comp.materialid = IA[0]
        try:
            cur_comp.fabricangle = np.mean(blade_internal_structure['layers'][i]['fiber_orientation']['values'])
        finally:
            pass
        if 'spar' in blade_internal_structure['layers'][i]['name'].lower():
            cur_comp.imethod = 'pchip'
        else:
            cur_comp.imethod = 'linear'
        # cur_comp.cp[:,0] = np.transpose(blade_internal_structure['layers'][i]['thickness']['grid'])
        cptemp1 = np.transpose(blade_internal_structure['layers'][i]['thickness']['grid'])
        temp_n_layer = np.multiply(np.transpose(blade_internal_structure['layers'][i]['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness
        I_round_up = np.flatnonzero((temp_n_layer > 0.05) & (temp_n_layer < 0.5))
        cptemp2 = np.round(np.multiply(np.transpose(blade_internal_structure['layers'][i]['thickness']['values']),1000.0) / blade.materials[cur_comp.materialid].layerthickness)
        cur_comp.cp = np.stack((cptemp1,cptemp2),axis=1)
        if I_round_up.size > 0:
            cur_comp.cp[I_round_up,1] = 1 # increase n_layers from 0 to 1 for 0.05<n_layers<0.5
        #     comp['cp'](:,2) = cell2mat(blade_internal_structure['layers']{i}['thickness']['values'])'.*1000;  # use when each material ply is 1 mm
        cur_comp.pinnedends = 0
        if blade.components:
            blade.components.append(cur_comp)
        else:
            blade.components = []
            blade.components.append(cur_comp)


    # Spar Caps (pressure and suction)
    blade.components[I_spar_hp].hpextents = ['b','c']
    blade.components[I_spar_lp].lpextents = ['b','c']
    
    for i in range(N_layer_comp):

        # uv coating
        if 'uv' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['le','te']
            blade.components[i].lpextents = ['le','te']
            blade.components[i].cp[:,1] = blade.components[i].cp[:,1]

        # Shell skin1
        if 'shell_skin_outer' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['le','te']
            blade.components[i].lpextents = ['le','te']
            # CK Change me when yaml is fixed!!!!
            blade.components[i].cp[:,1] = blade.components[i].cp[:,1]
        
        # LE Band
        if 'le_reinf' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['le','a']
            blade.components[i].lpextents = ['le','a']  
        
        # TE Band
        if 'te_reinf' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['d','te']
            blade.components[i].lpextents = ['d','te']
        
        # Trailing edge suction surface panel
        if 'te_ss' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].lpextents = ['c','d']
    
        # Leading edge suction surface panel
        if 'le_ss' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].lpextents = ['a','b']

        # Leading edge pressure surface panel)
        if 'le_ps' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['a','b']

        # Trailing edge pressure surface panel
        if 'te_ps' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = ['c','d']
    
        # Shell skin2
        if 'shell_skin_inner' in blade_internal_structure['layers'][i]['name'].lower():
            blade.components[i].hpextents = np.array(['le','te'])
            blade.components[i].lpextents = np.array(['le','te'])
            # CK Change me when yaml is fixed!!!!
            blade.components[i].cp[:,1] = blade.components[i].cp[:,1]

        # Forward Shear
        if 'web' in blade_internal_structure['layers'][i].keys():
            if 'fore' in blade_internal_structure['layers'][i]['web'].lower():
                # Web Skin1
                if 'skin_le' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z+' sw_offset]};
                    # comp[i].lpextents = {['z+' sw_offset]};
                    blade.components[i].hpextents = ['b']
                    blade.components[i].lpextents = ['b']
                    blade.components[i].group = 1
                    blade.components[i].name = blade_internal_structure['layers'][i]['name']
                    # CK Change me when yaml is fixed!!!!
                    blade.components[i].cp[:,1] = blade.components[i].cp[:,1]
                
                # Web Filler
                if 'filler' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z+' sw_offset]};
                    # comp[i].lpextents = {['z+' sw_offset]};
                    blade.components[i].hpextents = ['b']
                    blade.components[i].lpextents = ['b']
                    blade.components[i].group = 1
                    blade.components[i].name = blade_internal_structure['layers'][i]['name']

                # Web Skin2
                if 'skin_te' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z+' sw_offset]};
                    # comp[i].lpextents = {['z+' sw_offset]};
                    blade.components[i].hpextents = ['b']
                    blade.components[i].lpextents = ['b']
                    blade.components[i].group = 1
                    blade.components[i].name = np.array([blade_internal_structure['layers'][i]['name']])
                    # CK Change me when yaml is fixed!!!!
                    blade.components[i].cp[:,1] = blade.components[i].cp[:,1]

        # Rear Shear
        if 'web' in blade_internal_structure['layers'][i].keys():
            if 'rear' in blade_internal_structure['layers'][i]['web'].lower():
                # Web Skin1
                if 'skin_le' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z-' sw_offset]};
                    # comp[i].lpextents = {['z-' sw_offset]};
                    blade.components[i].hpextents = ['c']
                    blade.components[i].lpextents = ['c']
                    blade.components[i].group = 2
                    blade.components[i].name = blade_internal_structure['layers'][i]['name']
                    # CK Change me when yaml is fixed!!!!
                    blade.components[i].cp[:,1] = blade.components[i].cp[:,1]
                # Web Filler
                if 'filler' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z-' sw_offset]};
                    # comp[i].lpextents = {['z-' sw_offset]};
                    blade.components[i].hpextents = ['c']
                    blade.components[i].lpextents = ['c']
                    blade.components[i].group = 2
                    blade.components[i].name = blade_internal_structure['layers'][i]['name']
    
                # Web Skin2
                if 'skin_te' in blade_internal_structure['layers'][i]['name'].lower():
                    # comp[i].hpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].lpextents = {[num2str(xs,2) 'b-c']};
                    # comp[i].hpextents = {['z-' sw_offset]};
                    # comp[i].lpextents = {['z-' sw_offset]};
                    blade.components[i].hpextents = ['c']
                    blade.components[i].lpextents = ['c']
                    blade.components[i].group = 2
                    blade.components[i].name = blade_internal_structure['layers'][i]['name']
                    # CK Change me when yaml is fixed!!!!
                    blade.components[i].cp[:,1] = blade.components[i].cp[:,1]



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