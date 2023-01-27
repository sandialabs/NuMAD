import numpy as np
import pandas as pd
from os.path import abspath, dirname, join
from pynumad.objects.Component import Component
from pynumad.objects.Material import Material
from pynumad.utils.interpolation import interpolator_wrap

# DATA_DIR = os.path.dirname(os.path.abspath(__file__)) 
# todo: fix this to not be dependent on my proj structure
file_dir = dirname(abspath(str(__file__)))
data_dir = join(file_dir, "..","..","data")

def excel_to_blade(blade, filename = None):
    """
    xlsBlade  Construct BladeDef object with inputs from spreadsheet.  

    Parameters
    ----------
    filename : string
        path to Excelfile
    
    Returns
    -------
    blade : Blade
        blade object populated by Excel file
    
    Example
    -------
    blade = xlsBlade(FILENAME);

    See also BladeDef
    if not ('filename' is not None)  or len(filename)==0:
        raise Exception('xlsBlade:NoInput','No input file specified.')

    if not os.path.exist(str(filename)) :
        raise Exception('xlsBlade:FileNotFound','File {} not found.'.format(filename))
    """
    
    MPa_to_Pa = 1000000.0
    # Dictionary containing column indices
    xls_dict = {}

    xls_dict['geom'] = {}
    xls_dict['geom']['datarow1'] = 6
    xls_dict['geom']['span'] = 0
    xls_dict['geom']['twist'] = 1
    xls_dict['geom']['chord'] = 2
    xls_dict['geom']['thick'] = 3
    xls_dict['geom']['offset'] = 4
    xls_dict['geom']['aerocenter'] = 5
    xls_dict['geom']['afspan'] = 7
    xls_dict['geom']['afname'] = 8
    xls_dict['geom']['ispan'] = 10

    xls_dict['cmpt'] = {}
    xls_dict['cmpt']['paramcol'] = 2
    xls_dict['cmpt']['paramrow1'] = 1
    xls_dict['cmpt']['datarow1'] = 6
    xls_dict['cmpt']['group'] = 0
    xls_dict['cmpt']['name'] = 1
    xls_dict['cmpt']['matid'] = 2
    xls_dict['cmpt']['angle'] = 3
    xls_dict['cmpt']['hpext'] = 4
    xls_dict['cmpt']['lpext'] = 5
    xls_dict['cmpt']['cpspan'] = 6
    xls_dict['cmpt']['cpnlay'] = 7
    xls_dict['cmpt']['imethod'] = 8
    
    xls_dict['mtrl'] = {}
    xls_dict['mtrl']['datarow1'] = 3
    xls_dict['mtrl']['id'] = 0
    xls_dict['mtrl']['name'] = 1
    xls_dict['mtrl']['type'] = 2
    xls_dict['mtrl']['thickness'] = 3
    xls_dict['mtrl']['ex'] = 4
    xls_dict['mtrl']['ey'] = 5
    xls_dict['mtrl']['ez'] = 6
    xls_dict['mtrl']['gxy'] = 7
    xls_dict['mtrl']['gyz'] = 8
    xls_dict['mtrl']['gxz'] = 9
    xls_dict['mtrl']['prxy'] = 10
    xls_dict['mtrl']['pryz'] = 11
    xls_dict['mtrl']['prxz'] = 12
    xls_dict['mtrl']['density'] = 13
    xls_dict['mtrl']['drydensity'] = 14
    xls_dict['mtrl']['uts'] = 15
    xls_dict['mtrl']['ucs'] = 16
    xls_dict['mtrl']['reference'] = 17

    ## GEOMETRY

    # Read the Geometry tab of the xls file
    num = pd.read_excel(filename, sheet_name = 'Geometry', header = None)
    txt = pd.read_excel(filename, sheet_name = 'Geometry', dtype=str, header = None)

    if txt.iloc[1,1] == 'T':
        blade.naturaloffset = 1
    else:
        blade.naturaloffset = 0
    
    if txt.iloc[2,1] == 'CW':
        blade.rotorspin = 1
    else:
        blade.rotorspin = - 1
    
    if txt.iloc[3,1] == 'T':
        blade.swtwisted = 1
    else:
        blade.swtwisted = 0
    
    blade.span = np.array(num.iloc[xls_dict['geom']['datarow1']:,xls_dict['geom']['span']],dtype=float)
    # Next two lines handle the case where the spreadsheet tab has data
    # elsewhere on the tab below the last span value
    nans = np.nonzero(np.isnan(blade.span))
    if nans:
        lastrow = nans[0][0]
    else:
        lastrow = blade.span.shape[0]

    blade.span = blade.span[0:lastrow]
    if np.any(np.isnan(blade.span)):
        raise Exception('xlsBlade: span column must not have "holes" in data')
    
    lastrow = lastrow + xls_dict['geom']['datarow1']
    blade.degreestwist = np.array(num.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['twist']],dtype= float)
    blade.chord = np.array(num.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['chord']],dtype=float)
    blade.percentthick = np.array(num.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['thick']], dtype=float)
    blade.chordoffset = np.array(num.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['offset']], dtype=float)
    blade.aerocenter = np.array(num.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['aerocenter']], dtype=float)
    blade.sweep = np.zeros(blade.span.shape)
    blade.prebend = np.zeros(blade.span.shape)
    props = ['degreestwist','chord','percentthick','chordoffset','aerocenter']
    for prop in props:
        # For each of the input properties, interpolate where ever values
        # are missing.
        ind = np.isnan(getattr(blade,prop))
        if np.any(ind):
            if not prop == 'percentthick' :
                getattr(blade,prop)[ind] = interpolator_wrap(
                    np.delete(blade.span, ind),np.delete(getattr(blade,prop), ind),
                    blade.span[ind],'pchip'
                    )
            else:
                # jcb: note that blade.chord must be interpolated before
                # blade.percentthick
                absthick = np.multiply(blade.percentthick,blade.chord) / 100
                iabsthick = interpolator_wrap(
                    np.delete(blade.span, ind),np.delete(absthick,ind),blade.span[ind],'pchip'
                    )
                blade.percentthick[ind] = iabsthick / blade.chord[ind] * 100
            # The next two lines report when a property has been
            # interpolated.
            # rowstr = sprintf('#d,',find(ind==1));
            # fprintf('Interpolating "#s" on rows [#s]\n',props{k},rowstr(1:end-1))
    
    afspan = np.array(num.iloc[xls_dict['geom']['datarow1']:,xls_dict['geom']['afspan']],dtype=float)
    nans = np.nonzero(np.isnan(afspan))
    if nans:
        lastrow = nans[0][0]
    else:
        lastrow = blade.span.shape[0]
    afspan = afspan[0:lastrow]
    lastrow = lastrow + xls_dict['geom']['datarow1']
    afname = list(txt.iloc[xls_dict['geom']['datarow1']:lastrow,xls_dict['geom']['afname']])
    for k in range(len(afspan)):
        if afspan[k] < np.amin(blade.span) or afspan[k] > np.amax(blade.span):
            raise Exception('xlsBlade: location of airfoil #%d is outside given span distribution',k)
        affile = data_dir + '/airfoils/{}.txt'.format(afname[k])
        blade.addStation(affile,afspan[k])
        blade.stations[-1].airfoil.resample(175,'cosine')
    
    # afdb = np.array([blade.stations.airfoil])
    # afdb.resample(175,'cosine')
    blade.ispan = np.array(num.iloc[xls_dict['geom']['datarow1']:,xls_dict['geom']['ispan']],dtype=float)
    nans = np.nonzero(np.isnan(blade.ispan))
    try:
        lastrow = nans[0][0]
    except IndexError:
        lastrow = blade.ispan.shape[0]
    blade.ispan = blade.ispan[0:lastrow]

    ## COMPONENTS

    # Read the Components tab of the xls file
    num = pd.read_excel(filename, sheet_name = 'Components', header = None)
    txt = pd.read_excel(filename, sheet_name = 'Components', dtype=str, header = None)
    raw = pd.read_excel(filename, sheet_name = 'Components', dtype=object, header = None)

    blade.sparcapwidth = num.iloc[xls_dict['cmpt']['paramrow1'],xls_dict['cmpt']['paramcol']]
    blade.leband = num.iloc[xls_dict['cmpt']['paramrow1'] + 1,xls_dict['cmpt']['paramcol']]
    blade.teband = num.iloc[xls_dict['cmpt']['paramrow1'] + 2,xls_dict['cmpt']['paramcol']]
    # ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    blade.sparcapoffset = num.iloc[xls_dict['cmpt']['paramrow1'] + 3,xls_dict['cmpt']['paramcol']]
    if np.isnan(blade.sparcapoffset):
        blade.sparcapoffset = 0
    
    # ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    blade.components = []
    N = num.shape[0] - xls_dict['cmpt']['datarow1']
    for k in range(N):
        comp = Component()
        comp.group = num.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['group']]
        comp.name = txt.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['name']]
        comp.materialid = num.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['matid']]
        comp.fabricangle = readnumlist(raw.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['angle']])
        comp.hpextents = readstrlist(txt.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['hpext']])
        comp.lpextents = readstrlist(txt.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['lpext']])
        comp.cp = readnumlist(raw.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['cpspan']])
        comp.cp = np.stack((comp.cp,readnumlist(raw.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['cpnlay']])),axis=1)
        comp.imethod = txt.iloc[xls_dict['cmpt']['datarow1'] + k,xls_dict['cmpt']['imethod']]
        if not np.any(len(comp.hpextents) == np.array([0,1,2])) :
            raise Exception('xlsBlade: component #%d, length of hpextents must be 0, 1, or 2',k + 1)
        if not np.any(len(comp.lpextents) == np.array([0,1,2])) :
            raise Exception('xlsBlade: component #%d, length of lpextents must be 0, 1, or 2',k + 1)
        blade.components.append(comp)
    
    ## MATERIALS

    # Read the Materials tab of the xls file
    num = pd.read_excel(filename, sheet_name = 'Materials', header = None)
    txt = pd.read_excel(filename, sheet_name = 'Materials', dtype=str, header = None)

    blade.materials = []
    N = num.shape[0] - xls_dict['mtrl']['datarow1']
    for k in range(N):
        mat = Material()
        mat.name = txt.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['name']]
        mat.type = txt.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['type']]
        mat.layerthickness = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['thickness']]
        mat.ex = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['ex']]
        mat.prxy = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['prxy']]
        mat.density = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['density']]
        mat.drydensity = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['drydensity']]
        mat.uts = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['uts']]
        mat.ucs = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['ucs']]
        mat.reference = txt.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['reference']]
        if mat.type=='orthotropic':
            mat.ey = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['ey']]
            mat.ez = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['ez']]
            mat.gxy = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['gxy']]
            mat.gyz = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['gyz']]
            mat.gxz = MPa_to_Pa * num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['gxz']]
            mat.pryz = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['pryz']]
            mat.prxz = num.iloc[xls_dict['mtrl']['datarow1'] + k,xls_dict['mtrl']['prxz']]
        else:
            mat.ey = []
            mat.ez = []
            mat.gxy = []
            mat.gyz = []
            mat.gxz = []
            mat.pryz = []
            mat.prxz = []
        if not mat.ez:
            mat.ez = mat.ey
        if not mat.gyz:
            mat.gyz = mat.gxy
        if not mat.gxz:
            mat.gxz = mat.gyz
        if not mat.pryz:
            mat.pryz = mat.prxy
        if not mat.prxz:
            mat.prxz = mat.prxy
        blade.materials.append(mat)
    
    return blade
    
    
def readnumlist(strings): 
    # read a list of numeric values

    # check if numeric
    try:
        strings - 1
    # handle string list case
    except TypeError:
        #drop brackets on either end '[string]' -> 'string'
        if strings[0] == '[':
            strings = strings[1:]
        if strings[-1] == ']':
            strings = strings[:-1]
        # get list of numbers
        strings = strings.split(',')
        strings = [float(num) for num in strings]
        # could return as array - not sure.
    return np.array(strings)

def readstrlist(strings = None):
    if len(strings) == 0:
        strlist = strings
    else:
        strlist = strings.split(',')
    return strlist
    
    
"""
readnumlist scratch code

    # str = strreps(str,np.array(['[',']',',']),np.array(['','',' ']))
    # numv = cell2mat(textscan(str,'%f'))


def readstrlist(str = None): 
    # read a list of string values
    str = strreps(str,np.array(['[',']',',']),np.array(['','',' ']))
    if len(str)==0:
        str = ' '
    
    strvcell = textscan(str,'%s')
    strv = transpose(strvcell[0])
    return strv
    
    return blade

def strreps(strin = None,oldsubstrcell = None,newsubstrcell = None): 
    assert np.asarray(oldsubstrcell).size == np.asarray(newsubstrcell).size,'Lengths of substring cell arrays must be equal.'
    strout = strin
    for k in np.arange(1,np.asarray(oldsubstrcell).size+1).reshape(-1):
        strout = strout.replace(oldsubstrcell[k],newsubstrcell[k])
    
    return strout
"""