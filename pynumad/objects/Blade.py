########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import re, warnings
import numpy as np

from pynumad.io.yaml_to_blade import yaml_to_blade
from pynumad.io.excel_to_blade import excel_to_blade
from pynumad.utils.interpolation import interpolator_wrap
from pynumad.utils.affinetrans import rotation, translation
from pynumad.objects.Station import Station
from pynumad.objects.Airfoil import getAirfoilNormals, getAirfoilNormalsAngleChange
from pynumad.objects.Stack import Stack
from pynumad.objects.Subobjects import MatDBentry, Layer, Shearweb, BOM, Ply
from pynumad.shell.shell import shellMeshGeneral
# for type hints
from numpy import ndarray

class Blade():
    """BladeDef A class definition for wind & water turbine blades.
    
    Parameters
    ----------
    file : string
    write_afs : bool

    Attributes
    ----------
    aerocenter : array 
        Aerodynamic center of airfoil (used only by NuMAD->FAST)
    chord : array 
        Chord distribution [m] 
    chordoffset : array 		
        Chordwise offset (in addition to natural offset)
    components : list 
		Blade components such as spar, panels, etc., refer to ``ComponentDef``
    degreestwist : array 
		Twist distribution [degrees]
    ispan : array 
		Spanwise locations of interpolated output
    leband : float 
		Location of keypoint a
    materials : list 
		Material properties, refer to ``MaterialDef``
    mesh : float 
		Approximate element edge size for FE model [m]
    percentthick : array 
		Percent thickness of airfoil [%]
    prebend : array 
		Blade prebend, reference axis location along x2 [m]
    span : array 
		Spanwise location of distributed properties [m]
    sparcapoffset : array 
		(Does Nothing)
    sparcapwidth : array 
		Locations of keypoints b & c, defines distance 
        between keypoints b & c [mm]. First entry is the HP spar cap. 
        Second entry is the LP spar cap
    stations : list 
		Blade Stations, define the camber and thickness along the blade, 
        refer to ``StationDef``
    sweep : array 
		Blade Sweep, Reference axis location along x1 [m] 
    self.teband
    idegreestwist : array 
		interpolated twist
    ichord : array 
		interpolated chord
    ipercentthick : array 
		interpolated thickness
    self.ic = array
    self.icamber = array
    self.ithickness = array
    ichordoffset : array 
		interpolated offset
    iaerocenter : array 
		interpolated aerocenter
    isweep : array 
		interpolated sweep
    iprebend : array 
		interpolated prebend
    xoffset : array 
		natural offset
    profiles : array 
		normalized airfoil profiles
    geometry : array 
		actual x,y,z geometry
    arclength : array 
		surface distance from L.E.
    cpos : array 
		chordwise position
    LEindex : int  
    HParcx0 : array
    LParcx0 : array
    keylabels : list
    keypoints : array 
    keyarcs : array 
    keycpos : array 
    keyareas : array 
    LEbond : array 
    TEbond : array 
    webindices : list
    webpoints : list
    webarcs : list
    webcpos : list
    webareas : list
    webwidth : list
    webbonds : list
    bom : dict
    bomIndices : dict
    stacks : array 
		array of StackDef
    swstacks : list 
		contains StackDefs
    matdb : dict 
		Composite definition for each region at each station
    TEtype : list 
		trailing edge type; assigned in updateKeypoints
    shearweb : list
    bomPlot : dict
    hgGeometry : list
    hgKeypoints : list
    job_name : string
    paths : dict
    ansys : dict
        generate ANSYS settings
    write_airfoils : bool
    Example
    -------
    ``blade = BladeDef();``
    
    See also ``BladeDef_to_NuMADfile``, ``xlsBlade``, ``AirfoilDef``, 
    ``StationDef``, ``ComponentDef``, ``StackDef``
    """
    def __init__(self, filename=None, write_afs=False):

        self.aerocenter : ndarray = None
        self.chord : ndarray = None
        self.chordoffset : ndarray = None
        self.components : list = None
        self.degreestwist : ndarray = None
        self.ispan : ndarray = None
        self.leband : float = None
        self.materials : list = None
        self.mesh : float = 0.45
        self.percentthick : ndarray = None
        self.prebend : ndarray = None
        self.span : ndarray = None
        self.sparcapoffset : ndarray = None
        self.sparcapwidth : ndarray = None 
        self.stations : list = None 
        self.sweep : ndarray = None
        self.teband : float = None
        self.idegreestwist : ndarray = None
        self.ichord : ndarray = None
        self.ipercentthick : ndarray = None
        self.ic : ndarray = None
        self.icamber : ndarray = None
        self.ithickness : ndarray = None
        self.ichordoffset : ndarray = None
        self.iaerocenter : ndarray = None
        self.isweep : ndarray = None
        self.iprebend : ndarray = None
        self.xoffset : ndarray = None
        self.profiles : ndarray = None
        self.geometry : ndarray = None
        self.arclength : ndarray = None
        self.cpos : ndarray = None
        self.LEindex : int = None
        self.HParcx0 : ndarray = None
        self.LParcx0 : ndarray = None
        self.keylabels : list = None
        self.keypoints : ndarray = None
        self.keyarcs : ndarray = None
        self.keycpos : ndarray = None
        self.keyareas : ndarray = None
        self.LEbond : ndarray = None
        self.TEbond : ndarray = None
        self.webindices : list = None
        self.webpoints : list = None
        self.webarcs : list = None
        self.webcpos : list = None
        self.webareas : list = None
        self.webwidth : list = None
        self.webbonds : list = None
        self.bom : dict = None
        self.bomIndices : dict = None
        self.stacks : ndarray = None
        self.swstacks : list = None
        self.hgGeometry : list = None
        self.hgKeypoints : list = None
        self.matdb : dict = None
        self.TEtype : list = None 
        self.shearweb : list = None
        self.bomPlot : dict = {
            'kLayer':1, 'hgLinesHP':[], 'hgLinesLP':[],
            'hgPatchHP':[], 'hgPatchLP':[],
            'uisliderHP':[], 'uisliderLP':[],
            'hTitleHP':[], 'hTitleLP':[]
            }
        self.job_name : str = 'numad.nmd'
        self.paths : dict = {
            'job':'', 'numad':'', 'precomp':'', 'bmodes':'',
            'ansys':'','batch_run':0
            }
        self.ansys : dict = {
            'BoundaryCondition':'','ElementSystem':'',
            'MultipleLayerBehavior':'', 'meshing':'','smartsize':[],
            'elementsize':[], 'shell7gen':[],
            'dbgen':[], 'FailureCriteria':[]
            }
        
        #properties
        self._naturaloffset = 1  # 1 = offset by max thickness location, 0= do not offset to max thickness
        self._rotorspin = 1 # Rotor Spin, 1= CW rotation looking downwind, -1= CCW rotation
        self._swtwisted = 0 # Shear Web, 0 = planar shear webs, 1= shear webs twisted by blade twist
        self.write_airfoils = write_afs # To determine whether to write airfoil files when parsing yaml
        try:
            if 'yaml' in filename or 'yml' in filename:
                self.read_yaml(filename)
            elif 'xls' in filename or 'xlsx' in filename:
                self.read_excel(filename)
            else:
                raise Exception('Unknown filetype. Currently supported inputs are excel and yaml files.')
        # To handle when filename == None
        except TypeError:
            pass


    @property
    def naturaloffset(self):
        """
        TODO docstring
        """
        return self._naturaloffset


    @naturaloffset.setter
    def naturaloffset(self, new_naturaloffset):
        """
        TODO docstring
        """
        if not (new_naturaloffset==0 or new_naturaloffset==1):
            raise Exception('naturaloffset must be 0 or 1')
        else:
            self._naturaloffset = new_naturaloffset

    
    @property
    def rotorspin(self):
        """
        TODO docstring
        """
        return self._rotorspin


    @rotorspin.setter
    def rotorspin(self,new_rotorspin):
        """
        TODO docstring
        """   
        if not (new_rotorspin==1 or new_rotorspin==-1):
            raise Exception('rotorspin must be 1 (cw) or -1 (ccw)')
        else:
            self._rotorspin = new_rotorspin


    @property
    def swtwisted(self):
        """
        TODO docstring
        """
        return self._swtwisted

    @swtwisted.setter
    def swtwisted(self,new_swtwisted):
        """
        TODO docstring
        """  
        if not (new_swtwisted==0 or new_swtwisted==1):
            raise Exception('swtwisted must be 0 or 1')
        else:
            self._swtwisted = new_swtwisted

    
    def read_yaml(self, filename):
        """Populate blade attributes with yaml file data

        Extended description of function.

        Parameters
        ----------
        filename: str
            name of yaml file to be read

        Returns
        -------
        None

        """
        yaml_to_blade(self, filename)


    def read_excel(self, filename):
        """Populate blade attributes with excel file data

        Extended description of function.

        Parameters
        ----------
        filename: str
            name of excel file to be read

        Returns
        -------
        None

        """
        excel_to_blade(self, filename)


    def updateBlade(self):
        """
        TODO docstring
        """
        self.updateGeometry()
        self.updateKeypoints()
        self.updateBOM()
        return


    def updateGeometry(self):
        """This method updates the interpolated blade parameters
        """
        
        # update the interpolated station profiles
        nStations = len(self.stations)
        if nStations > 0:
            nPoints = len(self.stations[1].airfoil.c)
        else:
            raise Exception(
                'BladeDef must have at least one station before updating geometry.')
        
        # ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # add some error checking -- first station must be at blade
        # root to prevent extrapolation
        #   assert(self.stations(1).spanlocation==0,'first station must be at the blade root')
        # ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        # Collect parameter tables from the stations.
        spanlocation = np.array([self.stations[i].spanlocation for i in range(len(self.stations))])
        c = np.zeros((nPoints,nStations))
        camber = np.zeros((nPoints,nStations))
        thickness = np.zeros((nPoints,nStations))
        tetype = [None]*nStations
        for k in range(0,nStations):
            ck = self.stations[k].airfoil.c
            if len(ck) != nPoints:
                raise Exception(
                    'Station airfoils must have same number of samples.')
            c[:,k] = ck
            camber[:,k] = self.stations[k].airfoil.camber
            thickness[:,k] = self.stations[k].airfoil.thickness
            tetype[k] = self.stations[k].airfoil.TEtype
        
        # ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # fix numerical issue due to precision on camber calculation
        # camber should start and end at y-values of zero
        camber[0,:] = np.zeros((1,nStations))
        camber[-1,:] = np.zeros((1,nStations))
        # ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        # Interpolate the station parameter tables.
        # Each column corresponds to an interpolated station.
        ic = interpolator_wrap(spanlocation,c,self.ispan,'pchip', axis=1)
        icamber = interpolator_wrap(spanlocation,camber,self.ispan,'pchip', axis=1)
        ithickness = interpolator_wrap(spanlocation,thickness,self.ispan,'pchip', axis=1)
        self.ic = ic #Export for TE opening
        self.icamber = icamber
        self.ithickness = ithickness
        self.cpos = np.concatenate((-ic[-1,:].reshape(1,-1),-np.flipud(ic),
            ic[1:,:],ic[-1,:].reshape(1,-1)),axis=0)
        
        # figure(101); surf(repmat(spanlocation,nPoints,1),c,camber,'MeshStyle','column');
        # figure(102); surf(repmat(spanlocation,nPoints,1),c,thickness,'MeshStyle','column');
        # figure(103); surf(repmat(self.ispan,nPoints,1),ic,icamber,'MeshStyle','column');
        # figure(104); surf(repmat(self.ispan,nPoints,1),ic,ithickness,'MeshStyle','column');
        
        # Adjust the thickness profiles based on TEtype of stations.
        # This is mainly for transitions to flatbacks were the
        # interpolated airfoil needs to look like a round.
        for k in range(len(self.ispan)):
            try:
                ind = np.argwhere(self.ispan[k] < spanlocation)[0][0]
                #maybe better: ind = np.flatnonzero(self.ispan[k] < spanlocation)[0]
            except:
                continue
            else:
                if ind == 1:
                    continue
            if tetype[ind]=='flat' and tetype[ind - 1]=='round':
                ithickness[-1,k] = 0
        
        # Interpolate the blade parameter curves.
        # Results are row vectors.
        self.idegreestwist = interpolator_wrap(self.span,self.degreestwist,self.ispan,'pchip')
        self.ichord = interpolator_wrap(self.span,self.chord,self.ispan,'pchip')
        absolutethick = np.multiply(self.percentthick,self.chord) / 100
        iabsolutethick = interpolator_wrap(self.span,absolutethick,self.ispan,'pchip')
        self.ipercentthick = iabsolutethick / self.ichord * 100

        # ensure that the interpolation doesn't reduce the percent
        #   thickness beneath the thinnest airfoil
        self.ipercentthick[self.ipercentthick < np.amin(self.percentthick)] = np.amin(self.percentthick)
        self.ichordoffset = interpolator_wrap(self.span,self.chordoffset,self.ispan,'pchip')
        self.iaerocenter = interpolator_wrap(self.span,self.aerocenter,self.ispan,'pchip')
        if len(self.sweep)==0:
            self.sweep = np.zeros((self.span.shape,self.span.shape)) 
        if len(self.prebend)==0:
            self.prebend = np.zeros((self.span.shape,self.span.shape))
        self.isweep = interpolator_wrap(self.span,self.sweep,self.ispan,'pchip')
        self.iprebend = interpolator_wrap(self.span,self.prebend,self.ispan,'pchip')
        
        # Generate the blade surface geometry.
        N = np.asarray(self.ispan).size
        M = nPoints * 2 + 1
        self.profiles = np.zeros((M,2,N))
        self.geometry = np.zeros((M,3,N))
        self.xoffset = np.zeros((1,N))
        self.LEindex = nPoints
        
        for k in range(0,N):
            self.updateAirfoilProfile(k)
            mtindex = np.argmax(ithickness[:,k])
            self.xoffset[0,k] = ic[mtindex,k]
            self.updateOMLgeometry(k)
        
        # Calculate the arc length of each curve
        self.arclength = np.zeros((M,N))
        self.HParcx0 = np.zeros((1,N))
        self.LParcx0 = np.zeros((1,N))
            # ptind = 1:M;  # point indices
            # ptindover = 1:0.2:M;  # over-sample point indices
            # LE = find(ptindover==self.LEindex);
        LE = self.LEindex
        for k in range(N):
            # xx = spline(ptind,self.geometry(:,1,k),ptindover);
            # yy = spline(ptind,self.geometry(:,2,k),ptindover);
            # zz = spline(ptind,self.geometry(:,3,k),ptindover);
            xx = self.geometry[:,0,k]
            yy = self.geometry[:,1,k]
            zz = self.geometry[:,2,k]
            arclen = np.sqrt(np.diff(xx)**2 + np.diff(yy)**2 + np.diff(zz)**2)
            arclen = np.concatenate((np.array([0]),np.cumsum(arclen)),axis=0)
            # self.arclength(:,k) = interpolator_wrap(ptindover,arclenover,ptind);
            self.arclength[:,k] = arclen
            LEarcsum = self.arclength[self.LEindex,k]
            self.arclength[:,k] = self.arclength[:,k] - LEarcsum
            # find where x=0 intersects the surface
            self.HParcx0[0,k] = interpolator_wrap(xx[1:LE+1],arclen[1:LE+1],0) - LEarcsum
            self.LParcx0[0,k] = interpolator_wrap(xx[-2:LE-1:-1],arclen[-2:LE-1:-1],0) - LEarcsum
        
        return


    def updateKeypoints(self):
        """This method updates the keypoints (a,b,c,...) which define the blade
        regions.

        Returns
        -------

        None

        Example:
          ``blade.updateKeypoints``

        find the curves which bound each blade region
        """
        
        #NOTE this might not work right with ispan taken as a 2D array
        N = self.ispan.size # number of interpolated span stations
        M = 12 # number of areas around airfoil profile; must be even (see calc of web areas)

        self.keypoints = np.zeros((M - 2,3,N)) # keypoints in xyz geometry
        self.keyarcs = np.zeros((M + 1,N)) # surface arclength distance of keypoints from LE
        self.keycpos = np.zeros((M + 1,N)) # chordwise position of keypoints
        self.keyareas = np.zeros((M,N - 1)) # surface area of regions created by keypoints
        self.keylabels = ['te','e','d','c','b','a','le',
                           'a','b','c','d','e','te']
        self.LEbond = np.zeros((N - 1))
        self.TEbond = np.zeros((N - 1))
        mm_to_m = 0.001

        # start and finish indices in geometry/arcs
        ns = 1
        nf = self.geometry.shape[0] - 2

        n1 = mm_to_m * self.leband # no foam width
        n2 = mm_to_m * self.teband # no foam width

        #keypoints, keyarcs, keycpos
        for k in range(0,N):
            # allow for separate definitions of HP and LP spar cap
            # width and offset [HP LP]
            if len(self.sparcapwidth) > 2 or len(self.sparcapoffset) > 2:
                raise Exception('too many entries for spar cap definition')
            scwidth_hp = mm_to_m * self.sparcapwidth[0] #type: float
            scwidth_lp = mm_to_m * self.sparcapwidth[-1] #type: float
            scoffset_hp = mm_to_m * self.sparcapoffset[0] #type: float
            scoffset_lp = mm_to_m * self.sparcapoffset[-1] #type: float

            tempTE = self.getprofileTEtype(k)
            if self.TEtype:
                self.TEtype.append(tempTE)
            else:
                self.TEtype = []
                self.TEtype.append(tempTE)
            if self.swtwisted:
                # get angle of each xy pair w.r.t. pitch axis (0,0)
                xyangle = np.zeros(self.geometry.shape[0])
                for j in range(len(xyangle)):
                    xy = self.geometry[j,0:2,k]
                    xyangle[j] = np.arctan2(self.rotorspin * xy[1],xy[0])
                # unwrap and center around 0
                xyangle = np.unwrap(xyangle)
                xyangle = xyangle - np.pi * np.round(xyangle[self.LEindex] / np.pi)

            k_arclen = self.arclength[ns:nf+1,k]
            k_geom = self.geometry[ns:nf+1,:,k]
            k_cpos = self.cpos[ns:nf+1,k]

            # ==================== HP surface ====================
            if self.swtwisted:
                # find arclength where xyangle equals normal to chord
                twistnorm = np.pi / 180 * (- self.idegreestwist[k] - 90) # angle normal to chord line
                z = interpolator_wrap(xyangle[ns:nf+1],k_arclen,twistnorm)
            else:
                z = self.HParcx0[0,k]
            z0 = z
            z = z - scoffset_hp
            a = np.amax(((0 - n1),0.1 * self.arclength[ns,k])) #type: float
            a = np.amin((a,0.01 * self.arclength[ns,k]))
            b = np.amin(((z + 0.5 * scwidth_hp),0.15 * self.arclength[ns,k]))
            c = np.amax(((z - 0.5 * scwidth_hp),0.8 * self.arclength[ns,k]))
            d = np.amin(((self.arclength[0,k] + n2),0.85 * self.arclength[ns,k]))
            d = np.amax((d,0.98 * self.arclength[ns,k]))
            if str(self.TEtype[k]) == 'flat':
                e = self.arclength[ns,k]
                self.keypoints[0,:,k] = self.geometry[ns,:,k]
                self.keycpos[1,k] = -1
            else:
                # e = 0.5 * (d + self.arclength(ns,k));
                e = 0.99 * self.arclength[ns,k]
                self.keypoints[0,:,k] = interpolator_wrap(k_arclen,k_geom,e)
                self.keycpos[1,k] = interpolator_wrap(k_arclen,k_cpos,e)
            # 1 -> e
            self.keypoints[1,:,k] = interpolator_wrap(k_arclen,k_geom,d)
            self.keypoints[2,:,k] = interpolator_wrap(k_arclen,k_geom,c)
            # self.keypoints(  ,:,k) = interpolator_wrap(self.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            self.keypoints[3,:,k] = interpolator_wrap(k_arclen,k_geom,b)
            self.keypoints[4,:,k] = interpolator_wrap(k_arclen,k_geom,a)
            self.keyarcs[0,k] = self.arclength[ns,k]
            self.keyarcs[1,k] = e
            self.keyarcs[2,k] = d
            self.keyarcs[3,k] = c
            # self.keyarcs(  ,k)   = z;
            self.keyarcs[4,k] = b
            self.keyarcs[5,k] = a
            self.keyarcs[6,k] = 0 # le
            self.keycpos[0,k] = self.cpos[ns,k] #te, hp surface
            #            2   -> e
            self.keycpos[2,k] = interpolator_wrap(k_arclen,k_cpos,d)
            self.keycpos[3,k] = interpolator_wrap(k_arclen,k_cpos,c)
            #                 self.keycpos(  ,k) = interpolator_wrap(self.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            self.keycpos[4,k] = interpolator_wrap(k_arclen,k_cpos,b)
            self.keycpos[5,k] = interpolator_wrap(k_arclen,k_cpos,a)
            self.keycpos[6,k] = interpolator_wrap(k_arclen,k_cpos,0)

            # ==================== LP surface ====================
            if self.swtwisted:
                twistnorm = np.pi / 180 * (-self.idegreestwist[k] + 90) # angle normal to chord line
                z = interpolator_wrap(xyangle[ns:nf+1],k_arclen,twistnorm)
            else:
                z = self.LParcx0[0,k]
            z0 = z # ble: location where airfoil surface crosses Xglobal=0
            z = z + scoffset_lp # positive scoffset moves z toward t.e.
            a = np.amin(((0 + n1),0.1 * self.arclength[nf,k]))
            a = np.amax((a,0.01 * self.arclength[nf,k]))
            b = np.amax(((z - 0.5 * scwidth_lp),0.15 * self.arclength[nf,k]))
            c = np.amin((z + 0.5 * scwidth_lp,0.8 * self.arclength[nf,k]))
            d = np.amax((self.arclength[-1,k] - n2,0.85 * self.arclength[nf,k]))
            d = np.amin((d,0.96 * self.arclength[nf,k]))
            if str(self.TEtype[k]) == str('flat'):
                e = self.arclength[nf,k]
                self.keypoints[9,:,k] = self.geometry[nf,:,k]
                self.keycpos[11,k] = 1
            else:
                # e = 0.5 * (d + self.arclength(nf,k));
                e = 0.98 * self.arclength[nf,k]
                self.keypoints[9,:,k] = interpolator_wrap(k_arclen,k_geom,e)
                self.keycpos[11,k] = interpolator_wrap(k_arclen,k_cpos,e)
            self.keypoints[5,:,k] = interpolator_wrap(k_arclen,k_geom,a)
            self.keypoints[6,:,k] = interpolator_wrap(k_arclen,k_geom,b)
            # self.keypoints(  ,:,k) = interpolator_wrap(self.arclength(ns:nf,k),self.geometry(ns:nf,:,k),z);
            self.keypoints[7,:,k] = interpolator_wrap(k_arclen,k_geom,c)
            self.keypoints[8,:,k] = interpolator_wrap(k_arclen,k_geom,d)
            # 10   -> e
            self.keyarcs[7,k] = a
            self.keyarcs[8,k] = b
            # self.keyarcs( ,k)   = z;
            self.keyarcs[9,k] = c
            self.keyarcs[10,k] = d
            self.keyarcs[11,k] = e
            self.keyarcs[12,k] = self.arclength[nf,k]
            self.keycpos[7,k] = interpolator_wrap(k_arclen,k_cpos,a)
            self.keycpos[8,k] = interpolator_wrap(k_arclen,k_cpos,b)
            # self.keycpos(  ,k) = interpolator_wrap(self.arclength(ns:nf,k),self.cpos(ns:nf,k),z);
            self.keycpos[9,k] = interpolator_wrap(k_arclen,k_cpos,c)
            self.keycpos[10,k] = interpolator_wrap(k_arclen,k_cpos,d)
            # 12   -> e
            self.keycpos[12,k] = self.cpos[nf,k] # te, lp surface
        
        # find the points used by each shear web
        cmpt_groups = [self.components[i].group for i in range(len(self.components))]
        #NOTE I am not sure if I translated the next two lines correctly -kb
        uniq_groups = np.unique(cmpt_groups)
        uniq_groups = uniq_groups[uniq_groups != 0] # group "0" is the blade skins
        self.webindices = []
        self.webarcs = []
        self.webcpos = []
        self.webpoints = []
        self.webareas = []
        self.webwidth = []
        self.webbonds = []
        for ksw in uniq_groups - 1: # for each shear web
            # pre-allocating arrays
            self.webindices.append([])
            self.webarcs.append(np.ndarray((2,N)))
            self.webcpos.append(np.ndarray((2,N)))
            self.webpoints.append(np.ndarray((2,3,N)))
            self.webareas.append(np.ndarray((N-1)))
            self.webwidth.append(np.ndarray((N)))
            self.webbonds.append(np.ndarray((2,N-1)))

            ksw_cmpts = np.argwhere(ksw+1 == cmpt_groups).reshape(-1) # find the components that are part of the shear web
            hpextents = np.unique([self.components[i].hpextents for i in ksw_cmpts]).tolist() # get the hp extents
            lpextents = np.unique([self.components[i].lpextents for i in ksw_cmpts]).tolist() # get the lp extents
            assert len(hpextents) == 1,f'HP Extents for components in group {ksw} must be identical and contain no spaces or commas'
            assert len(lpextents) == 1,f'LP Extents for components in group {ksw} must be identical and contain no spaces or commas'
            # match extents that have form of either '0.5b-c' or
            # 'b+/-100' or 'b' or 'z+/-100'
            # pat = '(?<fraction>\d*[\.]?\d*)(?<pt1>[a-zA-Z]+)-(?<pt2>[a-zA-Z]+)|(?<pt3>[a-zA-Z]+)(?<mm_offset>[+-]\d+)|(?<pt>[a-zA-Z])'
            pat = '(?P<fraction>\d*[\.]?\d*)(?P<pt1>[a-zA-Z]+)-(?P<pt2>[a-zA-Z]+)|(?P<pt3>[a-zA-Z]+)(?P<mm_offset>[+-]\d+)|(?P<pt>[a-zA-Z])'
            
            hp = re.search(pat,hpextents[0]).groupdict()
            lp = re.search(pat,lpextents[0]).groupdict()
            try:
                le = self.keylabels.index('le')
            except:
                print(f"HP extent label \"{hp['pt']}\" not defined.")
            # get shear web placement on HP side
            if hp['pt']:
                try:
                    n = self.keylabels[0:le+1].index(hp['pt'])## EMA
                except:
                    print(f"HP extent label \"{hp['pt']}\" not defined.")
                self.webindices[ksw].append(n)
                self.webarcs[ksw][0,:] = self.keyarcs[n,:]
                self.webcpos[ksw][0,:] = self.keycpos[n,:]
                n = n - 1
                self.webpoints[ksw][0,:,:] = self.keypoints[n,:,:]
            elif hp['pt1']:
                #NOTE not sure what hp['fraction'] looks like, but hopefully this line works -kb
                f = float(hp['fraction'])
                if f <= 0 or f >= 1:
                    raise Exception(f'Component group {ksw}: HP extent fraction={f}, which is outside range (0..1)')
                try:
                    n1 = self.keylabels[0:le+1].index(hp['pt1'])
                except:
                    print(f"HP extent label \"{hp['pt1']}\" not defined.")
                try:
                    n2 = self.keylabels[0:le+1].index(hp['pt2'])
                except:
                    print(f"HP extent label \"{hp['pt2']}\" not defined.")
                self.webindices[ksw].append(np.nan)
                p1 = self.keyarcs[n1,:]
                p2 = self.keyarcs[n2,:]
                p = (1-f)*p1 + f*p2
                self.webarcs[ksw][0,:] = p
                for k in range(N):
                    self.webcpos[ksw][0,k] = interpolator_wrap(k_arclen,k_cpos,p[k])
                    self.webpoints[ksw][0,:,k] = interpolator_wrap(k_arclen,k_geom,p[k])
            elif hp['pt3']:
                try:
                    n3 = self.keylabels[0:le+1].index(hp['pt3'])
                except:
                    print(f"HP extent label \"{hp['pt3']}\" not defined.")
                self.webindices[ksw].append(np.nan)
                p3 = self.keycpos[n3,:]
                p = p3 - float(hp['mm_offset']) / 1000
                iMax = self.keylabels[0,le+1].index('d')
                #NOTE potential for error here - array shapes TBD -kb
                pMax = np.multiply(self.keycpos[iMax,:],np.transpose(self.ichord))
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = self.keylabels[0:le+1].index('a')
                #NOTE same issue here -kb
                pMin = np.multiply(self.keycpos[iMin,:],np.transpose(self.ichord))
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                self.webcpos[ksw][0,:] = p
                for k in range(N):
                    self.webarcs[ksw][0,k] = interpolator_wrap(self.cpos[ns:nf+1,:,k],self.arclength[ns:nf+1,:,k],p[k])
                    self.webpoints[ksw][0,:,k] = interpolator_wrap(k_cpos,k_geom,p[k])
                # ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                #NOTE confused about these lines.-kb
                # p3Save[ksw][0,:] = - 1.0 * p3
                # pSave[ksw][0,:] = - 1.0 * p
                # ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            else:
                raise Exception('Shear web geometry HP extents not defined correctly (e.g., 0.5b-c, b, b+200)')
            # get shear web placement on LP side
            if lp['pt']:
                try:
                    n = self.keylabels[le:].index(lp['pt']) + le
                    self.webindices[ksw].append(n)
                    self.webarcs[ksw][1,:] = self.keyarcs[n,:]
                    self.webcpos[ksw][1,:] = self.keycpos[n,:]
                    self.webpoints[ksw][1,:,:] = self.keypoints[n,:,:]
                except:
                    print(f"LP extent label \"{lp['pt']}\" not defined.")

            elif lp['pt1']:
                f = float(lp['fraction'])
                if f < 0 or f > 1:
                    raise Exception(f'Component group {ksw}: LP extent fraction={f}, which is outside range [0..1]')
                try:
                    n1 = self.keylabels[le:].index(lp['pt1']) + le
                except:
                    print(f"LP extent label \"{lp['pt1']}\" not defined.")
                try:
                    n2 = self.keylabels[le:].index(lp['pt2']) + le 
                except:
                    print(f"LP extent label \"{lp['pt2']}\" not defined.")
                self.webindices[ksw].append(np.nan)
                p1 = self.keyarcs[n1,:]
                p2 = self.keyarcs[n2,:]
                p = (1 - f) * p1 + f * p2
                self.webarcs[ksw][1,:] = p
                #NOTE stopped translating here
                for k in range(N):
                    self.webcpos[ksw][1,k] = interpolator_wrap(k_arclen,k_cpos,p[k])
                    self.webpoints[ksw][1,:,k] = interpolator_wrap(k_arclen,k_geom,p[k])
            elif lp['pt3']: 
                try:
                    n3 = self.keylabels[le:].index(lp['pt3']) + le
                except:
                    print(f"LP extent label \"{lp['pt3']}\" not defined.")       
                self.webindices[ksw].append(np.nan)
                p3 = self.keycpos[n3,:]
                p = p3 + float(lp['mm_offset']) / 1000
                iMax = self.keylabels[le:].index('d') + le
                pMax = np.multiply(self.keycpos[iMax,:],np.transpose(self.ichord))
                p[np.abs(p) > np.abs(pMax)] = pMax[np.abs(p) > np.abs(pMax)]
                iMin = self.keylabels[le:].index('a') + le
                pMin = np.multiply(self.keycpos[iMin,:],np.transpose(self.ichord))
                p[np.abs(p) < np.abs(pMin)] = pMin[np.abs(p) < np.abs(pMin)]
                self.webcpos[ksw][1,:] = p
                for k in range(N):
                    self.webarcs[ksw][1,k] = interpolator_wrap(k_cpos,k_arclen,p[k])
                    self.webpoints[ksw][1,:,k] = interpolator_wrap(k_cpos,k_geom,p[k])
                # ble <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                #NOTE again, not sure where these vars are coming from -kb
                # p3Save[ksw][1,:] = p3
                # pSave[ksw][1,:] = p
                # ble >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
            else:
                raise Exception('Shear web geometry LP extents not defined correctly (e.g., 0.5b-c, b, b+200)')
        
        # calculate shell areas
        for kc in range(N-1):
            for kr in range(M):
                # choose number of points to use in area calculation
                # jcb: I decided to base this on the number of points
                # in the interpolated station profile found within the region
                #  of interest.
                npts = sum(np.logical_and(self.arclength[:,kc]>= self.keyarcs[kr,kc],
                                          self.arclength[:,kc]<= self.keyarcs[kr+1,kc]))
                npts = np.amax((npts,2)) # need at least two points
                ibarc = np.linspace(self.keyarcs[kr,kc],self.keyarcs[kr+1,kc],npts) # inboard curve arclengths
                obarc = np.linspace(self.keyarcs[kr,kc+1],self.keyarcs[kr+1,kc+1],npts) # outboard curve arclengths
                ib = interpolator_wrap(self.arclength[ns:nf+1,kc],self.geometry[ns:nf+1,:,kc],ibarc) # inboard xyz
                ob = interpolator_wrap(self.arclength[ns:nf+1,kc+1],self.geometry[ns:nf+1,:,kc+1],obarc) # outboard xyz
                dspan = np.sqrt(np.sum((ob-ib)**2,1)) # "ds" in the span direction
                # treat each "rectangular" area as two triangles
                t1 = 0.5 * np.dot(np.sqrt(np.sum(np.diff(ib,1,axis=0) ** 2, 1)), dspan[0:-1])
                t2 = 0.5 * np.dot(np.sqrt(np.sum(np.diff(ob,1,axis=0) ** 2, 1)), dspan[1:])
                self.keyareas[kr,kc] = t1 + t2
                if kr==0:
                    self.TEbond[kc] = dspan[0]
                if (M / 2 + 1)==(kr + 1):
                    self.LEbond[kc] = dspan[0]
        
        # calculate areas used by shear webs
        # jcb: note that these areas come purely from the geometry and
        # do not take into account the thickness of the shell or
        # sparcap layup.
        for ksw in range(len(self.webpoints)):
            for kc in range(N-1):
                ib = self.webpoints[ksw][:,:,kc]
                ob = self.webpoints[ksw][:,:,kc+1]
                # treat each "rectangular" area as two triangles
                b1 = np.diff(ib,axis=0)
                b2 = np.diff(ob,axis=0)
                base1 = np.sqrt(np.sum(b1 ** 2, 1))[0]
                base2 = np.sqrt(np.sum(b2 ** 2, 1))[0]
                b1 = b1 / base1
                b2 = b2 / base2
                h1 = np.abs(np.dot((ob[0,:] - ib[0,:]), (1 - np.transpose(b1))))
                h2 = np.abs(np.dot((ib[1,:] - ob[1,:]), (1 - np.transpose(b2))))
                self.webareas[ksw][kc] = 0.5 * (base1*h1 + base2*h2)
                self.webwidth[ksw][kc] = base1
                # calculate edge (bond-line) lengths
                self.webbonds[ksw][0:2,kc] = np.sqrt(np.sum((ob - ib) ** 2, 1))
            self.webwidth[ksw][N-1] = base2
        
        return


    def updateBOM(self): 
        """This method updates the Bill-of-Materials
        See datatypes.BOM

        Returns
        -------

        None

        """

        # initialize structures
        self.bom = {'hp':[],'lp':[],'sw':[],'lebond':[],'tebond':[],
                    'swbonds':[],'dryweight':[]}
        self.bomIndices = {'hp':[],'lp':[],'sw':[]}
        # calculate non-dimensional span
        ndspan = (self.ispan - self.ispan[0]) / (self.ispan[-1] - self.ispan[0])
        hprow = 0
        lprow = 0
        swnum = None
        swrow = 0
        g_to_kg = 0.001
        m_to_mm = 1000.0
        mm_to_m = 0.001
        swBeginSta = []
        swEndSta = []
        for kc in range(len(self.components)):
            comp = self.components[kc]
            mat = self.materials[comp.materialid]
            hpRegion,lpRegion = self.findRegionExtents(comp)
            nlay = comp.getNumLayers(ndspan)
            nlay = np.round(nlay)
            # if 1
            layermult = 1 # multiplier=1 for individual layers
            # else
                # layermult = nlay(1);  # multiplies layer thickness
                # assert(all(nlay==nlay(1)),'updateBOM: component #d, ''T'' flag requires uniform layer thickness',kc);
                # nlay = ones(size(ndspan));
                # end
            for klay in range(1,int(np.max(nlay))+1):
                #NOTE: a little confused on beginSta endSta - are they
                # numbers, arrays? Assuming numbers for now...
                beginSta,endSta = self.findLayerExtents(nlay,klay)
                ## EMA original:
                # for ks = 1:length(beginSta)
                ## changed to :
                ksMax = np.amin((len(beginSta),len(endSta)))
                # situation that beginSta/endSta is longer than 1
                for ks in range(ksMax):
                    ## END
                    if comp.group == 0 and hpRegion:
                        areas = self.keyareas[hpRegion[0]:hpRegion[1],beginSta[ks]:endSta[ks]]
                        regionarea = sum(areas.flatten())
                        arcs = self.keyarcs[hpRegion[1],beginSta[ks]:endSta[ks]+1] - self.keyarcs[hpRegion[0],beginSta[ks]:endSta[ks]+1]
                        cur_bom = BOM()
                        cur_bom.layernum = hprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = self.ispan[beginSta[ks]]
                        cur_bom.endsta = self.ispan[endSta[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = layermult * mat.layerthickness
                        cur_bom.weight = layermult * mat.drydensity * regionarea
                        self.bomIndices['hp'].append([beginSta[ks],endSta[ks],*hpRegion])
                        self.bom['hp'].append(cur_bom)
                        hprow = hprow + 1
                        
                    if comp.group == 0 and lpRegion:
                        areas = self.keyareas[lpRegion[0]:lpRegion[1],beginSta[ks]:endSta[ks]]
                        regionarea = sum(areas.flatten())
                        arcs = self.keyarcs[lpRegion[1],beginSta[ks]:endSta[ks]+1] - self.keyarcs[lpRegion[0],beginSta[ks]:endSta[ks]+1]
                        cur_bom = BOM()
                        cur_bom.layernum = lprow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = self.ispan[beginSta[ks]]
                        cur_bom.endsta = self.ispan[endSta[ks]]
                        cur_bom.maxwidth = np.amax(arcs)
                        cur_bom.avgwidth = np.mean(arcs)
                        cur_bom.area = regionarea
                        cur_bom.thickness = layermult * mat.layerthickness
                        cur_bom.weight = layermult * mat.drydensity * regionarea
                        self.bomIndices['lp'].append([beginSta[ks],endSta[ks],*lpRegion])
                        self.bom['lp'].append(cur_bom)
                        lprow = lprow + 1
                    
                    if comp.group > 0:
                        if swnum != comp.group - 1:
                            swnum = comp.group - 1
                            swrow = 0
                            #NOTE unsure of how to handle swBeginSta -kb
                            swBeginSta.append(beginSta[0])
                            swEndSta.append(endSta[0])
                            self.bom['sw'].append([])
                            self.bomIndices['sw'].append([])
                        # EMA original:
                        # swBeginSta(swnum) = min(beginSta,swBeginSta(swnum));
                        # swEndSta(swnum) = max(endSta,swEndSta(swnum));
                        # changed to:
                        swBeginSta[swnum] = np.amin([*beginSta,swBeginSta[swnum]])
                        swEndSta[swnum] = np.amax([*endSta,swEndSta[swnum]])
                        # END
                        areas = self.webareas[swnum][beginSta[ks]:endSta[ks]]
                        regionarea = sum(areas.flatten())
                        cur_bom = BOM()
                        cur_bom.layernum = swrow
                        cur_bom.materialid = comp.materialid
                        cur_bom.name = comp.name
                        cur_bom.beginsta = self.ispan[beginSta[ks]]
                        cur_bom.endsta = self.ispan[endSta[ks]]
                        cur_bom.maxwidth = np.amax(self.webwidth[swnum])
                        cur_bom.avgwidth = np.mean(self.webwidth[swnum])
                        cur_bom.area = regionarea
                        cur_bom.thickness = layermult * mat.layerthickness
                        cur_bom.weight = layermult * mat.drydensity * regionarea
                        self.bom['sw'][swnum].append(cur_bom)
                        self.bomIndices['sw'][swnum].append([beginSta[ks],endSta[ks]])
                        swrow = swrow + 1
        
        self.bom['lebond'] = sum(self.LEbond) * m_to_mm
        self.bom['tebond'] = sum(self.TEbond) * m_to_mm
        hp_dw = sum([L.weight for L in self.bom['hp']])
        lp_dw = sum([L.weight for L in self.bom['lp']])
        self.bom['dryweight'] = g_to_kg * (hp_dw + lp_dw)
        
        nsw = len(self.bom['sw'])
        self.bom['swbonds'] = [None]*nsw
        for k in range(nsw):
            sw_dw = sum([L.weight for L in self.bom['sw'][k]])
            self.bom['dryweight'] = self.bom['dryweight'] + sw_dw
            C = self.webbonds[k][:,swBeginSta[k]:swEndSta[k]]
            self.bom['swbonds'][k] = m_to_mm * np.sum(C, 1)
        
        # build the material stack for each area
        nSegments = self.keyareas.shape[0]
        nStations = self.keyareas.shape[1]
        nWebs = len(self.bomIndices['sw'])
        segmentLabels = ['HP_TE_FLAT','HP_TE_REINF','HP_TE_PANEL',
                         'HP_SPAR','HP_LE_PANEL','HP_LE','LP_LE',
                         'LP_LE_PANEL','LP_SPAR','LP_TE_PANEL',
                         'LP_TE_REINF','LP_TE_FLAT']
        # self.stacks = [[None]*nStations]*nSegments
        self.stacks = np.empty(shape=(nSegments,nStations), dtype = object)
        for swsk1 in range(nSegments):
            for swsk2 in range(nStations):
                self.stacks[swsk1,swsk2] = Stack()

        for kr in range(nSegments):
            for kc in range(nStations):
                # name the stacks <mm_span_location>_<segmentLabel>
                stack_num = int(np.fix(m_to_mm * self.ispan[kc]))
                self.stacks[kr][kc].name = '{:06d}_{}'.format(stack_num,segmentLabels[kr])
                self.stacks[kr][kc].indices = [kc,kc + 1,kr,kr + 1]
        

        for k in range(len(self.bom['hp'])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = self.bom['hp'][k].name #parent component of ply
            cur_ply.materialid =  self.bom['hp'][k].materialid # materialid of ply
            cur_ply.thickness = self.bom['hp'][k].thickness # thickness [mm] of single ply
            cur_ply.angle = 0, #TODO, set to 0 for now, self.bom['lp'](k, );
            cur_ply.nPlies = 1 # default to 1, modified in addply() if necessary

            # ... and add the ply to every area that is part of the region
            ind = [self.bomIndices['hp'][k][i] for i in range(len(self.bomIndices['hp'][k]))]
            for kr in range(ind[2],ind[3]):
                for kc in range(ind[0],ind[1]):
                    self.stacks[kr][kc].addply(cur_ply)
        
        for k in range(len(self.bom['lp'])):
            # for each row in the BOM, get the ply definition ...
            cur_ply = Ply()
            cur_ply.component = self.bom['lp'][k].name #parent component of ply
            cur_ply.materialid =  self.bom['lp'][k].materialid # materialid of ply
            cur_ply.thickness = self.bom['lp'][k].thickness # thickness [mm] of single ply
            cur_ply.angle = 0, #TODO, set to 0 for now, self.bom['lp'](k, );
            cur_ply.nPlies = 1 # default to 1, modified in addply() if necessary
    
            # ... and add the ply to every area that is part of the region
            ind = [self.bomIndices['lp'][k][i] for i in range(len(self.bomIndices['lp'][k]))]
            for kr in range(ind[2],ind[3]):
                for kc in range(ind[0],ind[1]):
                    self.stacks[kr][kc].addply(cur_ply)
        self.swstacks = [None]*nWebs
        for kw in range(nWebs):
            self.swstacks[kw] = []
            for swsk in range(nStations):
                self.swstacks[kw].append(Stack())
            for kc in range(nStations):
                # name the stacks <mm_span_location>_SW#
                swstack_num = int(np.fix(m_to_mm * self.ispan[kc]))
                self.swstacks[kw][kc].name = '{:06d}_SW{}'.format(swstack_num,kw)
                ind = self.webindices[kw] # currently, the shearweb indices do not change down the span
                self.swstacks[kw][kc].indices = [kc,kc + 1,ind[0],ind[1]]
            for k in range(len(self.bom['sw'][kw])):
                # for each row in the BOM, get the ply definition ...
                cur_ply = Ply()
                cur_ply.component = self.bom['sw'][kw][k].name #parent component of ply
                cur_ply.materialid = self.bom['sw'][kw][k].materialid # materialid of ply
                cur_ply.thickness = self.bom['sw'][kw][k].thickness # thickness [mm] of single ply
                cur_ply.angle = 0, #TODO, set to 0 for now, self.bom['lp'](k, );
                cur_ply.nPlies = 1 # default to 1, modified in addply() if necessary
                # ... and add the ply to every area that is part of the region
                ind = [self.bomIndices['sw'][kw][k][i] for i in range(len(self.bomIndices['sw'][kw][k]))]

                for kc in range(ind[0],ind[1]):
                    self.swstacks[kw][kc].addply(cur_ply)
    # need to add the 'MatDB' information which stores composite stack
    # information in each region at each station
    # see datatypes.MatDBentry
    # prepare material database ==========================================
        self.matdb = []
        #NOTE stopped cleaning here -kb
        for k in range(len(self.materials)):
            cur_entry = MatDBentry()
            cur_entry.name = self.materials[k].name
            cur_entry.type = self.materials[k].type
            cur_entry.ex = self.materials[k].ex
            cur_entry.ey = self.materials[k].ey
            cur_entry.ez = self.materials[k].ez
            cur_entry.gxy = self.materials[k].gxy
            cur_entry.gyz = self.materials[k].gyz
            cur_entry.gxz = self.materials[k].gxz
            if cur_entry.type == 'isotropic':
                cur_entry.nuxy = self.materials[k].prxy
            else:
                cur_entry.prxy = self.materials[k].prxy
                cur_entry.pryz = self.materials[k].pryz
                cur_entry.prxz = self.materials[k].prxz
            cur_entry.dens = self.materials[k].density
            cur_entry.reference = self.materials[k].reference
            self.matdb.append(cur_entry)
        
        flat_stacks = self.stacks.flatten('F')
        for k in range(self.stacks.size):
            cur_entry = MatDBentry()
            cur_entry.name = flat_stacks[k].name
            cur_entry.type = 'composite'
            cur_entry.reference = 'Reference text'
            cur_entry.thicknessType = 'Constant'
            cur_entry.uniqueLayers = len(flat_stacks[k].plygroups)
            cur_entry.symmetryType = 'none'
            cur_entry.layer = [None]*cur_entry.uniqueLayers
            for j in range(cur_entry.uniqueLayers):
                cur_layer = Layer()
                matid = flat_stacks[k].plygroups[j].materialid
                cur_layer.layerName = self.matdb[matid].name
                cur_layer.thicknessA = mm_to_m * flat_stacks[k].plygroups[j].thickness
                cur_layer.thicknessB = cur_layer.thicknessA
                cur_layer.quantity = flat_stacks[k].plygroups[j].nPlies
                cur_layer.theta = flat_stacks[k].plygroups[j].angle
                cur_entry.layer[j] = cur_layer
            self.matdb.append(cur_entry)

        for kw in range(len(self.swstacks)):
            for k in range(len(self.swstacks[kw])):
                cur_entry = MatDBentry()
                cur_entry.name = self.swstacks[kw][k].name
                cur_entry.type = 'composite'
                cur_entry.reference = 'Reference text'
                cur_entry.thicknessType = 'Constant'
                try:
                    cur_entry.uniqueLayers = len(self.swstacks[kw][k].plygroups)
                except TypeError:
                    cur_entry.uniqueLayers = 0
                cur_entry.symmetryType = 'none'
                cur_entry.layer = [None]*cur_entry.uniqueLayers
                for j in range(cur_entry.uniqueLayers):
                    cur_layer = Layer()
                    matid = self.swstacks[kw][k].plygroups[j].materialid
                    cur_layer.layerName = self.matdb[matid].name
                    cur_layer.thicknessA = mm_to_m * self.swstacks[kw][k].plygroups[j].thickness
                    cur_layer.thicknessB = cur_layer.thicknessA
                    cur_layer.quantity = self.swstacks[kw][k].plygroups[j].nPlies
                    cur_layer.theta = self.swstacks[kw][k].plygroups[j].angle
                    cur_entry.layer[j] = cur_layer
                self.matdb.append(cur_entry)
        # shearweb information from NuMAD v1 is formatted in a specific
        # way, recreating that here
        # recreating data.shearweb ====================================
        ctr = 0
        self.shearweb = []
        for kw in range(len(self.swstacks)):
            ind = self.webindices[kw]
            for k in range(len(self.swstacks[kw])):
                if self.swstacks[kw][k].plygroups:
                    cur_sw = Shearweb()
                    cur_sw.Material = self.swstacks[kw][k].name
                    cur_sw.BeginStation = self.swstacks[kw][k].indices[0] # =k
                    cur_sw.EndStation = self.swstacks[kw][k].indices[1] # =k+1
                    cur_sw.Corner = [ind[1]-1,ind[0]-1,ind[0]-1,ind[1]-1]  # dp number is offset by 1 in NuMAD v1
                    self.shearweb.append(cur_sw)
                    ctr += 1
        return


    def updateAirfoilProfile(self,k):
        """
        TODO docstring 
        """
        thickness = self.ithickness[:,k]
        percentthick = self.ipercentthick[k]
        camber = self.icamber[:,k]
        c = self.ic[:,k]
        #jcb: note that I'm using max thickness about camber
        #instead of overall thickness of airfoil. We may need to
        #change this definition.
        maxthick = np.amax(thickness)
        tratio = percentthick / (maxthick * 100)
        thick = thickness * tratio
        hp = camber - 0.5 * thick
        lp = camber + 0.5 * thick
        profile1 = np.concatenate(([c[-1]],np.flipud(c),c[1:],[c[-1]]))
        profile2 = np.concatenate(([0],np.flipud(hp),lp[1:],[0]))
        profile = np.stack((profile1,profile2),axis=1)
        self.profiles[:,:,k] = profile

    
    def updateOMLgeometry(self,k):
        """
        TODO docstring
        """
        x = self.profiles[:,0,k]
        y = self.profiles[:,1,k]
        # self.xoffset[0,k] = c[mtindex]
        if self.naturaloffset:
            x = x - self.xoffset[0,k]
        x = x - self.ichordoffset[k] # apply chordwise offset
        x = x * self.ichord[k] * - 1 * self.rotorspin # scale by chord
        y = y * self.ichord[k] # scale by chord
        twist = - 1 * self.rotorspin * self.idegreestwist[k]
        # prepare for hgtransform rotate & translate
        coords = np.zeros((len(x),4))
        coords[:,0] = np.cos(np.deg2rad(twist)) * x - np.sin(np.deg2rad(twist)) * y
        coords[:,1] = np.sin(np.deg2rad(twist)) * x + np.cos(np.deg2rad(twist)) * y
        coords[:,2] = np.zeros(len(x))
        coords[:,3] = np.ones(len(x))
        # use the generating line to translate and rotate the coordinates
        # NOTE currently, rotation is not assigned from blade properties
        # and defaults to 0
        prebend_rot = 0
        sweep_rot = 0
        """
        jcb: This code, copied from NuMAD 2.0, causes each section to rotate out
        of plane so that its normal follows the generating line direction. Need
        to replace 'twistFlag' with '-1*self.rotorspin' and calculate the slopes
        based on the available data. For now, default to parallel sections.
                        if isequal(blade.PresweepRef.method,'normal')
                            sweep_slope = ppval(blade.PresweepRef.dpp,sta.LocationZ);
                            sweep_rot = atan(sweep_slope*twistFlag);
                        end
                        if isequal(blade.PrecurveRef.method,'normal')
                            prebend_slope = ppval(blade.PrecurveRef.dpp,sta.LocationZ);
                            prebend_rot = atan(-prebend_slope);
                        endc
        """
        transX = - 1 * self.rotorspin * self.isweep[k]
        transY = self.iprebend[k]
        transZ = self.ispan[k]
        Ry = rotation('y', sweep_rot)
        Rx = rotation('x', prebend_rot)
        R = Ry@Rx
        T = translation(transX,transY,transZ)
        coords = coords @ np.transpose(R) @ np.transpose(T)
        # save the transformed coordinates
        self.geometry[:,:,k] = coords[:,0:3]
        # self.geometry[:,0,k] = coords[:,0]
        # self.geometry[:,1,k] = coords[:,1]
        # self.geometry[:,2,k] = coords[:,2]
        pass


    def addStation(self,af = None,spanlocation = None): 
        """This method adds a station

        Specifically, the station object is created 
        and appended to self.stations.

        Parameters
        ----------
        af : airfoil
        spanlocation : float

        Returns
        -------
        None

        Example
        -------
        ``blade.addStation(af,spanlocation)`` where  ``af`` = airfoil filename 
        or ``AirfoilDef`` object
        """
        newStation = Station(af)
        newStation.spanlocation = spanlocation
        newStation.parent = self
        if self.stations:
            self.stations.append(newStation)
        else:
            self.stations = []
            self.stations.append(newStation)
        
        # N = np.asarray(self.stations).size
        # k = N + 1
        # if k > 1:
        #     self.stations[k] = StationDef(af)
        # else:
        #     self.stations = StationDef(af)
        
        # self.stations[k].spanlocation = spanlocation
        # self.stations[k].parent = self
        return

        
    def findLayerExtents(self,layerDist = None,layerN = None): 
        """
        TODO docstring
        """
        assert np.isscalar(layerN),'second argument "layerN" must be a scalar'
        staLogical = layerDist >= layerN
        prev = 0
        beginSta = []
        endSta = []
        for k in range(len(staLogical)):
            if staLogical[k] == 1 and prev == 0:
                beginSta.append(k)
            if staLogical[k] == 0 and prev == 1:
                endSta.append(k)
            elif k == len(staLogical)-1 and prev == 1:
                    endSta.append(k)
            prev = staLogical[k]
        
        return beginSta,endSta
        

    def findRegionExtents(self,comp = None): 
        """
        TODO docstring
        """
        le = self.keylabels.index('le')
        # "keylabels" is expected to wrap from te on hp side around to te on lp side
        try:
            if len(comp.hpextents) == 2:
                try:
                    hp1 = self.keylabels[0:le+1].index(comp.hpextents[0])
                except KeyError:
                    print(f'HP extent label "{comp.hpextents[0]}" not defined.')
                try:
                    hp2 = self.keylabels[0:le+1].index(comp.hpextents[1])
                except KeyError:
                    print(f'HP extent label "{comp.hpextents[1]}" not defined.')
                hpRegion = [hp1,hp2]
                hpRegion.sort()
            else:
                hpRegion = []
        except TypeError:
            hpRegion = []
        
        try:
            if len(comp.lpextents) == 2:
                try:
                    lp1 = self.keylabels[le:].index(comp.lpextents[0]) + le 
                except KeyError:
                    print(f'HP extent label "{comp.hpextents[0]}" not defined.')
                try:
                    lp2 = self.keylabels[le:].index(comp.lpextents[1]) + le
                except KeyError:
                    print(f'HP extent label "{comp.hpextents[1]}" not defined.')
                lpRegion = [lp1,lp2]
                lpRegion.sort()
            else:
                lpRegion = []
        except:
            lpRegion = []
        
        # if length(comp['hpextents'])==1 && length(comp['lpextents'])==1
            # sw1 = find(1==strcmpi(comp['hpextents']{1},keylabels(1:le)));
            # assert(~isempty(sw1),'HP extent label "#s" not defined.',comp['hpextents']{1});
            #  w2 = find(1==strcmpi(comp['lpextents']{1},keylabels(le:end))) + le-1;
            # assert(~isempty(sw2),'LP extent label "#s" not defined.',comp['lpextents']{1});
            # swRegion = [sw1 sw2];
        # else
        swRegion = []
        return hpRegion,lpRegion #,swRegion
        

    def getTEtype(self,xy):
        """
        TODO docstring

        Parameters
        ----------
        xy : ndarray

        Returns
        -------

        tetype : str

        """
        unitNormals = getAirfoilNormals(xy)
        angleChange = getAirfoilNormalsAngleChange(unitNormals)
        disconts = np.flatnonzero(angleChange>30)

        if np.std(angleChange) < 2:
            tetype = 'round'
        elif len(disconts) > 1:
            tetype = 'flat'
        else:
            tetype = 'sharp'
        return tetype


    def getprofileTEtype(self,k: int):
        """
        TODO docstring
        """
        xy = self.profiles[:,:,k]
        tetype = self.getTEtype(xy)       
        return tetype
    
    
    def expandBladeGeometryTEs(self): 
        """NOTE: is this deprecated? I don't see it used anywhere
        TODO: docstring
        """
        nStations = self.geometry.shape[2]
        minimumTEedgelength = 0.003

        for iStation in range(0,nStations):
            firstPoint = self.ichord[iStation] * self.profiles[-2,:,iStation]
            secondPont = self.ichord[iStation] * self.profiles[1,:,iStation]
            edgeLength = np.linalg.norm(secondPont - firstPoint)
            #fprintf('station #i, edgeLength: #f\n',iStation,edgeLength*1000)

            maxthick = np.amax(self.ithickness[:,iStation])
            mtindex = np.argmax(self.ithickness[:,iStation])
            tratio = self.ipercentthick[iStation] / (maxthick * 100)
            airFoilThickness = self.ithickness[:,iStation] * tratio
            onset = self.ic[mtindex,iStation]
            if edgeLength < minimumTEedgelength:
                tet = (minimumTEedgelength - edgeLength) / self.ichord[iStation]
                tes = 5 / 3 * tet # slope of TE adjustment; 5/3*tet is "natural"
                # continuous first & second derivatives at 'onset'
                # maintain second & third derivative at mc==1 (TE)
                # adjust slope at mc==1 (TE) by tes
                A = np.array([
                    [1,1,1,1],
                    [3,4,5,6],
                    [6,12,20,30],
                    [6,24,60,120]
                    ])
                d = np.array([[tet],[tes],[0],[0]])
                p = np.linalg.solve(A,d)
                #onset = self(k).maxthick;  # start of TE modification, measured from LE
                mc = np.amax((self.ic[:,iStation] - onset) / (1 - onset),0)
                temod = np.array([mc ** 3,mc ** 4,mc ** 5,mc ** 6]) * p
                airFoilThickness = airFoilThickness + temod

                self.ithickness[:,iStation] = airFoilThickness / tratio
                self.updateAirfoilProfile(iStation)

                mtindex = np.argmax(self.ithickness[:,iStation])
                self.xoffset[0,iStation] = self.ic[mtindex,iStation]
                self.updateOMLgeometry(iStation)
                # firstPoint=self.ichord(iStation)*self.profiles(end-1,:,iStation);
                # secondPont=self.ichord(iStation)*self.profiles(2,:,iStation);
                # edgeLength2=norm(secondPont-firstPoint);
                # fprintf('station #i, edgeLength: #f, New edgeLength=#f, percent diff: #f\n',iStation,edgeLength*1000,edgeLength2*1000,(edgeLength2-edgeLength)/edgeLength2*100)
        
        self.updateKeypoints


    def getShellMesh(self, includeAdhesive): 
        nodes,elements,outerShellElSets,shearWebElSets,adhesNds,adhesEls = shellMeshGeneral(self,0,includeAdhesive)
        return nodes,elements,outerShellElSets,shearWebElSets,adhesNds,adhesEls
    
    """
    #NOTE need team help here -kb
    # not converted
    def generateBeamModel(self): 
        #This method generates blade sectional properties used for
        #aeroelastic analyses
        #TODO extend docstring
        
        
        global precompPath
        global bmodesPath
        parID = gcp('nocreate')
        if not len(parID)==0 :
            batchRun = True
        else:
            batchRun = False
        
        # NOTES: ******************************************************
        # 1. FIXED -- needs to read MatDBsi.txt file, store this
        # internally (blade.matdb)
        # 2. FIX -- the code creates files prepmat -- why? is this
        # needed elsewhere?
        # 3. using blade.profiles instead of data.station(ii).coords
        # 4. FIX -- web data not saved in same format as NuMAD
        # (data.shearweb structure in NuMAD)
        # *************************************************************
            
        # ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        # code from BladeDef_to_NuMADfile and NuMAD2PreComp files
        # need to add material properties for end station
        if self.stacks.shape[2-1] < len(self.ispan):
            for ii in np.arange(1,self.stacks.shape[1-1]+1).reshape(-1):
                self.stacks(ii,len(self.ispan)).name = '**UNSPECIFIED**'
        
        # original code from NuMAD files ended ------------------------
        # ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        
        # create precomp input files
        precomp = Blade2PreComp(self,self.matdb)
        # run precomp
        PreComp_SectionData = runPreCompAnalysis(precomp,precompPath,batchRun)
        # Use BModes to calculate mode shapes based on section properties from PreComp
        bmodesFrequencies = BModes4Blade2PreComp2FASTBlade(self,bmodesPath,PreComp_SectionData.data)
        # fit polynomials to the modes
        if batchRun:
            modeShapes = polyfitmodes(np.array([1,3,2]))
        else:
            modeShapes = polyfitmodes
        
        # Generate FAST Blade file from this analysis
        warnings.warn('need PresweepRef and PrecurveRef variables from NuMAD')
        PreComp2FASTBlade_BladeDef(self,PreComp_SectionData.data,modeShapes,batchRun)
        print('FAST Blade file has been written: FASTBlade_precomp.dat')
        # Input file cleanup
        if batchRun:
            qu = 'Yes'
        else:
            qu = questdlg('Delete all miscellaneous input/output files?','File Cleanup...','Yes','No','No')
        
        # switch qu
    #     case 'Yes'
    #         delete('*.pci')
    #         delete('*.inp')
    #         delete PrepMat.txt
    #         delete bmodes*.bmi
    #         delete blade_sec_props.dat
    #         delete bmodes*.echo
    #     case 'No'
    # end
    # above lines caused problems with ML2PY
        
        # Estimate the flexural bending distance
        warnings.warn('this code is not correct -- should be based on distance to neutral axis, not purely geometry')
        make_c_array_BladeDef(self)
        return bmodesFrequencies

    def old_getTEtype(self,xy = None):
        
        if np.abs(xy[1,1] - xy[-2,1]) > 1e-5:
            # y-diff of second and end-1 points is non-zero for flatback
            tetype = 'flat'
            print('FLATBACK AIRFOIL')
        else:
            # y-diff of first two points is zero otherwise (point
            # is duplicated)
            hp_angle = np.arctan2(xy[1,1] - xy[2,1],xy[1,0] - xy[2,0])
            lp_angle = np.arctan2(xy[-3,1] - xy[-2,1],xy[-2,0] - xy[-3,0])
            if (hp_angle + lp_angle) > 0.8 * np.pi:
                # if angle is approaching 180deg, then treat as
                # 'round'
                # jcb: it may be better to base this decision on
                # continuity of slope or curvature
                tetype = 'round'
            else:
                tetype = 'sharp'
        return tetype

    def generateFEA(self): 
        # This method generates FEA
        
        # 1. check that functionality from original code is not needed
        # 2. FIXED -- shell7 needs to read MatDBsi.txt file, store this
        # internally if the file isn't used elsewhere (blade.matdb)
        # 3. FIX -- in shell7 file, need to define PrecurveRef and
        # PresweepRef (currently set to zero)
        # 4. FIX -- web data not saved in same format as NuMAD
        # (data.shearweb structure in NuMAD)
        
        # NOTE:can add flags into the call -- e.g., element type, ...
        global ansysPath
        # define ANSYS model settings (can be options in generateFEA)
        config.ansys.BoundaryCondition = 'cantilered'
        config.ansys.ElementSystem = '181'
        config.ansys.MultipleLayerBehavior = 'distinct'
        config.ansys.meshing = 'elementsize'
        config.ansys.smartsize = 5
        config.ansys.elementsize = 0.2026
        config.ansys.shell7gen = 1
        config.ansys.dbgen = 1
        fcopts = np.array(['EMAX','SMAX','TWSI','TWSR','HFIB','HMAT','PFIB','PMAT','L3FB','L3MT','L4FB','L4MT','USR1','USR2','USR3','USR4','USR5','USR6','USR7','USR8','USR9'])
        config.ansys.FailureCriteria = cell(np.asarray(fcopts).size,2)
        config.ansys.FailureCriteria[:,1] = np.transpose(fcopts)
        config.ansys.FailureCriteria[:,2] = deal(np.array([False]))
        # Generate a mesh using shell elements -- can add options here
        shell7_name = 'shell7.src'
        ansys_product = 'ANSYS'
        self.paths.job = pwd
        
        filename = fullfile(self.paths.job,shell7_name)
        develop__write_shell7(self,filename)
        if self.ansys.dbgen:
            if len(ansysPath)==0:
                errordlg('Path to ANSYS not specified. Aborting.','Operation Not Permitted')
                return
            try:
                #tcl: exec "$ANSYS_path" -b -p $AnsysProductVariable -I shell7.src -o output.txt
                ansys_call = sprintf('"#s" -b -p #s -I #s -o output.txt',ansysPath,ansys_product,shell7_name)
                status,result = dos(ansys_call)
                if status==0:
                    # dos command completed successfully; log written to output.txt
                    if 1:
                        print('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.')
                    else:
                        helpdlg('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.','ANSYS Call Completed')
                if status==7:
                    # an error has occured which is stored in output.txt
                    if 1:
                        print('Could not complete ANSYS call. See "output.txt" for details.')
                    else:
                        warndlg('Could not complete ANSYS call. See "output.txt" for details.','Error: ANSYS Call')
            finally:
                pass
        
        return
        
    #NOTE ask team about this
    # ignore for now...
    def writeBOMxls(self,file = None): 
        # This method writes the bill-of-materials out to a spreadsheet.
        
        # Example:
        
        #   ``bladeDef.writeBOMxls('bom.xlsx')``
        
        m_to_mm = 1000.0
        if os.path.exist(str('BOM_template.xlsx')):
            copyfile('BOM_template.xlsx',file)
        
        header = np.array([['Layer #','Material ID','Component','Begin Station','End Station','Max width','Mean width','3D area','Layer Thickness','Computed layer weight'],['','','','(m)','(m)','(m)','(m)','(m^2)','(mm)','(g)']])
        # LP skin table
        array_ = np.array([[header],[self.bom['lp']]])
        xlswrite(file,array_,'LP skin')
        # HP skin table
        array_ = np.array([[header],[self.bom['hp']]])
        xlswrite(file,array_,'HP skin')
        # shear web table
        array_ = np.array([np.array([['SW #'],['']]),header])
        
        for k in np.arange(1,np.asarray(self.bom.sw).size+1).reshape(-1):
            nr = self.bom.sw[k].shape[1-1]
            array_ = np.array([[array_],[np.array([np.matlib.repmat(np.array([k]),nr,1),self.bom.sw[k]])]])
        
        xlswrite(file,array_,'shear webs')
        # bond line lengths
        array_ = np.array([['','Length'],['','(mm)'],['Root diameter',self.ichord(1) * m_to_mm],['LE bond',np.round(self.bom.lebond)],['TE bond',np.round(self.bom.tebond)]])
        for r in np.arange(1,2+1).reshape(-1):
            for k in np.arange(1,np.asarray(self.bom.swbonds).size+1).reshape(-1):
                surfs = np.array(['HP','LP'])
                str = sprintf('#s bond, SW #d',surfs[r],k)
                cellrow = np.array([str,np.round(self.bom.swbonds[k](r))])
                array_ = np.array([[array_],[cellrow]])
        
        xlswrite(file,array_,'lengths')
        return
        
     
    def writePlot3D(self,file = None,breakpoints = None):
        #NOTE ask team about this
        # ignore for now  

        # Write the current blade geometry in Plot3D format.
        # breakpoints is a list of chord fractions at which the
        # surface geometry is divided into blocks
        
        # Examples:
        
        #   ``BladeDef.writePlot3D(filename,[breakpoints])``
        
        #   ``BladeDef.writePlot3D('file.p3d',[-.3, .3]);``
        
        if not ('breakpoints' is not None) :
            breakpoints = []
        
        indicesOfBreakpoints = np.zeros((1,np.asarray(breakpoints).size))
        # get the chordwise spacing of points, assuming identical
        # spacing for all stations
        chordspacing = self.cpos[:,1]
        for kBreakpoint in np.arange(1,np.asarray(breakpoints).size+1).reshape(-1):
            bp = breakpoints(kBreakpoint)
            __,ind = np.amin(np.sqrt((chordspacing - bp) ** 2))
            indicesOfBreakpoints[kBreakpoint] = ind
        
        N,M = self.cpos.shape
        INCLUDE_REPEATS = False
        if INCLUDE_REPEATS:
            indicesOfBreakpoints = unique(np.array([1,indicesOfBreakpoints,N]))
        else:
            indicesOfBreakpoints = unique(np.array([2,indicesOfBreakpoints,N - 1]))
        
        fid = open(file,'wt')
        
        if (fid == - 1):
            raise Exception('Could not open file "#s"',file)
        
        # output the data in Plot3d format
        #TODO
        # try:
        #     nBlocks = np.asarray(indicesOfBreakpoints).size - 1
        #     fid.write('#d\n' # (nBlocks))
        #     for kblock in np.arange(1,nBlocks+1).reshape(-1):
        #         a = indicesOfBreakpoints(kblock)
        #         b = indicesOfBreakpoints(kblock + 1)
        #         fid.write('#d  #d  #d\n' # (1 + b - a,M,1))
        #     columnsPerLine = 5
        #     for kblock in np.arange(1,nBlocks+1).reshape(-1):
        #         a = indicesOfBreakpoints(kblock)
        #         b = indicesOfBreakpoints(kblock + 1)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),1,:),columnsPerLine)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),2,:),columnsPerLine)
        #         self.fprintf_matrix(fid,self.geometry(np.arange(a,b+1),3,:),columnsPerLine)
        # finally:
        #     pass
        
        # fid.close()
        # return
    
    def fprintf_matrix(self,fid = None,matrixData = None,columnsPerLine = None): 
        kColumn = 1
        for kData in np.arange(1,np.asarray(matrixData).size+1).reshape(-1):
            fid.write('%g ', (matrixData(kData)))
            kColumn = kColumn + 1
            if kColumn > columnsPerLine:
                fid.write('\n')
                kColumn = 1
        
        return
    def downsampleProfile(self,k = None,n_points = None): 
        #NOTE doesn't seem to be relevant
        # ignoring for now -kb
        # ble: can this be deleted?? Doesn't precisely control the
        # number of points. replace with resampleAirfoil_ble
        # currently used by BladeDef_to_NuMADfile and develop__write_shell7
        N = self.profiles.shape[1-1]
        assert (np.isscalar(k),'Profile index "k" must be scalar')
        assert (n_points < N,'n_points must be less than length of profile')
        dk = np.round((N - 2) / n_points)
        LE = self.LEindex
        tetype = self.getprofileTEtype(k)
        if 'flat' == tetype[0]:
            ind = unique(np.array([1,np.arange(2,LE+dk,dk),LE,fliplr(np.arange(N - 1,LE+- dk,- dk))]))
        else:
            if np.array(['sharp','round']) == tetype[0]:
                ind = unique(np.array([np.arange(2,LE+dk,dk),LE,fliplr(np.arange(N - 2,LE+- dk,- dk))]))
        
        coords = self.profiles[ind,:,k]
        return coords

        
    #     function delete(self)
    #     try ##ok<TRYNC>
    #         delete(self.bomPlot.hgLinesHP);
    #         delete(self.bomPlot.hgLinesLP);
    #         delete(self.bomPlot.hgPatchHP);
    #         delete(self.bomPlot.hgPatchLP);
    #         delete(self.bomPlot.uisliderHP);
    # #                 delete(self.bomPlot.uisliderLP);
    #         delete(self.bomPlot.hTitleHP);
    #         delete(self.bomPlot.hTitleLP);
    #     end
    # end
        
    
    def surf(self):
        #NOTE not used
        #ignoring for now -kb
        h = surf(np.squeeze(self.geometry[:,3,:]),np.squeeze(self.geometry[:,1,:]),np.squeeze(self.geometry[:,2,:]),'MeshStyle','column')
        if nargout > 0:
            varargout = np.array([h])
        
        return varargout
        

    # Ignoring all plot stuff for now  
    def plotregions(self): 
        # try ##ok<TRYNC>
    #     delete(self.hgKeypoints)
    # end
    # COMMENTED FOR ML2PY PROBLEMS
        M = self.keypoints.shape[1-1]
        for km in np.arange(1,M+1).reshape(-1):
            self.hgKeypoints[km] = line(np.squeeze(self.keypoints[km,3,:]),np.squeeze(self.keypoints[km,1,:]),np.squeeze(self.keypoints[km,2,:]))
        
        return
        
        # 
    def plotgeom(self): 
        # try ##ok<TRYNC>
    #     delete(self.hgGeometry)
    # COMMENTED FOR ML2PY PROBLEMS
    # end
        N = self.geometry.shape[3-1]
        for k in np.arange(1,N+1).reshape(-1):
            self.hgGeometry[k] = line(self.geometry[:,3,k],self.geometry[:,1,k],self.geometry[:,2,k])
        # 
        return
        # 
        # 
    def plotbom(self,k = None): 
        if not ('k' is not None)  or len(k)==0:
            k = 1
        # 
        if k=='cb':
            k = np.rint(get(self.bomPlot.uisliderHP,'Value'))
        # 
        if len(self.bomPlot.hgLinesHP)==0 or not np.all(ishandle(self.bomPlot.hgLinesHP)) :
            clf
            self.bomPlot.axHP = axes('Position',np.array([0.1,0.6,0.8,0.3]))
            self.bomPlot.hgLinesHP = plt.plot(self.ispan,self.keyarcs(1,:) - self.HParcx0,'k-.',self.ispan,self.keyarcs(7,:) - self.HParcx0,'k-.')
            self.bomPlot.hTitleHP = plt.title('','Interpreter','none')
            self.bomPlot.axLP = axes('Position',np.array([0.1,0.2,0.8,0.3]))
            self.bomPlot.hgLinesLP = plt.plot(self.ispan,self.keyarcs(13,:) - self.LParcx0,'k-.',self.ispan,self.keyarcs(7,:) - self.LParcx0,'k-.')
            self.bomPlot.hTitleLP = plt.title('','Interpreter','none')
            n = self.bom['hp'].shape[1-1]
            self.bomPlot.uisliderHP = uicontrol('Style','slider','Min',1,'Max',n,'Value',1,'SliderStep',np.array([1 / (n - 1),10 / (n - 1)]),'Units','normalized','Position',np.array([0.02,0.02,0.96,0.04]),'Callback',lambda src = None,evt = None: self.plotbom('cb'))
            #             self.bomPlot.uisliderLP = uicontrol('Style','slider',...
    #                 'Min',1,'Max',n,'Value',1,'SliderStep',[1/(n-1) 10/(n-1)],...
    #                 'Units','normalized','Position',[.02 .02 .96 .04],...
    #                 'Callback',@(src,evt) self.plotbom('cb'));
        # 
        k = np.amin(k,self.bom['hp'].shape[1-1])
        self.bomPlot.kLayer = k
        set(self.bomPlot.uisliderHP,'Value',k)
        str = sprintf('HP Layer #d: #s',k,self.bom['hp'][k,3])
        set(self.bomPlot.hTitleHP,'String',str)
        str = sprintf('LP Layer #d: #s',k,self.bom['lp'][k,3])
        set(self.bomPlot.hTitleLP,'String',str)
        hp = self.bomIndices.hp(k,:)
        x = self.ispan(np.arange(hp(1),hp(2)+1))
        y1 = self.keyarcs(hp(3),np.arange(hp(1),hp(2)+1)) - self.HParcx0(np.arange(hp(1),hp(2)+1))
        y2 = self.keyarcs(hp(4),np.arange(hp(1),hp(2)+1)) - self.HParcx0(np.arange(hp(1),hp(2)+1))
        if ishandle(self.bomPlot.hgPatchHP):
            os.delete(self.bomPlot.hgPatchHP)
        # 
        axes(self.bomPlot.axHP)
        self.bomPlot.hgPatchHP = patch(np.array([x,fliplr(x)]),np.array([y1,fliplr(y2)]),'b')
        lp = self.bomIndices.lp(k,:)
        x = self.ispan(np.arange(lp(1),lp(2)+1))
        y1 = self.keyarcs(lp(3),np.arange(lp(1),lp(2)+1)) - self.LParcx0(np.arange(lp(1),lp(2)+1))
        y2 = self.keyarcs(lp(4),np.arange(lp(1),lp(2)+1)) - self.LParcx0(np.arange(lp(1),lp(2)+1))
        if ishandle(self.bomPlot.hgPatchLP):
            os.delete(self.bomPlot.hgPatchLP)
        # 
        axes(self.bomPlot.axLP)
        self.bomPlot.hgPatchLP = patch(np.array([x,fliplr(x)]),np.array([y1,fliplr(y2)]),'b')
        return
        
        
    def plotprofile(self,k = None): 
        # This method plots profiles
        
        # Examples:
        
        #   ``blade.plotprofile(1);``
        
        #   ``blade.plotprofile(1:N);``
        
        # plt.plot(np.squeeze(self.profiles(:,1,k)),np.squeeze(self.profiles(:,2,k)),'.-')
        return
    
    """
