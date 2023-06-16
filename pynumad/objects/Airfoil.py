########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import re
import numpy as np

import os
import matplotlib.pyplot as plt
import warnings
import scipy as sp
from typing import Optional

from pynumad.io.xml_to_airfoil import xml_to_airfoil
from pynumad.utils.interpolation import interpolator_wrap

from numpy import ndarray
    
class Airfoil():
    """Airfoil object

    Parameters
    ----------
    coords : array
    ref : string
        Name of airfoil reference
    filename : string
        path to filename of airfoil (xml somewhat supported)

    Attributes
    ----------
    name : str
        User selected airfoil name
    reference : string
        Header info in file
    coordinates : array
        Profile data in two columns
    c : array
        Computed by NuMAD
    camber : array
        Camber line as a function of x. 
        Distance in percent chord between LP and HP curves. 
        Computed by NuAMD.
    thickness : float
        Relative thickness as a function of the x coordinate. 
        Values between 0 and 1, where 1 corresponds to maximum thickness. 
        Computed by NuAMD.
    percentthick : float
        Maximum airfoil thickness as a percentage of chord length [#]
    maxthick : float
        Airfoil thickness as a percentage of chord. Values between 0 and 1.
    TEtype : str
        Options, 'round', 'sharp', or 'flat'
    """
    def __init__(self, filename = None, coords = None, ref = None):
        self.name : str = None
        self.reference : str = None
        self.coordinates : ndarray = None
        self.c : ndarray = None
        self.camber : ndarray = None 
        self.thickness : float = None
        self.percentthick : float = None
        self.maxthick : float = None
        self.TEtype : str = None

        if filename:
            # currently assuming XML format
            # Use the base filename as the airfoil name
            #TODO fix finding fn
            __,fn,__ = os.path.split(filename)[0],os.path.splitext(os.path.split(filename)[1])[0],os.path.splitext(os.path.split(filename)[1])[1]
            self.name = fn
            # Open the file and read the entire contents
            with open(filename,'r') as f:
                file_contents = f.read().splitlines()
            self.read_xml(file_contents)

        #   TODO decide if AirfoilColumns is used or if we deprecate it
        #   coords,ref = readAirfoilColumns(fid)
        #   self.reference = ref
        #   self.coordinates = coords
        else:
            try:
                # check if ref is a string and coords an array
                ref*5
                coords*5
                # if so, assign as attributes
                self.name = ref
                self.reference = ref
                self.coordinates = coords
            except TypeError:
                # otherwise, provide default AF profile
                self.name = 'circular'
                self.reference = ''
                theta = np.linspace(0,np.pi,50)
                theta = np.concatenate((theta,[np.pi],theta[1:-1]))
                xcoord = 0.5 * np.cos(-1 * theta) + 0.5
                ycoord = 0.5 * np.sin(-1 * theta)
                self.coordinates = np.stack((xcoord,ycoord),axis=1)
        self.manageTE()


    # Properties
    @property
    def x(self):
        """Horizontal axis of Airfoil shape coordinates Working 
        clockwise starting from the TE to the LE and back to the TE. 
        LE must be at (1,0) and TE at (0,0). 
        Needed only by ``AirfoilDef.plot``

        TODO docstring
        """
        cc = self.c
        xcoord = np.concatenate([[cc[-1]], np.flipud(cc), cc[1:], [cc[-1]]])
        return xcoord

    
    #TODO check this func
    @property
    def y(self):
        """Vertical axis of Airfoil shape coordinates 
        Working clockwise starting from the TE to the LE and back to the TE. 
        LE must be at (1,0) and TE at (0,0). 
        Needed only by ``AirfoilDef.plot``

        TODO docstring
        """
        lp = self.camber + (self.thickness / 2)
        hp = self.camber - (self.thickness / 2)
        ycoord = np.concatenate(([0], np.flipud(hp), lp[1:], [0]))
        return ycoord


    ### IO

    def read_xml(self, filename):
        """
        TODO docstring
        """
        xml_to_airfoil(self,filename)
        return self


    def manageTE(self) -> None:
        """TODO docstring
        
        Parameters
        ----------

        Returns
        -------
        """

        unitNormals=getAirfoilNormals(self.coordinates)
        angleChange=getAirfoilNormalsAngleChange(unitNormals)
        discontinuities = np.flatnonzero(angleChange>45)
                    
        if discontinuities.shape[0] == 2:
            #Flatback piece in airfoil. delete for resampling
            if (discontinuities[0] == 0) & (discontinuities[1] == angleChange.shape[0] - 1):
                indexToDelete = 0
            else:
                indexToDelete = angleChange.shape[0] - 1
            newcoords = np.delete(self.coordinates,indexToDelete,0)
            self.coordinates = newcoords
        
        # ddof set to 1 to match default matlab behavior
        if np.std(angleChange,ddof=1) < 1:
            self.TEtype='round'
        return self


    ### Geometry
    
    def resample(self,n_samples: int = 150,spacing: str = 'cosine') -> None: 
        """TODO docstring
        
        Parameters
        ----------

        Returns
        -------

        Example
        -------
        AirfoilDef.resample
        af.resample(n_samples,spacing)
            spacing = 'auto' | 'half-cosine' | 'cosine'
        
        af.resample(200,'half-cosine');
        """
        af_out = resampleAirfoil(self.coordinates,n_samples,spacing)
        xcoord = af_out[:,0]
        ycoord = af_out[:,1]
        # self(k).percentthick = (max(ycoord) - min(ycoord))*100;
        self.c,self.camber,self.thickness = computeCamberAndThickness(xcoord,ycoord)
        m = np.max(self.thickness)
        i = np.argmax(self.thickness)
        self.percentthick = m * 100
        self.maxthick = self.c[i]
        if not self.TEtype or ('round' not in self.TEtype):
            if np.abs(self.thickness[-1]) < 1e-4:
                self.TEtype = 'sharp'
            else:
                self.TEtype = 'flatback'
        return self
        

    #currently unused
    def adjustTE(self, tet, tes, onset) -> None:
        """TODO docstring
        
        Parameters
        ----------

        Returns
        -------

        Example
        -------
        AirfoilDef.adjustTE
        af.adjustTE(TE_thick,[TE_slope],[onset])
        TE_thick is the amount of TE thickness to add
        TE_slope is the slope of the added thickness profile at TE,
                    defaults to 5/3 * TE_thick
        onset    is the chord fraction where adjustment begins,
                    defaults to location of max thickness
        
         af.adjustTE(0.02);
         af.adjustTE(0.02,0);
         af.adjustTE(0.02,[],0.8);
        """

        if not tes:
            tes = 5 / 3 * tet # slope of TE adjustment; 5/3*tet is "natural"
        
        if not onset:
            USEMAXTHICK = True
        else:
            USEMAXTHICK = False # use the given 'onset' instead
        # continuous first & second derivatives at 'onset'
        # maintain second & third derivative at mc==1 (TE)
        # adjust slope at mc==1 (TE) by tes
        A = np.array([[1,1,1,1],[3,4,5,6],[6,12,20,30],[6,24,60,120]])
        d = np.array([[tet],[tes],[0],[0]])
        p = np.linalg.solve(A,d)
        if USEMAXTHICK:
            onset = self.maxthick
        mc = np.amax((self.c - onset) / (1 - onset),0)
        temod = np.array([mc ** 3,mc ** 4,mc ** 5,mc ** 6]) * p
        self.thickness = self.thickness + temod
        pass


    ### Plotting 

    def plotAirfoil(self): 
        """ Plot airfoil
        """
        fig, ax = plt.subplots()
        # ax[0].plot(self.x,self.y,'.-')
        ax.plot(self.coordinates[:,0],self.coordinates[:,1],'.-')
        ax.plot(self.c,self.camber)
        # mtx = self.maxthick * np.array([1,1])
        # kn = find(self.c >= self.maxthick,1)
        # mty = self.camber(kn) + self.thickness(kn) * np.array([0.5,- 0.5])
        # line(mtx,mty,'LineStyle',':','Color','k')
        # else:
        fig.show()
        return fig, ax


### Helper functions

def resampleAirfoil(af_in, n_samples, spacing): 
    """Resample airfoil coordinates

    Parameters
    ----------
    af_in : array
        Nx2 array containing N normalized xy airfoil points
    n_samples : int
        number of points to be created around surface
    spacing : string
        spacing routine to be used: 
        'cosine', 'half-cosine', 'constant', 'auto'
    
    Returns
    -------
    af_out : array 
        array containing n_samples+1 airfoil points
       
    Cosine spacing: puts higher density of points at both LE and TE;
        constant arc length point spacing around a perfect circle.
    Half-cosine spacing: puts higher density of points at LE and lesser
        density of points at TE
    constant spacing: constant spacing of points along chord line
    auto: choose between Cosine and Half-cosine based on steepness of TE
       
    Assumes coordinates begin at trailing edge (x,y)=(1,0) and trace out the HP
    surface, then the LP surface.
    
    Flatback airfoil inputs are designated by ensuring the following:
        Point #1 = (1,0)   (center of trailing edge)
        Point #2 = (1,-y) where y~=0  (HP corner of flatback trailing edge)
        Point #N = (1,y)  where y~=0  (LP corner of flatback trailing edge)
    Notes:

    * This routine enforces LE point is at (0,0) - Warning, this may have
      complications when it comes ot CFD analyses!!
    * spline_type = spline algorithm to be used for oversampling:
      'linear', 'pchip', 'spline'; right now, the function is hard-coded
      for 'spline', but others can be used by changing the te_type setting
      in the code.

    Assumes leading edge at (0,0) and trailing edge at (1,0)

    Example
    -------
    af_out = resampleAirfoil_nmd(af_in, n_samples, spacing)

    JP: initial creation
    BRR: modified for release 1/25/2011
    """

    # Error checking
    # if airfoil coordinates are not in Nx2 array
    if af_in.shape[1] != 2:
        tmpN = af_in.shape[0]
        tmpM = af_in.shape[1]
        warnings.warn('af_in array was defined in '+ str(tmpN)+'x'+str(tmpM)+' array.  Automatically changing it to an '+str(tmpM)+'x'+str(tmpN)+' array.')
        af_in = np.transpose(af_in)
    
    # End error checking routines
    xy = af_in
    
    #Calculate arc length of xy points clockwise from trailing edge
    n_points = xy.shape[0]
    t = np.zeros(n_points)
    for i in range(1,n_points):
        #formula: t(i) = hypot( x(i)-x(i-1) , y(i)-y(i-1) ) + t(i-1);
        t[i] = np.hypot(xy[i,0] - xy[i-1,0],xy[i,1] - xy[i-1,1]) + t[i-1]
    
    #Get total arc length
    arc_length = t[-1]

    #Spline airfoil with many points
    oversample = 10000
    delta = arc_length / (oversample - 1)

    #  The manypoints range from 0 to total arc_length, adding a bit on each
    #  side so that flatbacks extend past x=1 after rotation corrections.
    manypoints = np.linspace(- delta,arc_length + delta,num = oversample+2)
    spline_type = 'pchip'
    if (np.array(['linear','pchip','spline']) == spline_type).any():
        xxyy = interpolator_wrap(t,xy,manypoints,spline_type)
    else:
        print('Airfoil oversampling algorithm specified is not an available option. Defaulting to "spline".' % ())
        xxyy = interpolator_wrap(t,xy,manypoints,'spline')
    
    # Normalize the airfoil:
    #   correct rotation so that LE is at (0,0) and TE is at (1,0).
    # jcb: Technially, the LE is at the point of max curvature, but that
    # definition can produce situations that break the interpolation step.
    # Instead, we define the LE as the point that is the furthest distance from
    # the TE.
    xyTE = np.array((np.mean([xxyy[0,0],xxyy[-1,0]]), np.mean([xxyy[0,1],xxyy[-1,1]])))
    xxyy = xxyy - np.tile(xyTE,(xxyy.shape[0],1))
    rays = np.hypot(xxyy[:,0],xxyy[:,1]) # distance of each point from the TE
    max_ray = np.max(rays)
    max_point = np.argmax(rays)
    ray_angle = np.arctan2(xxyy[max_point,1],- xxyy[max_point,0])
    xxyy = rotate2d(xxyy,ray_angle)
    xxyy = xxyy / max_ray + np.tile(np.array([1,0]),(xxyy.shape[0],1))

    #Separate into high and low pressure surfaces
    HP = xxyy[0:max_point+1,:] # HP points progress from TE (x=1) to LE (x=0)
    LP = xxyy[max_point:,:] # LP points progress from LE (x=0) to TE (x=1)
    
    # if 'auto', determine which spacing algorithm to use
    if spacing == 'auto':
        dx = xxyy[1,0] - xxyy[2,0]
        # If x-spacing of the oversampled data at the TE is below a threshold,
        # assume that cosine spacing would be best choice for the profile.
        if dx < 1 / (10 * oversample):
            spacing = 'cosine'
        else:
            spacing = 'half-cosine'
    
    #Calculate x points based on spacing algorithm specified
    n_panels = int(np.trunc(n_samples / 2) - 1)
    #NOTE might need to workshop switch cases here
    if 'cosine' == spacing:
        beta = np.linspace(0,np.pi,n_panels + 1)
        x_fwd = 0.5 * (1 - np.cos(beta))
    elif 'half-cosine' == spacing:
        beta = np.linspace(0,np.pi / 2,n_panels + 1)
        x_fwd = (1 - np.cos(beta))
    elif 'constant' == spacing:
        x_fwd = np.linspace(0,1,n_panels + 1)
    else:
        raise Exception('Resampling algorithm specified is not an available option')
    
    x_fwd = x_fwd # make x_fwd a column vector
    x_rev = np.flipud(x_fwd) # x_rev values are nominally 1 to 0
    
    #Calculate interpolated airfoil points. For sharp trailing edge airfoils,
    #the trailing edge point is not repeated
    #NOTE need to address 'extrap' option used in matlab code  below
    LP_new = np.stack((x_fwd,interpolator_wrap(LP[:,0],LP[:,1],x_fwd)),axis=1)
    HP_new = np.stack((x_rev,interpolator_wrap(HP[:,0],HP[:,1],x_rev)),axis=1)
    
    # Make sure that LE point is at (0,0)
    HP_new[-1,:] = np.array([0,0])
    xyTE = np.array([1,0])

    #Assemble the two curves into a continuous line
    af_out = np.concatenate((xyTE.reshape(1,-1),HP_new,
        LP_new[1:,:],xyTE.reshape(1,-1)),axis=0)

    return af_out
    

def getAirfoilNormals(coordinates):
    """Method finds which airfoil is flatback. 
    
    If points are placed
    in flatback region, they are removed for good resampling
    results. Currently this removal only works
    for one point located on TE region. Method also gets the TE tpye for round sections.
    coordinates m by 2 matrix where m is the number of points about
    airfoil. Coordiantes are 2D.

    Parameters
    ----------
    coordinates

    Returns
    -------
    """
    nPoints = coordinates.shape[0]
    unitNormals = np.zeros((nPoints - 1,2))
    for iPoint in range(0,nPoints - 1):
        currentPoint = coordinates[iPoint,:]
        nextPoint = coordinates[iPoint + 1,:]
        r = nextPoint - currentPoint # Postion vector from currentPoint to nextPoint
        if (np.abs(r[0]) + np.abs(r[1])) != 0: # Skip if points are coincedint
            unitNorm = np.transpose(sp.linalg.null_space(r.reshape(1,-1)))
            crossProduct = np.cross(np.concatenate((r,[0])),np.concatenate((unitNorm.reshape(-1),[0])))
            if crossProduct[2] < 0:
                unitNorm = - unitNorm
            unitNormals[iPoint,:] = unitNorm
        else:
            unitNormals[iPoint,:] = np.array([np.nan,np.nan])
    
    return unitNormals
    
    # for iPoint=1:nPoints
    # text(coordinates(iPoint,1),coordinates(iPoint,2),num2str(iPoint),'Color','b')
    # end


def getAirfoilNormalsAngleChange(unitNormals):
    """
    TODO: Docstring
    TODO: Test
    """
    #Find the angle changes between adjacent unit vectors
    nPoints = unitNormals.shape[0]
    angleChange = np.zeros(nPoints)
    for iVector in range(0,nPoints-1):
        currentVector = unitNormals[iVector,:]
        nextVector = unitNormals[iVector + 1,:]
        idotted = np.dot(currentVector,nextVector)
        angleChange[iVector] = np.rad2deg(np.arccos(idotted))
    
    #angle change between last point and first point
    currentVector = unitNormals[-1,:]
    nextVector = unitNormals[0,:]
    dotted = np.dot(currentVector, nextVector)
    angleChange[-1] = np.rad2deg(np.arccos(dotted))
    return angleChange


def rotate2d(xyin, angle): 
    """
    NOTE: might be able to use affinetrans module here
    TODO: Docstring
    TODO: Test
    """
    xyout1 = np.cos(angle) * xyin[:,0] - np.sin(angle) * xyin[:,1]
    xyout2 = np.sin(angle) * xyin[:,0] + np.cos(angle) * xyin[:,1]
    xyout = np.stack((xyout1,xyout2),axis=1)
    return xyout
    

def computeCamberAndThickness(x, y):
    """
    TODO: Docstring
    TODO: Test
    """
    n_samples = len(x)
    LE = int(np.trunc(n_samples / 2))
    xhp = x[LE:0:-1]
    xlp = x[LE:-1]
    assert len(xhp) == len(xlp),'Error computing camber and thickness.'
    #NOTE unsure of how to translate `eps` from matlab -kb
    #assert(sum(xhp - xlp) < np.finfo(1.0),'Upper and lower surface x-coordinates must align.')
    yhp = y[LE:0:-1]
    ylp = y[LE:-1]
    c = xlp
    camber = (yhp + ylp) / 2
    thickness = np.abs(ylp - yhp)
    return c,camber,thickness
        

def readAirfoilColumns(filecontents):
    """
    """
    # All of these file formats assume that the
    # LE is at (0,0) and the TE is at (1,0)
    raw = re.findall('[^\n\r]*',filecontents) # get lines
    
    Nraw = np.asarray(raw).size
    kh = 1  # index counter for header lines
    kt = 1  # index counter for tables
    kr = 1  # index counter for rows
    header = []
    table = []
    for k in range(0,Nraw+1):
        # try to read pairs of coordinates
        pair = cell2mat(textscan(raw[k],'%f %f'))
        if len(pair)==0:
            if kt > 1 or kr > 1:
                # then move to a new table
                kt = kt + 1
                kr = 1
            else:
                # otherwise keep reading the header
                header[kh] = raw(k)
                kh = kh + 1
        else:
            # place coordinate pair in table
            table[kt][kr,:] = pair
            kr = kr + 1
    
    if np.asarray(table).size == 1:
        # assume points wrap around either LE or TE
        coords = table[0]
        if coords[0,0] < 0 or coords[0,0] > 1:
            raise Exception('First x-coordinate not in range 0..1')
        dc = np.diff(coords[:,1])
        dc = dc * np.sign(dc[1])
        k = np.find(dc < 0,1)
        sideA = coords[0:k+1:2]
        sideB = coords[k:-1:2]
        if np.mean(sideA) > np.mean(sideB):
            # LP (upper) surface given first, so flipud
            #             disp('LP first');
            coords = np.flipud(coords)
            k = coords.shape[1-1] - k + 1
        if (1 - coords[0,0]) > 0.5:
            # coordinates begin at LE and wrap around TE
            #             disp('TE wrap');
            if coords[0,:]==coords[-1,:]:
                coords[-1,:] = []
            coords = np.concatenate([
                [coords[k,0,-1,:]],
                [coords[-1,k,-1,:]]])
    else:
        if np.asarray(table).size == 3:
            # assume "Lednicer's" format
            #(see http://www.ae.illinois.edu/m-selig/ads.html)
            npoints = table[0]
            lp = table[2]
            hp = table[3]
            if npoints.shape[1-1] != 1:
                raise Exception('Format similar to "Lednicers", but more than one row found for table sizes')
            if hp[0,:]==lp[0,:]:
                lp[1,:] = []
            coords = np.concatenate([[np.flipud(hp)],[lp]])
        else:
            raise Exception('File format not recognized')
    
    if np.asarray(header).size >= 1:
        reference = header[0]
        for k in np.arange(2,np.asarray(header).size+1).reshape(-1):
            reference = '%s\n%s' % (reference,header[k])
    else:
        reference = ''
    
    return coords,reference
