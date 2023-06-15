########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

import numpy as np
import os
from pynumad.objects.Airfoil import Airfoil
from pynumad.utils.interpolation import interpolator_wrap

class Station():
    """Station object
    
    Parameters
    ----------
    af : AirfoilDef or string
        airfoil

    Attributes
    ----------
    airfoil
        ``AirfoilDef`` object
    spanlocation
        Spanwise location where station is defined [m]
    parent
    hgProfiile
    """
    def __init__(self,af = None):
        self.airfoil = None
        self.spanlocation = None
        self.parent = None
        self.hgProfile = None

        if isinstance(af,str):
            self.airfoil = Airfoil(filename = af)            
        elif isinstance(af, Airfoil):
            self.airfoil = af
        else:
            self.airfoil = Airfoil()        
    
    @property
    def degreestwist(self):
        """
        TODO docstring
        """
        _degreestwist = interpolator_wrap(self.parent.span,self.parent.degreestwist,self.spanlocation)
        return _degreestwist
        
    @property  
    def chord(self): 
        """
        TODO docstring
        """
        _chord = interpolator_wrap(self.parent.span,self.parent.chord,self.spanlocation)
        return _chord
        
    @property    
    def percentthick(self):
        """
        TODO docstring
        """
        _percentthick = interpolator_wrap(self.parent.span,self.parent.percentthick,self.spanlocation)
        return _percentthick
        
    @property    
    def coffset(self):
        """
        TODO docstring
        """
        #             coffset = interp1(self.parent.span,self.parent.coffset,self.spanlocation);
        _coffset = self.airfoil.maxthick
        return _coffset
        
    @property    
    def xyz(self):
        """
        TODO docstring
        """
        twistFlag = - 1
        tratio = self.percentthick / self.airfoil.percentthick
        thick = self.airfoil.thickness * tratio
        hp = self.airfoil.camber - 0.5 * thick
        lp = self.airfoil.camber + 0.5 * thick
        c = self.airfoil.c
        x = np.concatenate((c[-1],np.flipud(c),c[1:],c[-1]),axis = 0)
        y = np.concatenate((hp[-1],np.flipud(hp),lp[1:],lp[-1]),axis = 0)
        x = (x - self.coffset) * self.chord * twistFlag
        y = (y) * self.chord
        twist = twistFlag * self.degreestwist * np.pi / 180
        xyz = np.zeros((len(x),3))
        xyz[:,1] = np.cos(twist) * x - np.sin(twist) * y
        xyz[:,2] = np.sin(twist) * x + np.cos(twist) * y
        xyz[:,3] = self.spanlocation
        return xyz

    
    # NOTE: Not finished. Not sure what this is used for
    def updateProfile(self): 
        xyz = self.xyz
        if len(self.hgProfile)==0:
            self.hgProfile = line(0,0,0)
        set(self.hgProfile,'XData',xyz[:,0],'YData',xyz[:,1],'ZData',xyz[:,2])
        return
        
        
    # def delete(self = None): 
    #     if ishandle(self.hgProfile):
    #         os.delete(self.hgProfile)
        
    #     return