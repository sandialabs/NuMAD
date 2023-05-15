########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

from matplotlib import pyplot as plt
import numpy as np
import matplotlib.pyplot as plt

from pynumad.utils.interpolation import interpolator_wrap


class Component:
    """ComponentDef:  A class definition for blade components.

    Parameters
    ----------
    None

    Attributes
    ----------
    group : int
        0 = blade, 1 = first shear web, 2 = second shear web, etc.
    name : str
        Name, such as 'spar'
    materialid : str
        Material id number from blade.materials
    fabricangle : float
        Fiber angle
    hpextents : list
        Array of keypoints such as {'b','c'}
    lpextents : list
        String Array: Array of keypoints such as {'b','c'}
    cp : np
        control points defining layer distribution
    imethod: str
        imethod
    pinnedends
    hCtrl
    hLine

    Examples: 
     
    	``comp_self = ComponentDef();``
     
    	``comp_self = ComponentDef(comp_struct);``

    """
    
    def __init__(self):
        self.group: int = None
        self.name: str = None
        self.materialid: str = None
        self.fabricangle: float = None
        self.hpextents: list = None
        self.lpextents: list = None
        self.cp: np.ndarray = None
        self.imethod: str = 'linear'
        self.pinnedends: bool = None
        self.hCtrl = None
        self.hLine = None


    def getcp(self): 
        if self.pinnedends:
            if np.any(self.cp[:,0] < 0) or np.any(self.cp[:,0] > 1):
                raise Exception('ComponentDef: first coordinate of control points must be in range [0,1] when using "pinned" ends')
            cpx = np.concatenate(([-0.01],self.cp[:,0],[1.01]))
            cpy = np.concatenate(([0],self.cp[:,1],[0]))
        else:
            cpx = self.cp[:,0]
            cpy = self.cp[:,1]
        
        return cpx,cpy
        
        
    def getNumLayers(self,span): 
        cpx,cpy = self.getcp()
        nLayers = interpolator_wrap(cpx,cpy,span,self.imethod,0)
        return nLayers

    # TODO translate
    def plotcp(self):
        """
        TODO docstring
        """
        cpx,cpy = self.getcp()
        fig, ax = plt.subplots()
        ax.plot(cpx,cpy)
        x = np.linspace(0,1,100)
        y = np.round(interpolator_wrap(cpx,cpy,x,'pchip',0))
        ax.plot(x,y)
        plt.title(self.name)
        return