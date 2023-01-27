########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

from matplotlib import pyplot as plt
import numpy as np

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

    See also: ``xlsBlade``, ``BladeDef``, ``BladeDef.addComponent`` 
    """
    def __init__(self):
        self.group: int = None
        self.name: str = None
        self.materialid: str = None
        self.fabricangle: float = None
        self.hpextents: list = None
        self.lpextents: list = None
        self.cp: np.ndarray = None
        self.imethod: str = None
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
        try:
            nLayers = interpolator_wrap(cpx,cpy,span,self.imethod,0)
        finally:
            nLayers = interpolator_wrap(cpx,cpy,span,'linear',0)
        
        return nLayers

    # TODO translate
    def plotcp(self = None):
        """
        TODO docstring
        """
        cpx,cpy = self.getcp()
        self.hCtrl = line(cpx,cpy,'Marker','s','LineStyle','none')
        #             x = linspace(min(cpx),max(cpx),100);
        x = np.linspace(0,1,100)
        y = np.round(interpolator_wrap(cpx,cpy,x,'pchip',0))
        self.hLine = line(x,y,'LineStyle','--')
        set(self.hCtrl,'ButtonDownFcn',bdf_movepts,'UserData',self)
        set(self.hLine,'HitTest','off')
        plt.title(self.name)
        return
    

# non methods
    

#NOTE IGNORED 
#is this gui stuff?
# def bdf_movepts(cbo = None,__ = None): 
#     click = get(gca,'CurrentPoint')
    
#     self = get(cbo,'UserData')
#     pt = click(1,np.arange(1,2+1))
    
#     cpx,cpy = getcp(self)
#     dist = hypot(cpx - pt(1),cpy - pt(2))
#     m,i = np.amin(dist)
    
#     if self.pinnedends and np.any(i == np.array([1,np.asarray(cpx).size])):
#         return
    
#     ax = axis
#     if m < 0.05 * np.amax((ax(2) - ax(1)),(ax(4) - ax(3))):
#         rbbox
#         release = get(gca,'CurrentPoint')
#         pt = release(1,np.arange(1,2+1))
#         if self.pinnedends:
#             pt[1] = np.amin(np.amax(pt(1),0),1)
#         pt[2] = np.amax(pt(2),0)
#         if i > 1:
#             pt[1] = np.amax(pt(1),cpx(i - 1) + eps(cpx(i - 1)))
#         if i < np.asarray(cpx).size:
#             pt[1] = np.amin(pt(1),cpx(i + 1) - eps(cpx(i + 1)))
#         cpx[i] = pt(1)
#         cpy[i] = pt(2)
#         set(self.hCtrl,'XData',cpx,'YData',cpy)
#         x = get(self.hLine,'XData')
#         y = np.round(interp1(cpx,cpy,x,'pchip',0))
#         set(self.hLine,'XData',x,'YData',y)
#         if self.pinnedends:
#             self.cp[i - 1,1] = pt(1)
#             self.cp[i - 1,2] = pt(2)
#         else:
#             self.cp[i,1] = pt(1)
#             self.cp[i,2] = pt(2)
    
#     return