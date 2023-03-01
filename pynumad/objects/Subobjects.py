########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################

class MatDBentry:
    """A simple class to organize the attributes of a material
    """
    def __init__(self):
        self.type: str = None
        self.name:str = None
        self.reference:str = None
        self.dens:list = None
        self.nuxy:list = None
        self.ex:list = None
        self.ey:list = None
        self.ez:list = None
        self.gxy:list = None
        self.gyz:list = None
        self.gxz:list = None
        self.prxy:list = None
        self.pryz:list = None
        self.prxz:list = None
        self.xten:list = None
        self.xcmp:list = None
        self.yten:list = None
        self.ycmp:list = None
        self.zten:list = None
        self.zcmp:list = None
        self.xy:list = None
        self.yz:list = None
        self.xz:list = None
        self.xycp:list = None
        self.yzcp:list = None
        self.xzcp:list = None
        self.xzit:list = None
        self.xzic:list = None
        self.yzit:list = None
        self.yzic:list = None
        self.g1g2:list = None
        self.etal:list = None
        self.etat:list = None
        self.alp0:list = None
        self.thicknessType:list = None
        self.uniqueLayers:list = None
        self.symmetryType:list = None
        self.layer:list = None


class Layer:
    """A simple class to organize the attributes of a material layer.
    """
    def __init__(self):

        self.layerName:str = None
        self.thicknessA:float = None
        self.thicknessB:float = None
        self.quantity: int = None
        self.theta: float = None


class Shearweb:
    """A simple class to organize the attributes of a Shearweb
    """
    def __init__(self):

        self.Material: str = None
        self.BeginStation:int = None
        self.EndStation:int = None
        self.Corner:list = None


class BOM:
    """A simple class to organize the attributes of a Bill of Materials
    
    Attributes
    ----------

    layernum : int
        Layer
    materialid : int
        Material ID
    name : str
        Component or region name
    beginsta : float
        Begin station (m)
    endsta : float
        End station (m)
    maxwidth : float
        Max width (m)
    avgwidth : float
        Average width (m)
    area : float
        3D area (m^2)
    thickness : float
        Layer thickness (mm)
    weight : float
        Computed dry layer weight (g)
    """
    def __init__(self):

        self.layernum:int = None
        self.materialid:int = None
        self.name:str = None
        self.beginsta:float = None
        self.endsta:float = None
        self.maxwidth:float = None
        self.avgwidth:float = None
        self.area:float = None
        self.thickness:float = None
        self.weight:float = None


class Ply:
    """A simple class to organize the attributes of a ply

    Attributes
    ----------

    component : str
        parent component
    materialid : str
        Material id of ply
    thickness : float
        thickness of single ply (mm)
    angle : float 
        ply angle
    nPlies : int
        number of plies
    """
    def __init__(self):
        self.component: str = None # parent component``
        self.materialid: str = None # materialid of ply``
        self.thickness: float = None # thickness [mm] of single ply``
        self.angle: float = None # ply angle``
        self.nPlies: int = None # number of plies``


class SkinArea:
    def __init__(self):
        self.startIB = []
        self.endIB = []
        self.startOB = []
        self.endOB = []
        self.Material = []
