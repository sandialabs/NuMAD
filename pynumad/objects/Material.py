########################################################################
#                    Part of the SNL NuMAD Toolbox                     #
#  Developed by Sandia National Laboratories Wind Energy Technologies  #
#              See license.txt for disclaimer information              #
########################################################################


class Material():
    """MaterialDef:  A class for blade materials.
    
    Parameters
    ----------

    Attributes
    ----------
    name : str
		  User selected name of the material
    type : str
		  Two options: isotropic or orthotropic
    layerthickness : float
		  Layer thickness [mm]
    ex : float
		  Longitudinal elastic modulus [Pa]
    ey : float
		  Transverse elastic modulus [Pa]
    ez : float
		Through-the-thickness elastic modulus in the 
      principal material coordinates [Pa]
    gxy : float
		  In-plane shear modulus [Pa]
    gyz : float
		  Transverse shear modulus [Pa]
    gxz : float
		  Transverse shear modulus [Pa]
    prxy : float
		  In-plane Poisson ratio [ ]
    pryz : float
		  Transverse Poisson ratio [ ]
    prxz : float
		  Transverse Poisson ratio [ ]
    density : float
		  Cured mass density [kg/m2]
    drydensity : float
		  Density of fabric
    uts : float
		1 x 3 array of ultimate tensile strength design values.
        Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    ucs : float
		1 x 3 array of ultimate compressive strength design values.
        Sequence: SL , ST, Sz, 1 x 1 for isotropic.
    uss : float
		1 x 3 array of ultimate shear strength design values.
        Sequence: SLT , STz, SLz, 1 x 1 for isotropic.
    xzit : float
		  Lz tensile inclination parameter for Puck failure index
    xzic : float
		  Lz compressive inclination parameter for Puck failure index
    yzit : float
		  Tz tensile inclination parameter for Puck failure index
    yzic : float
		  Tz compressive inclination parameter for Puck failure index
    g1g2 : float
		  Fracture toughness ratio between GI (mode I) and GII (mode II) [ ]
    alp0 : float
		  Fracture angle under pure transverse compression [degrees]
    etat : float
		  Transverse friction coefficient for Larc [ ]
    etal : float
		  Longitudinal friction coefficient for Larc [ ]
    m : list
		  Fatigue slope exponent [ ]
    gamma_mf : list
		  from DNL-GL standard, fatigue strength reduction factor
    gamma_ms : list
		  from DNV-GL standard, short term strength reduction factor
    reference : str = None

    Examples
    --------
    """
    def __init__(self):
        self.name: str = None # User selected name of the material
        self.type: str = None # Two options: ‘isotropic’ or ‘orthotropic’
        self.layerthickness: float = None # Layer thickness [mm]
        self.ex: float = None # Longitudinal elastic modulus [Pa]
        self.ey: float = None # Transverse elastic modulus [Pa]
        self.ez: float = None # Through-the-thickness elastic modulus in the principal material coordinates [Pa]
        self.gxy: float = None # In-plane shear modulus [Pa]
        self.gyz: float = None # Transverse shear modulus [Pa]
        self.gxz: float = None # Transverse shear modulus [Pa]
        self.prxy: float = None # In-plane Poisson ratio [ ]
        self.pryz: float = None # Transverse Poisson ratio [ ]
        self.prxz: float = None # Transverse Poisson ratio [ ]
        self.density: float = None # Cured mass density [kg/m2]
        self.drydensity: float = None # Density of fabric
        self.uts: float = None # 1 × 3 array of ultimate tensile strength design values. Sequence: SL , ST, Sz, 1 × 1 for isotropic.
        self.ucs: float = None # 1 × 3 array of ultimate compressive strength design values. Sequence: SL , ST, Sz, 1 × 1 for isotropic.
        self.uss: float = None # 1 × 3 array of ultimate shear strength design values. Sequence: SLT , STz, SLz, 1 × 1 for isotropic.
        self.xzit: float = None # Lz tensile inclination parameter for Puck failure index
        self.xzic: float = None # Lz compressive inclination parameter for Puck failure index
        self.yzit: float = None # Tz tensile inclination parameter for Puck failure index
        self.yzic: float = None # Tz compressive inclination parameter for Puck failure index
        self.g1g2: float = None # Fracture toughness ratio between GI (mode I) and GII (mode II) [ ]
        self.alp0: float = None # Fracture angle under pure transverse compression [degrees]
        self.etat: float = None # Transverse friction coefficient for Larc [ ]
        self.etal: float = None # Longitudinal friction coefficient for Larc [ ]
        self.m: list = None # Fatigue slope exponent [ ]
        self.gamma_mf: list = None # from DNL-GL standard, fatigue strength reduction factor
        self.gamma_ms: list = None # from DNV-GL standard, short term strength reduction factor
        self.reference: str = None