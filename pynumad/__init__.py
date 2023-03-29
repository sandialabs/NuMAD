from pynumad.objects.Blade import Blade
from pynumad.objects.Airfoil import Airfoil
from pynumad.objects.Component import Component
from pynumad.objects.Material import Material
from pynumad.objects.Stack import Stack
from pynumad.objects.Station import Station
from pynumad.objects.Subobjects import MatDBentry, BOM, Ply, Layer, Shearweb

from os.path import abspath, dirname, join


__version__ = '0.0.1'

__copyright__ = """Copyright 2019 National Technology & Engineering 
Solutions of Sandia, LLC (NTESS). Under the terms of Contract DE-NA0003525 
with NTESS, the U.S. Government retains certain rights in this software."""

__license__ = "Revised BSD License"

import json

dir = dirname(abspath(str(__file__)))
pathdir = join(dir,"..","software_paths.json")
with open(pathdir) as json_file:
    path_data = json.load(json_file)

import platform
path_data["OS"] = platform.system()