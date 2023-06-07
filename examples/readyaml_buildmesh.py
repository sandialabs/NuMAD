import pynumad as pynu
import numpy as np
import pickle
from pprint import pprint

from pynumad.shell.shell import getShellMesh

blade = pynu.Blade()
fileName = 'example_data/myBlade.yaml'
blade.read_yaml(fileName)

blade.mesh = 0.2

# with open('example_data/myBlade.obj','rb') as file:
#     blade = pickle.load(file)

adhes = 1

meshData = getShellMesh(blade, adhes)

# Print first 10 nodes coordinates

pprint(meshData['nodes'][:10,:])

# Print first 10 element connectivities

pprint(meshData['elements'][:10,:])

# Print first 10 elements in section 4, spanwise segment 10 (pressure side spar cap)

# meshData['OSSets'][3,9].elementList[:10]

# Print first 10 elements in spanwise segment 10 of first shear web

# meshData['SWSets'][1][10].elementList[10]

foo