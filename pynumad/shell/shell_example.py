import pynumad as pynu
import numpy as np
import pickle

# blade = pynu.Blade()
# fileName = 'myBlade.yaml'
# blade.read_yaml(fileName)

# blade.mesh = 0.2

adhes = 1

with open('myBlade.obj','rb') as file:
    blade = pickle.load(file)


[nodes,elements,OSSets,SWSets,adNds,adEls] = blade.getShellMesh(adhes)

print(5) #for breakpoint
np.savtext("save_data/nodes.csv", nodes, delimiter=",")
np.savtext("save_data/elements.csv", elements, delimiter=",")
np.savtext("save_data/nodes.csv", nodes, delimiter=",")
