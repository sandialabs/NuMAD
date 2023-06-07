import numpy as np
from os.path import join
from pynumad.analysis.ansys.write import writeAnsysShellModel
from pynumad.shell.shell import shellMeshGeneral
from pynumad.objects.Blade import Blade

"""
This example loads in a blade, gets the shell mesh, and writes the shell
model to an ansys .src file
"""

yamlpath = join("..","data","blade_yamls","myBlade.yaml")
blade = Blade(yamlpath)
blade.mesh = 0.2
adhes = 1

meshData = shellMeshGeneral(blade, forSolid=0, includeAdhesive=adhes)
config = dict()
config["BoundaryCondition"] = 'cantilevered'
config["elementType"] = '181'
config["MultipleLayerBehavior"] = 'multiply'
config["dbgen"] = 1
config["dbname"] = 'master'

filename = 'buildAnsysShell.src'
includeAdhesive = 1

writeAnsysShellModel(blade,filename,meshData,config,includeAdhesive)