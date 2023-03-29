import numpy as np
from os.path import join

from pynumad.analysis.ansys.ansys import writeANSYSshellModel
from pynumad.analysis.ansys.mainAnsysAnalysis import mainAnsysAnalysis
from pynumad.objects.Blade import Blade


yamlpath = join("..","data","blade_yamls","myBlade.yaml")
blade = Blade(yamlpath)
blade.mesh = 0.2

adhes = 1

meshData = blade.getShellMesh(includeAdhesive=adhes)
config = dict()
config["BoundaryCondition"] = 'cantilevered'
config["elementType"] = '181'
config["MultipleLayerBehavior"] = 'multiply'
config["dbgen"] = 1
config["dbname"] = 'master'
analysisConfig = dict()
# Set up configuration for deflection run
analysisConfig['meshFile'] = 'master.db'
analysisConfig['analysisFileName'] = 'bladeAnalysis'
analysisConfig['np'] = 1
analysisConfig['analysisFlags'] = dict()
analysisConfig['analysisFlags']['mass'] = 0
analysisConfig['analysisFlags']['deflection'] = 1
defLoadsTable = []
ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,analysisConfig)

