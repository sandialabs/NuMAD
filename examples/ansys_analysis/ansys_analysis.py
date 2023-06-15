import pynumad as pynu
import json

from pynumad.analysis.ansys.mainAnsysAnalysis import mainAnsysAnalysis

# Create a blade object from a yaml file
blade = pynu.Blade()
fileName = 'example_data/myBlade.yaml'
blade.read_yaml(fileName)
blade.mesh = 0.2


# create a shell model in ansys w/o adhesive
meshData=blade.generateShellModel('ansys',includeAdhesive=1)

# Load a previously built loadsTable
with open('myBlade_loadsTable.json','r') as fid:
    loadsTable = json.load(fid)

# Set up configuration for deflection run
analysisConfig = dict()
analysisConfig['meshFile'] = 'master.db'
analysisConfig['analysisFileName'] = 'bladeAnalysis'
analysisConfig['np'] = 1
analysisConfig['analysisFlags'] = dict()
analysisConfig['analysisFlags']['mass'] = 0
analysisConfig['analysisFlags']['deflection'] = 1

ansysResult = mainAnsysAnalysis(blade,meshData,[loadsTable],analysisConfig)