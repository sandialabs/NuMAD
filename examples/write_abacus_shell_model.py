import pynumad as pynu
import numpy as np
import os

from pynumad.shell.shell import getShellMesh

## Read blade data from yaml file
blade = pynu.objects.Blade.Blade()
fileName = 'example_data/myBlade.yaml'
blade.read_yaml(fileName)

## Set the airfoil point resolution
for stat in blade.stations:
    stat.airfoil.resample(n_samples=300)
    
blade.updateGeometry()
nStations = blade.geometry.shape[2]
minTELengths = 0.001*np.ones(nStations)
blade.expandBladeGeometryTEs(minTELengths)

## Set the target element size for the mesh
elementSize = 0.2

## Generate mesh
adhes = 1
bladeMesh = getShellMesh(blade, adhes, elementSize)

## Write mesh to yaml
meshFile = 'shellMeshData.yaml'
pynu.io.mesh_to_yaml.mesh_to_yaml(bladeMesh,meshFile)

## Write Abaqus input

outFile = open('shellBlade.inp','w')

outFile.write('*Part, name=Blade\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['nodes']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=S4\n')
i = 1
for el in bladeMesh['elements']:
    if(el[3] != -1):
        el = el + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=S3\n')
i = 1
for el in bladeMesh['elements']:
    if(el[3] == -1):
        el[0:3] = el[0:3] + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
for es in bladeMesh['sets']['element']:
    ln = '*Elset, elset=set' + es['name'] + '\n'
    outFile.write(ln)
    for el in es['labels']:
        ln = '  ' + str(el+1) + '\n'
        outFile.write(ln)

outFile.write('*End Part\n')

outFile.write('*Part, name=Adhesive\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['adhesiveNds']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D8\n')
i = 1
for el in bladeMesh['adhesiveEls']:
    if(el[7] != -1):
        el = el + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),str(el[6]),str(el[7]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D6\n')
i = 1
for el in bladeMesh['adhesiveEls']:
    if(el[7] == -1):
        el[0:6] = el[0:6] + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1

es = bladeMesh['adhesiveElSet']
ln = '*Elset, elset=set' + es['name'] + '\n'
outFile.write(ln)
for el in es['labels']:
    ln = '  ' + str(el+1) + '\n'
    outFile.write(ln)

outFile.write('*End Part\n')    

outFile.close()    

## End write Abaqus
