import pynumad as pynu
import numpy as np
import pickle
import os

blade = pynu.objects.Blade.Blade()
fileName = 'myBlade.yaml'
blade.read_yaml(fileName)

for stat in blade.stations:
    stat.airfoil.resample(n_samples=300)
    
blade.updateGeometry()

blade.expandBladeGeometryTEs()
blade.mesh = 0.2

# with open('myBlade.obj','rb') as file:
    # blade = pickle.load(file)

bladeMesh = pynu.shell.shell.getSolidMesh(blade, layerNumEls=[1,1,1])

## Write Abaqus input

outFile = open('solidBlade.inp','w')

outFile.write('*Part, name=Blade\n')
outFile.write('*Node\n')

i = 1
for nd in bladeMesh['nodes']:
    lst = [str(i),str(nd[0]),str(nd[1]),str(nd[2]),'\n']
    ln = ', '.join(lst)
    outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D8\n')
i = 1
for el in bladeMesh['elements']:
    if(el[7] != -1):
        el = el + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),str(el[6]),str(el[7]),'\n']
        ln = ', '.join(lst)
        outFile.write(ln)
    i = i + 1
    
outFile.write('*Element, type=C3D6\n')
i = 1
for el in bladeMesh['elements']:
    if(el[7] == -1):
        el[0:6] = el[0:6] + 1
        lst = [str(i),str(el[0]),str(el[1]),str(el[2]),str(el[3]),str(el[4]),str(el[5]),'\n']
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