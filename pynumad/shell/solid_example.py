import pynumad as pynu
import numpy as np
import pickle
import os

blade = pynu.objects.Blade.Blade()
fileName = 'myBlade.yaml'
blade.read_yaml(fileName)

blade.expandBladeGeometryTEs()
blade.mesh = 0.2

# with open('myBlade.obj','rb') as file:
    # blade = pickle.load(file)
    
## reduce number of plies
# numSec,numStat = blade.stacks.shape

# for stat in range(0,numStat):
    # for sec in range(0,numSec):
        # for pgi in range(0,len(blade.stacks[sec,stat].plygroups)):
            # newNp = int(np.ceil(0.02*blade.stacks[sec,stat].plygroups[pgi].nPlies))
            # numP = blade.stacks[sec,stat].plygroups[pgi].nPlies
            # if(numP > 50):
                # numP = int(np.ceil(0.02*numP))
                # blade.stacks[sec,stat].plygroups[pgi].nPlies = numP
            # blade.stacks[sec,stat].plygroups[pgi].nPlies = newNp
            

# for webi in range(0,len(blade.swstacks)):
    # for stki in range(0,len(blade.swstacks[webi])):
        # for pgi in range(0,len(blade.swstacks[webi][stki].plygroups)):
            # newNp = int(np.ceil(0.02*blade.swstacks[webi][stki].plygroups[pgi].nPlies))
            # blade.swstacks[webi][stki].plygroups[pgi].nPlies = newNp
            # numP = blade.swstacks[webi][stki].plygroups[pgi].nPlies
            # if(numP > 50):
                # numP = int(np.ceil(0.02*numP))
                # blade.swstacks[webi][stki].plygroups[pgi].nPlies = numP
##

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
    ln = '*Elset, elset=' + es['name'] + '\n'
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
ln = '*Elset, elset=' + es['name'] + '\n'
outFile.write(ln)
for el in es['labels']:
    ln = '  ' + str(el+1) + '\n'
    outFile.write(ln)

outFile.write('*End Part\n')    

outFile.close()    

## End write Abaqus