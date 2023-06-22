import yaml
import numpy as np
import os

def mesh_to_yaml(meshData,fileName):
    """
    TODO docstring
    """
    mDataOut = dict()
    nodes = list()
    for nd in meshData['nodes']:
        ndstr = str(nd)
        nodes.append(ndstr)
    elements = list()
    for el in meshData['elements']:
        elstr = str(el)
        elements.append(elstr)
    esList = list()
    for es in meshData['sets']['element']:
        newSet = dict()
        newSet['name'] = es['name']
        labels = list()
        for el in es['labels']:
            labels.append(int(el))
        newSet['labels'] = labels
        esList.append(newSet)
    sections = list()
    for sec in meshData['sections']:
        newSec = dict()
        newSec['type'] = sec['type']
        newSec['elementSet'] = sec['elementSet']
        if(sec['type'] == 'shell'):
            newLayup = list()
            for lay in sec['layup']:
                laystr = str(lay)
                newLayup.append(laystr)
            newSec['layup'] = newLayup
        else:
            newSec['material'] = sec['material']
        sections.append(newSec)
    mDataOut['nodes'] = nodes
    mDataOut['elements'] = elements
    mDataOut['sets'] = dict()
    mDataOut['sets']['element'] = esList
    mDataOut['sections'] = sections
    try:
        adNds = list()
        for nd in meshData['adhesiveNds']:
            ndstr = str(nd)
            adNds.append(ndstr)
        adEls = list()
        for el in meshData['adhesiveEls']:
            elstr = str(el)
            adEls.append(elstr)
        mDataOut['adhesiveNds'] = adNds
        mDataOut['adhesiveEls'] = adEls
        mDataOut['adhesiveElSet'] = meshData['adhesiveElSet']
    except:
        pass

    outStream = open('temp.yaml','w')
    yaml.dump(mDataOut,stream=outStream,sort_keys=False)
    outStream.close()
    
    inFile = open('temp.yaml','r')
    outFile = open(fileName,'w')

    fLine = inFile.readline()
    while(fLine != ''):
        newSt = fLine.replace("'","")
        newSt = newSt.replace('"','')
        outFile.write(newSt)
        fLine = inFile.readline()

    inFile.close()
    outFile.close()
    
    os.remove('temp.yaml')
