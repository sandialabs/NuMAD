from cubit import *
from PyCubed_Main import * 
from pynumad.analysis.cubit.cubitUtils import printSineCurveBetweenTwoVerts
import numpy as np
import re

def getOrderedList(partName):


    orderedList=[]
    surfacesToConnect=[1] #Initialize to enter loop
    iSurface=-1 #Initialize 
    while surfacesToConnect:
        iSurface+=1
        parseString=f'with name "*{partName}*surface{iSurface+1}"'
        surfacesToConnect=parse_cubit_list('surface', parseString)
        
        if surfacesToConnect:
            orderedList.append(surfacesToConnect)

    return orderedList
def makeSpanwiseSplines(surfaceDict,orderedList):
    spanwiseSplines=[]
    for alignedSurfaces in orderedList:
        tempList=[]
        for iPoint in range(4):
            vertexList=[]
            for index, iSurface in enumerate(alignedSurfaces):
                vertexID=surfaceDict[iSurface]['verts'][iPoint]
                vertexList.append(vertexID)
                vertexName=cubit.get_entity_name("vertex",vertexID)

            curve=cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
            tempList.append(get_last_id("curve"))
        spanwiseSplines.append(tempList)
    return spanwiseSplines
    
def makeOneVolume(currentSurfaceID,nextSurfaceID,spanwiseSplinesForAvolume,surfaceDict,iStationEnd):
    cubit.cmd(f'surface {currentSurfaceID} copy')
    currentSurfaceIDcopy=get_last_id("surface")
    
    cubit.cmd(f'surface {nextSurfaceID} copy')
    nextSurfaceIDcopy=get_last_id("surface")
    
    currentSurface=cubit.surface(currentSurfaceID)
    nextSurface=cubit.surface(nextSurfaceID)
    
    currentSurfaceCurves=surfaceDict[currentSurfaceID]['curves']
    nextSurfaceCurves=surfaceDict[nextSurfaceID]['curves']
    
    currentSurfaceVerteces=surfaceDict[currentSurfaceID]['verts']
    nextSurfaceVerteces=surfaceDict[nextSurfaceID]['verts'] 
    
    spanwiseSplinesForAvolume.append(spanwiseSplinesForAvolume[0]) #Make list circle back
    
    transverseSurfaceIDs=[]
    for iCurve in range(len(currentSurfaceCurves)):
        cubit.cmd(f'create surface curve {currentSurfaceCurves[iCurve]} {spanwiseSplinesForAvolume[iCurve]} {nextSurfaceCurves[iCurve]} {spanwiseSplinesForAvolume[iCurve+1]}')
        transverseSurfaceIDs.append(get_last_id("surface"))
        
    surfName=cubit.get_entity_name("surface", currentSurface.id()).split('_')
    regex = re.compile('layer')
    layerName = [string for string in surfName if re.match(regex, string)][0]
    stringName=layerName+'_bottomFace'
    cubit.cmd(f'surface {transverseSurfaceIDs[0]} rename "{stringName}"')
    stringName=layerName+'_topFace'
    cubit.cmd(f'surface {transverseSurfaceIDs[2]} rename "{stringName}"')
    
    #cubit.cmd(f'save as "python1.cub" overwrite') 
    #raise Exception(f'Volume "{volumeName}" creation failed')
    #Create Volume
    #nStart=get_last_id("volume")
    cubit.cmd(f'create volume surface {currentSurfaceIDcopy} {nextSurfaceIDcopy} {l2s(transverseSurfaceIDs)} noheal')
    #nEnd=get_last_id("volume")
    #print(f'nStart: {nStart}, nEnd: {nEnd}')

    if 'Station'+str(iStationEnd) in cubit.get_entity_name("surface",nextSurfaceID): #This if statement is needed for componets that may have been droped between the last station and the second to last station
        volumeName=cubit.get_entity_name("surface",nextSurfaceID)
    else:
        volumeName=cubit.get_entity_name("surface",currentSurfaceID)
    if len(cubit.volume(get_last_id("volume")).surfaces())<6:
        print(f'\n\n ERROR with:\n\n create volume surface {currentSurfaceIDcopy} {nextSurfaceIDcopy} {l2s(transverseSurfaceIDs)} ')
        print(f'currentSurfaceIDcopy: {currentSurfaceIDcopy}')
        print(f'nextSurfaceIDcopy: {nextSurfaceIDcopy}')
        print(f'spanwiseSplinesForAvolume: {spanwiseSplinesForAvolume}')
        cubit.cmd(f'save as "python1.cub" overwrite') 
        raise Exception(f'Volume "{volumeName}" creation failed')
    
    volumeName=volumeName.replace('surface','volume')
    cubit.cmd(f'volume {get_last_id("volume")} rename "{volumeName}"')
    
def getspanwiseSplinesForAvolume(iSpan,nCrossSections,spanwiseSplinesForOneSurface,nextSurfaceVerteces):
    #Split off spanwise curves for a single volume and store them
    if iSpan<nCrossSections-2:
        spanwiseSplinesForAvolume=[]
        temp=[]
        for iCurve,curveID in enumerate(spanwiseSplinesForOneSurface):
            cubit.cmd(f'split curve {curveID} at vertex {nextSurfaceVerteces[iCurve]}')
            temp.append(get_last_id("curve"))
            spanwiseSplinesForAvolume.append(get_last_id("curve")-1)
        spanwiseSplinesForOneSurface=temp
    else:
        spanwiseSplinesForAvolume=spanwiseSplinesForOneSurface
    return spanwiseSplinesForAvolume,spanwiseSplinesForOneSurface
# def assignIntervals(volID,nIntervals):
#     thicknessCurveID=cubit.volume(volID).curves()[1].id()
#     #cubit.cmd(f'locate curve {thicknessCurveID} ')
#     cubit.cmd(f'curve {thicknessCurveID} interval {nIntervals}')
    
def makeAeroshell(surfaceDict,orderedList,meshVolList,iStationEnd):
    #nIntervals=3
    spanwiseSplines=makeSpanwiseSplines(surfaceDict,orderedList)
    nCrossSections=len(orderedList[0])
    nPartSurfaceIDs=len(orderedList)
    if nCrossSections>1:
        for iSpan in range(nCrossSections-1):
            for partSurfaceIDs in range(nPartSurfaceIDs):
                currentSurfaceID=orderedList[partSurfaceIDs][iSpan]
                nextSurfaceID=orderedList[partSurfaceIDs][iSpan+1]
                spanwiseSplinesForAvolume,spanwiseSplines[partSurfaceIDs]=getspanwiseSplinesForAvolume(iSpan,nCrossSections,spanwiseSplines[partSurfaceIDs],surfaceDict[nextSurfaceID]['verts'])
                makeOneVolume(currentSurfaceID,nextSurfaceID,spanwiseSplinesForAvolume,surfaceDict,iStationEnd)
                meshVolList.append(get_last_id("volume"))
                #assignIntervals(get_last_id("volume"),nIntervals)


    return meshVolList


def verifyWebCuttingAmplitude(blade,amplitude,tolerance,iStationFirstWeb,iStationLastWeb):
    #Check to make sure that the amplitude does not result sharp volumes by cutting near a station location
    for iStationCheck in range(iStationFirstWeb+1,iStationLastWeb+1):
        bladeSegmentLength=blade.ispan[iStationCheck]-blade.ispan[iStationFirstWeb]
        gap=bladeSegmentLength-amplitude
        #print(f'bladeSegmentLength: {bladeSegmentLength}\ngap {gap}')

        if abs(gap) > tolerance:
            break
        else:
            if gap > 0:
                amplitude=bladeSegmentLength-tolerance
            else:
                amplitude=bladeSegmentLength+tolerance
            break
    #print(f'new amplitude {amplitude} \nnew gap = {bladeSegmentLength-amplitude}')
    return amplitude

def makeBirdsMouth(blade,birdsMouthVerts,amplitudeFraction,iStationFirstWeb,iStationLastWeb):
    ### Make birds mouth volume that will cut the web volumes ###
    #############################################################
    
    #This function must be ran before merging the volumes since "birdsMouthVerts" will change during mergeing

    v1=cubit.vertex(birdsMouthVerts[0])
    v2=cubit.vertex(birdsMouthVerts[1])
    distance=getDist(v1.coordinates(),v2.coordinates())[0] 
    create_curve(v1,v2)

    #Make the birds mouth cut-out start 5% from where the web meets the aeroshell
    cubit.cmd(f'create vertex on curve {get_last_id("curve")}  distance {0.05*distance} from start')
    cubit.cmd(f'create vertex on curve {get_last_id("curve")}  distance {0.05*distance} from end')
    v1=cubit.vertex(get_last_id("vertex")-1)
    v2=cubit.vertex(get_last_id("vertex"))
    straightLine=create_curve(v1,v2)



    amplitude=amplitudeFraction*distance
    tolerance=distance*0.05

    amplitude=verifyWebCuttingAmplitude(blade,amplitude,tolerance,iStationFirstWeb,iStationLastWeb)



    curvedLine=cubit.curve(printSineCurveBetweenTwoVerts(v1.id(),v2.id(),amplitude,'z'))
    cubit.cmd(f'create surface skin curve {curvedLine.id()} {straightLine.id()}')
    baseSurface=get_last_id("surface")

    midPoint=list(curvedLine.position_from_fraction(0.5))
    tangent=straightLine.tangent(midPoint)

    #Get the cross-section normal 
    parseString=f'in surface with name "*webStation{iStationFirstWeb}*"'
    surfaceID=parse_cubit_list('surface', parseString)[0] #Pick the first surface in this list since all on same plane
    coords=cubit.get_center_point("surface", surfaceID)
    surfaceNormal=cubit.surface(surfaceID).normal_at(coords)
    cutBlockLength=5*max(blade.ichord)
    sweepDirection=np.array(vectNorm(crossProd(list(tangent),list(surfaceNormal))))

    cubit.cmd(f'sweep surface {baseSurface} direction {l2s(sweepDirection)} distance {cutBlockLength}')
    cubit.cmd(f'move volume {get_last_id("volume")} x {-cutBlockLength/2*sweepDirection[0]} y {-cutBlockLength/2*sweepDirection[1]} z {-cutBlockLength/2*sweepDirection[2]}')


    cuttingVolume=get_last_id("volume")


    parseString=f'with name "*webStation*"'
    webVolumes=parse_cubit_list('volume', parseString)

    cubit.cmd(f'subtract volume {cuttingVolume} from volume {l2s(webVolumes)}')

    return



#cubit.cmd('open "/home/ecamare/myprojects/bar/cubitDev/python/python0.cub"')

def getApproximateThicknessDirectionForVolume(volumeID):
    #This function is used when assigning material orientations
    
    #Get thickness direction tangents
    approximateThicknessDirection=[]
    for currentCurve in cubit.volume(volumeID).curves():
        curveName=cubit.get_entity_name("curve", currentCurve.id())
        if 'layerThickness' in curveName:
            coords=currentCurve.position_from_fraction(0.5)
            approximateThicknessDirection.append(currentCurve.tangent(coords))
    approximateThicknessDirection=np.array(approximateThicknessDirection)
    nThicknessCurves,_ =approximateThicknessDirection.shape

    if nThicknessCurves==4:  #All other cases
        return np.mean(approximateThicknessDirection,0)
    elif nThicknessCurves ==8: #LE adhesive case and round TE adhesive
        return 0
    elif nThicknessCurves ==6: #Web overwrap 
        # Take the mean of all curves with name 'layerThickness'
        mean=np.mean(approximateThicknessDirection,0)
       
        errorList=[]
        for i in range(nThicknessCurves):
            diff=approximateThicknessDirection[i]-mean

            errorList.append(sqrt(dotProd(diff,diff)))
        sortIndex=np.argsort(errorList)[:4] #Take the first four. This discards the two directions with the largest deviation from the average

        return np.mean(approximateThicknessDirection[sortIndex,:],0)
    else:
        raise ValueError('The number of thickness curves in volume is unexpected. Cannot assign material orientation' ) 
 
    return

def getMatOriSurface(volumeID,spanwiseMatOriCurve):
    #This function is used when assigning material orientations
    #This gets returns the surface within a volume that will be used to get surface normals. 
    #The sign +-1 is also returned since some of the surfaces are oriented the wrong way
    
    approximateThicknessDirection=getApproximateThicknessDirectionForVolume(volumeID)
       
    #Create a list of surface IDs in the given volume
    surfaceIDs=[]
    volumeSurfaces=cubit.volume(volumeID).surfaces()
    for currentSurface in volumeSurfaces:
        surfaceIDs.append(currentSurface.id())
        
    #Eliminate surfaces that have two curves named thickness:
    surfaceCT=0
    for currentSurface in volumeSurfaces:
        curveCT=0 #Counts the number of curves in the surface with name 'layerThickness'
        for currentCurve in currentSurface.curves():
            curveName=cubit.get_entity_name("curve", currentCurve.id())
            if 'layerThickness' in curveName:
                curveCT+=1

        if curveCT>=2:
            surfaceCT+=1
            surfaceIDs.remove(currentSurface.id())
    

    #surfaceIDs now has the list of surfaces w/o thickness curves
    if len(surfaceIDs)==2 or len(surfaceIDs)==1:
        if len(surfaceIDs)==2:
            surfaceName=cubit.get_entity_name("surface", surfaceIDs[0])
            if 'topFace' in surfaceName:
                surfaceID=surfaceIDs[0]
            else:
                surfaceID=surfaceIDs[-1]
        elif len(surfaceIDs)==1: #Web overwrap 
            surfaceID=surfaceIDs[0]
            
        coords=cubit.get_center_point("surface", surfaceID)
        surfaceNormal=cubit.surface(surfaceID).normal_at(coords)

        if dotProd(surfaceNormal,approximateThicknessDirection) >0:
            sign=1.0
        else:
            sign=-1.0
    elif len(surfaceIDs)==0: #LE adhesive and/or TE adhesive for round cross-sections
        #print(f'~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~volumeID {volumeID}')
        surfaceID=False
        sign=1.0

    else:
        raise ValueError('The number of thickness curves in volume is unexpected. Cannot assign material orientation' ) 


    return surfaceID, sign
