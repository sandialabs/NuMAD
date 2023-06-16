from cubit import *
from PyCubed_Main import * 
import numpy as np
import os

def addColor(blade,volumeOrSurface):
    
    #Adds color to volume or surfaces by material
    colorDict={}
    colorDict['adhesive']='yellow'
    colorDict['carbon']='grey'  
    colorDict['uni']='seagreen'  
    colorDict['triax']='lightgreen' 
    colorDict['biax']='greenyellow' 
    colorDict['foam']='khaki'

    for matName in blade.materials:
        for color in colorDict.keys():
            if color in matName.lower():
                parseString=f'with name "*{matName}*"'

                volIDs=parse_cubit_list(volumeOrSurface, parseString)
                cubit.cmd(f'color {volumeOrSurface} {l2s(volIDs)}  mesh {colorDict[color]}')
                cubit.cmd(f'color {volumeOrSurface} {l2s(volIDs)}  geometry {colorDict[color]}')

                break

def surfaceFromTwoCurves(topCurve,bottomCurve):
    v2Left,v2Right=selCurveVerts(topCurve)
    v1Left,v1Right=selCurveVerts(bottomCurve)
    cubit.cmd(f'create curve vertex {v1Left} {v2Left}')
    cubit.cmd(f'create curve vertex {v1Right} {v2Right}')
    cubit.cmd(f'create surface curve {l2s([get_last_id("curve")-1, bottomCurve,topCurve,get_last_id("curve")])}')


def getCrossSectionNormalVector(xyz):
    npts,_=xyz.shape


    #Create Referece line as a spline
    vertexList=[]
    for kcp in range(npts):
        vertexList.append(create_vertex(xyz[kcp,0],xyz[kcp,1],xyz[kcp,2]))

    c1=create_curve(vertexList[0],vertexList[1])
    c2=create_curve(vertexList[0],vertexList[2])

    midPoint=list(c1.position_from_fraction(0.5))
    tangent1=c1.tangent(midPoint)
    midPoint=list(c2.position_from_fraction(0.5))
    tangent2=c2.tangent(midPoint)
    crossSectionNormal=vectNorm(crossProd(list(tangent1),list(tangent2)))
    cubit.cmd(f'delete curve {get_last_id("curve")}')
    cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    return crossSectionNormal
def getBladeGeometryForStation(blade,iStation):
    return np.array([blade.geometry[:,0,iStation],blade.geometry[:,1,iStation],blade.geometry[:,2,iStation]]).transpose()
    # xcoords=blade.profiles[:,0,iStation]*blade.ichord[iStation]
    # ycoords=blade.profiles[:,1,iStation]*blade.ichord[iStation]
    # zcoord=blade.ispan[iStation]
    # zcoords=[zcoord]*len(xcoords)    
    # return np.array([xcoords,ycoords,zcoords]).transpose()


def writeSplineFromCoordinatePoints(cubit,xyz):
    
    #xyz is npts by 3 array holding the coordinates of the points
    npts,_=xyz.shape
    nStart=get_last_id("vertex")+1
    for iPoint in range(npts):
        coords=xyz[iPoint,:]
        create_vertex(coords[0],coords[1],coords[2])
    nEnd=get_last_id("vertex")
    vertexList=list(range(nStart,nEnd+1))
    cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')

def extendCurveAtVertexToLength(curveToExtendID,extensionLength,curveStartOrEnd): 
    #Extend all offset curves
    if curveStartOrEnd.lower() == 'start':
        tempVertID,_=selCurveVerts(curveToExtendID)
        vectorDirectionSign=-1
    else:
        _,tempVertID=selCurveVerts(curveToExtendID)
        vectorDirectionSign=1 
        
    tangentLocationCoords=cubit.vertex(tempVertID).coordinates()
        
    x=cubit.curve(curveToExtendID).tangent(tangentLocationCoords)[0]
    y=cubit.curve(curveToExtendID).tangent(tangentLocationCoords)[1]
    z=cubit.curve(curveToExtendID).tangent(tangentLocationCoords)[2]
    tangentDirection=vectorDirectionSign*extensionLength*np.array(vectNorm([x,y,z]))  #Unit vector of tangent *Scaled by offset curve length
    newVertexCoords=np.array(tangentLocationCoords)+tangentDirection
    
    v1=cubit.create_vertex(newVertexCoords[0],newVertexCoords[1],newVertexCoords[2])
    
    #if statement needed to maintain original curve sense
    if curveStartOrEnd.lower() == 'start':
        c1=cubit.create_curve(v1,cubit.vertex(tempVertID))
        #Combine offset with extenstion
        cubit.cmd(f'create curve combine curve {get_last_id("curve")} {curveToExtendID}')
    else:
        c1=cubit.create_curve(cubit.vertex(tempVertID),v1)
        #Combine offset with extenstion
        cubit.cmd(f'create curve combine curve {curveToExtendID} {get_last_id("curve")}')  
        
    cubit.cmd(f'delete curve {curveToExtendID} {get_last_id("curve")-1}')
    return get_last_id("curve")

def removeBadTEgeometry(blade,iStation,curveID,flatBackCurveID):
    cubit.cmd(f'split curve {curveID} distance {blade.chord[iStation]*0.002} from start ')
    cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    curveID=get_last_id("curve")

    curveStartOrEnd='start'
    extensionLength=1*cubit.curve(curveID).length()
    curveID=extendCurveAtVertexToLength(curveID,extensionLength,curveStartOrEnd)

    _,v1=selCurveVerts(curveID)
    cubit.cmd(f'trim curve {curveID} atintersection curve {flatBackCurveID} keepside vertex {v1}')
    return get_last_id("curve")

def printOffsetDirectionCheck(curveID,lpHpside,crossSectionNormal):
    #This function is used to determine which way a curve offset will go. This is needed, for example, 
    #to make sure the outer mold line curve is being offset towrads the interior of the blade. 
    tol=0.01
    cubit.cmd(f'create curve offset curve {curveID} distance 0.0001')
    tempCurveID=get_last_id("curve")
    cubit.cmd(f'create surface skin curve {curveID} {tempCurveID}')
    cubit.cmd(f'delete curve {tempCurveID}')
    n=get_surface_normal(get_last_id("surface"))
    if abs(n[0]-crossSectionNormal[0]) < tol and abs(n[1]-crossSectionNormal[1]) < tol and abs(n[2]-crossSectionNormal[2]) < tol:
        offsetSign=1
    else:
        offsetSign=-1
    cubit.cmd(f'delete body {get_last_id("body")}')
    
    #Need to flip direction since LP curve is written to run clockwise
    if lpHpside.lower() == 'lp':
        offsetSign=-1*offsetSign    
    return offsetSign

def offsetCurveAndCombineFragmentsIfNeeded(curveID,offsetDistance):
    #Sometimes when offseting a curve, the offset is broken into multiple curves. 
    #This happens when the curvature is to high for the offset distance. If this happens
    #this fuction combines the fragmented curves into one spline.
    nStart=get_last_id("curve")
    cubit.cmd(f'create curve offset curve {curveID} distance {offsetDistance} extended')
    nEnd=get_last_id("curve")
    # print(f'curveID {curveID}')
    # print(f'nStart {nStart}')
    # print(f'nEnd {nEnd}')
    # print(f'offsetDistance {nEnd}')

    # cubit.cmd(f'save as "debug.cub" overwrite')        

    if nEnd-nStart > 1:
        curveList=[]
        curveList+=list(range(nStart+1,nEnd+1))
        vertexList=[]
        v1,_=selCurveVerts(curveList[0])
        vertexList.append(v1)
        for curve in curveList:
            nStart=get_last_id("vertex")+1
            cubit.cmd(f'create vertex on curve {curve} segment 200')
            nEnd=get_last_id("vertex")
            vertexList+=list(range(nStart,nEnd+1))
        _,v1=selCurveVerts(curveList[-1])
        cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
        cubit.cmd(f'delete vertex {l2s(vertexList[1:-1])}')
        cubit.cmd(f'delete curve {l2s(curveList)}')
        


def streamlineCurveIntersections(firstCurveID,secondCurveID,keepCurve):
    #keepCurve is either 1 or 2
    nStart=get_last_id("vertex")
    cubit.cmd(f'create vertex atintersection curve {firstCurveID} {secondCurveID}')
    nEnd=get_last_id("vertex")
    vertexList=list(range(nStart+1,nEnd+1))
    if len(vertexList)>0:
        cubit.cmd(f'split curve {firstCurveID} at vertex {vertexList[-1]}')
        firstCurveID=get_last_id("curve")-1
        cubit.cmd(f'split curve {secondCurveID} at vertex {vertexList[-1]}')
        secondCurveID=get_last_id("curve")
        tempID=spliceTwoCurves(firstCurveID,secondCurveID,keepCurve)
        return tempID, vertexList[-1]
    else:
        if keepCurve==1:
            return firstCurveID,None
        else:
            return secondCurveID,None
def spliceTwoCurves(firstCurveID,secondCurveID,keepCurve):
    #Given two curves splice them into one curve. This should even work for curves that make corners.
    #Curve sence matters. The first curve's sense (tangent) should point towards the second curve.
    cubit.cmd(f'split curve {firstCurveID} distance 0.005 from end ')
    firstCurveID=get_last_id("curve")-1
    cubit.cmd(f'split curve {secondCurveID} distance 0.005 from start ')
    secondCurveID=get_last_id("curve")
    vertexList=[]
    v1,_=selCurveVerts(firstCurveID)
    vertexList.append(v1)
    nStart=get_last_id("vertex")+1
    cubit.cmd(f'create vertex on curve {firstCurveID} segment 200')
    nEnd=get_last_id("vertex")
    vertexList+=list(range(nStart,nEnd+1))
    _,v1=selCurveVerts(firstCurveID)
    vertexList.append(v1)
    v2,_=selCurveVerts(secondCurveID)
    vertexList.append(v2)
    nStart=get_last_id("vertex")+1
    cubit.cmd(f'create vertex on curve {secondCurveID} segment 200')
    nEnd=get_last_id("vertex")
    vertexList+=list(range(nStart,nEnd+1))
    _,v2=selCurveVerts(secondCurveID)
    vertexList.append(v2)
    cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
    cubit.cmd(f'delete curve {firstCurveID} {secondCurveID}')
    secondCurveID=get_last_id("curve")
    cubit.cmd(f'delete vertex {l2s(vertexList[1:-1])}')
    if keepCurve==1:
        return firstCurveID
    elif keepCurve==2:
        return secondCurveID
def extendCurvePastCurveAndTrim(curveToExtendID,curveStartOrEnd,curveIDThatCutsExtendedCurve):
    #Given two curves that are not necessarily intersecting extend {curveToExtendID} then trim it at
    #{curveIDThatCutsExtendedCurve}. {curveStartOrEnd} defines which side of the curve to extend so 
    #you need to know the curve sense

    extensionLength=2*cubit.curve(curveToExtendID).length()
    curveToExtendID=extendCurveAtVertexToLength(curveToExtendID,extensionLength,curveStartOrEnd)  
    
    nStart=get_last_id("vertex")
    cubit.cmd(f'create vertex AtIntersection curve {curveIDThatCutsExtendedCurve}  {curveToExtendID}')
    splitVertexID=get_last_id("vertex")
    if nStart == splitVertexID:
        print(f'curveToExtendID {curveToExtendID} curveIDThatCutsExtendedCurve {curveIDThatCutsExtendedCurve}')

        cubit.cmd(f'save as "debug.cub" overwrite')
        raise Exception(f'Curve {curveToExtendID} was not able to be extended to curve {curveIDThatCutsExtendedCurve} because their intersection was not found.')

    cubit.cmd(f'split curve {curveToExtendID} at vertex {splitVertexID}')
    
    if curveStartOrEnd.lower() == 'start':
        curveToExtendID=get_last_id("curve")
        cubit.cmd(f'delete curve {get_last_id("curve")-1}')
    else:
        curveToExtendID=get_last_id("curve")-1
        cubit.cmd(f'delete curve {get_last_id("curve")}')   
    return curveToExtendID
def renameLastSurface(partName, iStation,iModeledLayers,materialName,partNameID):
    #Every cross sectional surface that is created must be followed by a call to this function
    partNameID+=1
    surfaceName=partName+'Station'+str(iStation)+'_layer'+str(iModeledLayers)+'_'+materialName+'_surface'+str(partNameID)
    cubit.cmd(f'surface {get_last_id("surface")} rename "{surfaceName}"')
    return partNameID

def addSurfaceDictEntry(surfaceDict,surfaceObject,myCurveOrder, myVertOrder,materialName,plyAngle):
    surfaceDict[surfaceObject.id()]={}

    #Curves:
    idList=[]
    for curveObject in surfaceObject.curves():
        idList.append(curveObject.id())
    idList = [idList[i] for i in myCurveOrder]

    surfaceDict[surfaceObject.id()]['curves']=idList
    # getSurfaceCurves[surfaceObject.id()]=idList
    
    #Verts:
    idList=[]
    for vertObject in surfaceObject.vertices():
        idList.append(vertObject.id())
    idList = [idList[i] for i in myVertOrder]
    surfaceDict[surfaceObject.id()]['verts']=idList
    surfaceDict[surfaceObject.id()]['materialName']=materialName
    surfaceDict[surfaceObject.id()]['plyAngle']=plyAngle
    # getSurfaceVerts[surfaceObject.id()]=idList
    
def makeCrossSectionSurface(surfaceDict,iStation,partName,topCurve,bottomCurve,materialName,plyAngle,partNameID,iModeledLayers,materialsUsed):  
    #Given two curves in a cross section, create a surface by connecting the end points then
    #rename the surface and add to the surface dictionary
    surfaceFromTwoCurves(topCurve,bottomCurve)
    materialsUsed.add(materialName)
    
    partNameID=renameLastSurface(partName,iStation,iModeledLayers,materialName,partNameID)
    addSurfaceDictEntry(surfaceDict,cubit.surface(get_last_id("surface")),[0,1,2,3],[0,1,2,3],materialName,plyAngle) 
    
    cubit.cmd(f'curve {surfaceDict[get_last_id("surface")]["curves"][1]} rename "layerThickness"')
    cubit.cmd(f'curve {surfaceDict[get_last_id("surface")]["verts"][-1]} rename "layerThickness"')

    
    return partNameID,materialsUsed


def writeLEAdhesiveCurves(HPstackThickness,LPstackThickness,adhesiveThickness,hpKeyCurve,lpKeyCurve,crossSectionNormal): 
    def uniteLEgeomWithAdhesiveCurves(lpHpside,keyCurve):
        if lpHpside.lower() == 'hp':
            #keyCurve=hpKeyCurve
            adhesiveCurveID=LEadhesiveCurveIDs[0][0]
            offsetCurve=hpOffset
        else:
            #keyCurve=lpKeyCurve
            adhesiveCurveID=LEadhesiveCurveIDs[1][0]
            offsetCurve=lpOffset
            
        print(f'lpHpside: {lpHpside} keyCurve {keyCurve}')

            
        v1,_=selCurveVerts(keyCurve)

        cubit.cmd(f'trim curve {keyCurve} atintersection curve {adhesiveCurveID} keepside vertex {v1}')
        keyCurve=get_last_id("curve")
        
        v1,_=selCurveVerts(offsetCurve)
        cubit.cmd(f'trim curve {offsetCurve} atintersection curve {adhesiveCurveID} keepside vertex {v1}')
        
        offsetCurve=get_last_id("curve")
        _,v1=selCurveVerts(keyCurve)
        _,v2=selCurveVerts(offsetCurve)
        
        cubit.cmd(f'create curve spline vertex {v1} {v2}')
        cubit.cmd(f'delete curve {adhesiveCurveID}')
        cubit.cmd(f'delete curve {offsetCurve}')
        
        return get_last_id("curve"), keyCurve
    
    LEadhesiveCurveIDs=[[],[]]
    
    
    _,p1=selCurveVerts(hpKeyCurve)
    _,p2=selCurveVerts(lpKeyCurve)
    
    coords=[]
    coords.append(list(cubit.vertex(p1).coordinates()))
    coords.append(list(cubit.vertex(p2).coordinates()))
    coords=np.array(coords)
    coords=np.mean(coords,0)
    v1=cubit.create_vertex(coords[0],coords[1],coords[2])
    
    cubit.cmd(f'vertex {get_last_id("vertex")} copy')
    midPointOML=get_last_id("vertex")
    
    #Offset OML curves to final layer offset
    cubit.cmd(f'curve {hpKeyCurve} copy')
    cubit.cmd(f'split curve {get_last_id("curve")} fraction 0.05 from end')
    offsetCurveAndCombineFragmentsIfNeeded(get_last_id("curve"),-1*HPstackThickness)
    #cubit.cmd(f'create curve offset curve {get_last_id("curve")} distance {-1*HPstackThickness} extended')
    _,p1=selCurveVerts(get_last_id("curve"))
    hpOffset=get_last_id("curve")
    
    cubit.cmd(f'curve {lpKeyCurve} copy')
    cubit.cmd(f'split curve {get_last_id("curve")} fraction 0.05 from end')
    offsetCurveAndCombineFragmentsIfNeeded(get_last_id("curve"),-1*LPstackThickness)
    #cubit.cmd(f'create curve offset curve {get_last_id("curve")} distance {-1*LPstackThickness} extended')
    _,p2=selCurveVerts(get_last_id("curve"))
    lpOffset=get_last_id("curve")

    coords=[]
    coords.append(list(cubit.vertex(p1).coordinates()))
    coords.append(list(cubit.vertex(p2).coordinates()))
    coords=np.array(coords)
    coords=np.mean(coords,0)
    v2=cubit.create_vertex(coords[0],coords[1],coords[2])
    c1=cubit.create_curve(v1,v2)
    adhesiveMidLine=get_last_id("curve")
    
    #Extend midline on both sides to make sure other curves eventually intersect with it 
    curveStartOrEnd='start'
    extensionLength=2*cubit.curve(adhesiveMidLine).length()
    adhesiveMidLine=extendCurveAtVertexToLength(adhesiveMidLine,extensionLength,curveStartOrEnd)
    
    curveStartOrEnd='end'
    extensionLength=1*cubit.curve(adhesiveMidLine).length()
    adhesiveMidLine=extendCurveAtVertexToLength(adhesiveMidLine,extensionLength,curveStartOrEnd)
    
    #Copy and move since offset does not seem to work with strait lines
    
    #get offset vector
    
    axialDirection=crossSectionNormal
    position = cubit.curve(adhesiveMidLine).position_from_fraction(1.0)
    tangentDirection=cubit.curve(adhesiveMidLine).tangent(position)
 
    normalDirection=crossProd(axialDirection,tangentDirection)
    normalDirection=adhesiveThickness/2*np.array(vectNorm([normalDirection[0],normalDirection[1],normalDirection[2]]))
    cubit.cmd(f'curve {adhesiveMidLine} copy move x {normalDirection[0]} y {normalDirection[1]} z {normalDirection[2]} nomesh')
    
    LEadhesiveCurveIDs[0].append(get_last_id("curve"))
    normalDirection=-1*normalDirection
    cubit.cmd(f'curve {adhesiveMidLine} copy move x {normalDirection[0]} y {normalDirection[1]} z {normalDirection[2]} nomesh')
    LEadhesiveCurveIDs[1].append(get_last_id("curve"))
    cubit.cmd(f'delete curve {adhesiveMidLine}')

    keyCurves=[hpKeyCurve,lpKeyCurve]
    for iSide, lpHpside in enumerate(['HP','LP']):
        
        ###HP###
        LEadhesiveCurveIDs[iSide][0],keyCurves[iSide]=uniteLEgeomWithAdhesiveCurves(lpHpside,keyCurves[iSide])

        #Make Copies
        cubit.cmd(f'curve {LEadhesiveCurveIDs[iSide][0]} copy')
        LEadhesiveCurveIDs[iSide].append(get_last_id("curve"))

        cubit.cmd(f'curve {LEadhesiveCurveIDs[iSide][1]} copy')
        LEadhesiveCurveIDs[iSide].append(get_last_id("curve"))
        #Extend
        curveStartOrEnd='end'


        extensionLength=1*cubit.curve(LEadhesiveCurveIDs[iSide][1]).length()
        LEadhesiveCurveIDs[iSide][1]=extendCurveAtVertexToLength(LEadhesiveCurveIDs[iSide][1],extensionLength,curveStartOrEnd)



    return keyCurves[0],keyCurves[1],LEadhesiveCurveIDs
def splitCurveAtCoordintePoints(coordinatesToSplitCurve,curveIDToSplit): 
    cubit.cmd(f'curve {curveIDToSplit} copy')
    tempCurveID=get_last_id("curve")
    

    nDPs,_=coordinatesToSplitCurve.shape
    idStart=get_last_id("vertex")+1
    for kcp in range(nDPs):
        temp=cubit.curve(tempCurveID).closest_point([coordinatesToSplitCurve[kcp,0],coordinatesToSplitCurve[kcp,1],coordinatesToSplitCurve[kcp,2]])
        create_vertex(temp[0],temp[1],temp[2])
    idEnd=get_last_id("vertex")
    
    DPverticies = [i for i in range(idStart,idEnd+1)]
    
    idStart=get_last_id("curve")+1
    cubit.cmd(f'split curve {tempCurveID} at vertex {l2s(DPverticies)}')
    idEnd=get_last_id("curve")
    keyCurves = [i for i in range(idStart,idEnd+1)]
    return keyCurves

def splitKeyCurves(keyCurves,aftWebStack,foreWebStack,web_adhesive_width):
    ###Do not split TE reinf
    tempBaseCurveIDs=[keyCurves[0]]

    ###split TE panel in half
    cubit.cmd(f'split curve {keyCurves[1]} fraction 0.5')

    tempBaseCurveIDs.append(get_last_id("curve")-1)
    tempBaseCurveIDs.append(get_last_id("curve"))

    ###Partition sparcap curve
    vertexList=[]
    webLayerThickness=0
    nStart=get_last_id("vertex")+1
    for iLayer in reversed(range(len(aftWebStack.plygroups))):
        webLayerThickness+=aftWebStack.plygroups[iLayer].thickness*aftWebStack.plygroups[iLayer].nPlies/1000
        cubit.cmd(f'create vertex on curve {keyCurves[2]}  distance {webLayerThickness} from start')
    cubit.cmd(f'create vertex on curve {keyCurves[2]}  distance {webLayerThickness+web_adhesive_width} from start')

    #get total foreweb thickness
    webLayerThickness=sum(foreWebStack.layerThicknesses())/1000
    cubit.cmd(f'create vertex on curve {keyCurves[2]}  distance {webLayerThickness+web_adhesive_width} from end')
    for iLayer in reversed(range(len(foreWebStack.plygroups))):
        cubit.cmd(f'create vertex on curve {keyCurves[2]}  distance {webLayerThickness} from end')
        webLayerThickness-=foreWebStack.plygroups[iLayer].thickness*foreWebStack.plygroups[iLayer].nPlies/1000


    nEnd=get_last_id("vertex")
    vertexList+=list(range(nStart,nEnd+1))

    nStart=get_last_id("curve")+1
    cubit.cmd(f'split curve {keyCurves[2]} at vertex {l2s(vertexList)}')
    nEnd=get_last_id("curve")
    tempBaseCurveIDs.append(nStart)
    tempBaseCurveIDs.append(nEnd)
    sparCapBaseCurves=list(range(nStart+1,nEnd))

    ###split LE panel in half
    cubit.cmd(f'split curve {keyCurves[3]} fraction 0.5')
    tempBaseCurveIDs.append(get_last_id("curve")-1)
    tempBaseCurveIDs.append(get_last_id("curve"))

    ###Do not split LE reinf
    tempBaseCurveIDs.append(keyCurves[-1])

    return tempBaseCurveIDs,sparCapBaseCurves
def getMidLine(blade,iLE,iStation,geometryScaling):
    X=blade.geometry[:,0,iStation]* geometryScaling
    Y=blade.geometry[:,1,iStation]* geometryScaling
    Z=blade.geometry[:,2,iStation]* geometryScaling


    ###### Get averge line
    xHP = X[1:iLE]
    xLP = np.flip(X[iLE-1:-1])
    yHP = Y[1:iLE]
    yLP = np.flip(Y[iLE-1:-1])
    zHP = Z[1:iLE]
    zLP = np.flip(Z[iLE-1:-1])
    midline = np.zeros((len(xHP),3))
    for iPoint in range(len(xHP)):
        midline[iPoint,0] = (xHP[iPoint] + xLP[iPoint]) / 2
        midline[iPoint,1] = (yHP[iPoint] + yLP[iPoint]) / 2
        midline[iPoint,2] = (zHP[iPoint] + zLP[iPoint]) / 2

    return midline

def getAdjustmentCurve(curveIDs,layerOffsetDist,curveStartOrEnd,endLayerTaperCurve):

    nStart=get_last_id("vertex")+1
    curveFraction=1.0/3
    for iCurve,curveID in enumerate(curveIDs):
        curveLength=cubit.curve(curveID).length()
        if endLayerTaperCurve is not None and iCurve < endLayerTaperCurve-1:
            if curveLength*curveFraction < layerOffsetDist:
                cubit.cmd(f'create vertex on curve {curveID} fraction {curveFraction} from {curveStartOrEnd}')
            else:
                cubit.cmd(f'create vertex on curve {curveID} distance {layerOffsetDist} from {curveStartOrEnd}')

        else:
            cubit.cmd(f'create vertex on curve {curveID} distance {layerOffsetDist} from {curveStartOrEnd}')

        
    nEnd=get_last_id("vertex")
    vertexList=list(range(nStart,nEnd+1))
    cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
    adjustmentCurve=get_last_id("curve")
    cubit.cmd(f'delete vertex {l2s(vertexList[1:-1])}')
    return adjustmentCurve
def makeCrossSectionLayerAreas_perimeter(surfaceDict,iStation,stationStacks,params,thicknessScaling,lpHpside,isFlatback,TEangle,lastRoundStation,partNameID, nModeledLayers,crossSectionNormal,lpHpCurveDict,materialsUsed):
    partName=lpHpside+'shell'

    #Assumes that #HP side is made first
    if lpHpside.lower() == 'hp':
        lpHpsideIndex=0
        camberOffsetSign=1
    else:
        lpHpsideIndex=1
        camberOffsetSign=-1

    offsetSign_camberID=printOffsetDirectionCheck(lpHpCurveDict['camberID'],'LP',crossSectionNormal)
    
    baseCurveIndexCT=0
    stationStacks.shape
    nStationLayups=len(stationStacks)
    stationStacks.shape

    lastPerimeter=nStationLayups-2

    for iPerimeter in range(nStationLayups-1): #Skip the last stack since the current and the next stack are generated at the same time.
        with open('cubitBlade.log', 'a') as logFile:
            logFile.write(f'\tlpHpside {lpHpside}, iPerimeter={iPerimeter}\n')

        currentStack = stationStacks[iPerimeter]
        nextStack = stationStacks[iPerimeter+1]

        currentStackLayerThicknesses=np.array(currentStack.layerThicknesses())/1000
        nextStackLayerThicknesses=np.array(nextStack.layerThicknesses())/1000

        cubit.cmd(f'curve {lpHpCurveDict["baseCurveIDs"][lpHpsideIndex][baseCurveIndexCT]} copy')
        currentBaseCurveID=get_last_id("curve")
        baseCurveIndexCT+=1
        cubit.cmd(f'curve {lpHpCurveDict["baseCurveIDs"][lpHpsideIndex][baseCurveIndexCT]} copy')
        nextBaseCurveID=get_last_id("curve")
        baseCurveIndexCT+=1

        currentStackSurfaceList=[]
        transitionStackSurfaceList=[]
        nextStackSurfaceList=[]

        currentStackLayerOffset = 0
        nextStackLayerOffset = 0
        layerThicknessTransitionLengths = []
        #Get all offsets and layerThicknessTransitionLengths
        thinestLayerThicknessCurrentStack = 1e+22 #initialize to a large value
        thinestLayerThicknessNextStack = 1e+22
        for iModeledLayers in range(nModeledLayers):
            currentStackLayerOffset+=currentStackLayerThicknesses[iModeledLayers]
            nextStackLayerOffset+=nextStackLayerThicknesses[iModeledLayers]

            adjacentLayerMissmatch=abs(currentStackLayerOffset-nextStackLayerOffset)


            if adjacentLayerMissmatch > params['minimumTransitionLength'][iStation]:
                layerThicknessTransitionLengths.append(adjacentLayerMissmatch / tan(math.radians(params['transitionTaperAngle'])))
            else:
                layerThicknessTransitionLengths.append(params['minimumTransitionLength'][iStation])

            #Also find the thinest layer in stack for meshing purposes
            if currentStackLayerThicknesses[iModeledLayers]>0 and currentStackLayerThicknesses[iModeledLayers] < thinestLayerThicknessCurrentStack:
                thinestLayerThicknessCurrentStack=currentStackLayerThicknesses[iModeledLayers]

            if nextStackLayerThicknesses[iModeledLayers]>0 and nextStackLayerThicknesses[iModeledLayers] < thinestLayerThicknessNextStack:
                thinestLayerThicknessNextStack=nextStackLayerThicknesses[iModeledLayers] 

        maxLayerThicknessTransitionLength=max(layerThicknessTransitionLengths)


        if iPerimeter in [0,2]:
            leftBottomCurve=currentBaseCurveID
            cubit.cmd(f'split curve {nextBaseCurveID} distance {maxLayerThicknessTransitionLength} from start ')
            transitionBottomCurve=get_last_id("curve")-1
            rightBottomCurve=get_last_id("curve")
            transitionStack = nextStack
        elif iPerimeter in [1,3]:
            rightBottomCurve=nextBaseCurveID
            cubit.cmd(f'split curve {currentBaseCurveID} distance {maxLayerThicknessTransitionLength} from end ')
            leftBottomCurve=get_last_id("curve")-1
            transitionBottomCurve=get_last_id("curve")
            transitionStack = currentStack
        else:
            raise ValueError(f'iPerimeter {iPerimeter} not recognized')

        bottomLeftVertexCurveLeft,bottomRightVertexCurveLeft=selCurveVerts(leftBottomCurve)
        bottomLeftVertexCurveRight,bottomRightVertexCurveRight=selCurveVerts(rightBottomCurve)



        #This if statement prepares all layer curves such that they taper at the TE
        if iPerimeter==0:
            currentStackRightCurves=[]
            currentStackLeftCurves=[]

            #Base curve copy
            cubit.cmd(f'curve {leftBottomCurve} copy')
            baseCurveIDCopy=get_last_id("curve")
            offsetSign_baseCurveIDCopy=printOffsetDirectionCheck(baseCurveIDCopy,lpHpside,crossSectionNormal)
            
            #offset camber to make gap
            cubit.cmd(f'create curve offset curve {lpHpCurveDict["camberID"]} distance {camberOffsetSign*offsetSign_camberID*params["TE_adhesive"][iStation]/2} extended')
            camberOffset=get_last_id("curve")

            #Top Bounding Curve
            offsetDistance=1*offsetSign_baseCurveIDCopy*sum(currentStackLayerThicknesses)
            offsetCurveAndCombineFragmentsIfNeeded(baseCurveIDCopy,offsetDistance)
            topBoundingCurve=get_last_id("curve")
               
            if isFlatback:

                curveStartOrEnd='start'
                extensionLength=1*cubit.curve(topBoundingCurve).length()
                topBoundingCurve=extendCurveAtVertexToLength(topBoundingCurve,extensionLength,curveStartOrEnd)
                keepCurve=2

                topBoundingCurve,beginLayerTaperVertexID=streamlineCurveIntersections(camberOffset,topBoundingCurve,keepCurve)

            else:
                lpHpCurveDict['flatBackCurveID']=camberOffset

                curveStartOrEnd='start'
                extensionLength=1*cubit.curve(baseCurveIDCopy).length()
                baseCurveIDCopy=extendCurveAtVertexToLength(baseCurveIDCopy,extensionLength,curveStartOrEnd)

                curveStartOrEnd='start'
                extensionLength=1*cubit.curve(topBoundingCurve).length()
                topBoundingCurve=extendCurveAtVertexToLength(topBoundingCurve,extensionLength,curveStartOrEnd)  

                _,v1=selCurveVerts(baseCurveIDCopy)
                cubit.cmd(f'trim curve {baseCurveIDCopy} atintersection curve {lpHpCurveDict["flatBackCurveID"]} keepside vertex {v1}')
                baseCurveIDCopy=get_last_id("curve")

            #Trim curve at TE.adhesive
            _,v1=selCurveVerts(topBoundingCurve)
            cubit.cmd(f'trim curve {topBoundingCurve} atintersection curve {lpHpCurveDict["flatBackCurveID"]} keepside vertex {v1}')
            topBoundingCurve=get_last_id("curve")
            offsetSign_topBoundingCurve=printOffsetDirectionCheck(topBoundingCurve,lpHpside,crossSectionNormal)

            if isFlatback:# and beginLayerTaperVertexID is not None:

                #Make list of curves that will be used to taper each layer
                npts=30
                nStart=get_last_id("curve")+1
                
                if lpHpside.lower()=='hp':
                    signCorrection=1
                else:
                    signCorrection=-1

                for iPoint in range(npts):
                    if iPoint==0:
                        cubit.cmd(f'create vertex on curve {baseCurveIDCopy} start')
                        cubit.cmd(f'create vertex on curve {topBoundingCurve} start')
                    else:
                        cubit.cmd(f'create vertex on curve {baseCurveIDCopy} fraction {(iPoint)/(npts-1)} from start')
                        cubit.cmd(f'create vertex on curve {topBoundingCurve} fraction {(iPoint)/(npts-1)} from start')
                    cubit.cmd(f'create curve vertex {get_last_id("vertex")-1} {get_last_id("vertex")}')

                nEnd=get_last_id("curve")
                curveIDs=list(range(nStart,nEnd+1))
                
                #If layers are tapeded toward TE find the index in the list of curves (curveIDs) marks the end of the tapering
                if beginLayerTaperVertexID is not None:
                    foundFlag=False
                    for curveID in curveIDs:
                        temp=cubit.curve(curveID).curve_center()
                        tangentDirection=cubit.curve(get_last_id("curve")).tangent(temp)

                        print(f'beginLayerTaperVertexID {beginLayerTaperVertexID}')
                        print(f'baseCurveIDCopy{baseCurveIDCopy} topBoundingCurve {topBoundingCurve}')

                        momentArm=np.array(temp)-np.array(cubit.vertex(beginLayerTaperVertexID).coordinates())

                        normalDirection=signCorrection*np.array(crossProd(tangentDirection,momentArm))

                        if not foundFlag and dotProd(normalDirection,crossSectionNormal)<0:
                            foundFlag=True
                            endLayerTaperCurve=iPoint
                else:
                    endLayerTaperCurve=None


            #TE Adhesive curve
            if isFlatback:
                cubit.cmd(f'curve {topBoundingCurve} copy')
                cubit.cmd(f'split curve {get_last_id("curve")} distance {params["TE_adhesive_width"][iStation]} from start')
                #HPTEadhesiveCurveID=get_last_id("curve")-1
                lpHpCurveDict['TEadhesiveCurveList'][lpHpsideIndex].append(get_last_id("curve")-1)
                cubit.cmd(f'delete curve {get_last_id("curve")}')
            else:
                v1,_=selCurveVerts(baseCurveIDCopy)
                v2,_=selCurveVerts(topBoundingCurve)
                cubit.cmd(f'split curve {lpHpCurveDict["flatBackCurveID"]} at vertex {v1} {v2}')
                cubit.cmd(f'delete curve {get_last_id("curve")} {get_last_id("curve")-2}')
                lpHpCurveDict['flatBackCurveID']=get_last_id("curve")-1
                cubit.cmd(f'curve {lpHpCurveDict["flatBackCurveID"]} copy')
                #HPTEadhesiveCurveID=get_last_id("curve")

            iLayer=0
            offsetDistance=1*offsetSign_baseCurveIDCopy*currentStackLayerThicknesses[iLayer]
            offsetCurveAndCombineFragmentsIfNeeded(baseCurveIDCopy,offsetDistance)
            firstLayerOffset=get_last_id("curve")

            curveStartOrEnd='start'
            firstLayerOffset=extendCurvePastCurveAndTrim(firstLayerOffset,curveStartOrEnd,lpHpCurveDict['flatBackCurveID'])
            
            #Only do the following if all layer thicknesses are unequal

            if isFlatback and abs(min(currentStack.layerThicknesses())-max(currentStack.layerThicknesses())) > 0.0001: 

                layerOffsetDist=currentStackLayerThicknesses[0]
                curveStartOrEnd='start'
                firstLayerOffset=getAdjustmentCurve(curveIDs,layerOffsetDist,curveStartOrEnd,endLayerTaperCurve)



            cubit.cmd(f'create curve offset curve {lpHpCurveDict["camberID"]} distance {camberOffsetSign*offsetSign_camberID*params["TE_adhesive"][iStation]/2} extended')
            camberOffset=get_last_id("curve")
            
            offsetSign_topBoundingCurve=printOffsetDirectionCheck(topBoundingCurve,lpHpside,crossSectionNormal)
            
            offsetDistance=-1*offsetSign_topBoundingCurve*currentStackLayerThicknesses[-1]
            offsetCurveAndCombineFragmentsIfNeeded(topBoundingCurve,offsetDistance)
            lastLayerOffset=get_last_id("curve")


            curveStartOrEnd='start'
            lastLayerOffset=extendCurvePastCurveAndTrim(lastLayerOffset,curveStartOrEnd,lpHpCurveDict['flatBackCurveID'])

            if isFlatback and abs(min(currentStack.layerThicknesses())-max(currentStack.layerThicknesses())) > 0.0001: 

                layerOffsetDist=currentStackLayerThicknesses[-1]
                curveStartOrEnd='end'
                lastLayerOffset=getAdjustmentCurve(curveIDs,layerOffsetDist,curveStartOrEnd,endLayerTaperCurve)

            cubit.cmd(f'split curve {firstLayerOffset} distance {params["TE_adhesive_width"][iStation]} from start')
            currentStackLeftCurves.append(get_last_id("curve")-1)
            currentStackRightCurves.append(get_last_id("curve"))
            curveLen=cubit.curve(currentStackLeftCurves[0]).length()
            cubit.cmd(f'split curve {baseCurveIDCopy} distance {curveLen} from start')
            currentStackLeftCurves.insert(0,get_last_id("curve")-1)
            currentStackRightCurves.insert(0,get_last_id("curve"))
            cubit.cmd(f'split curve {lastLayerOffset} distance {curveLen} from start')
            currentStackLeftCurves.append(get_last_id("curve")-1)
            currentStackRightCurves.append(get_last_id("curve"))
            cubit.cmd(f'split curve {topBoundingCurve} distance {curveLen} from start')
            currentStackLeftCurves.append(get_last_id("curve")-1)
            currentStackRightCurves.append(get_last_id("curve"))


            #### Next Stack (the panel might intersect the camberline so the following is needed 
            nextStackCurves=[]
            cubit.cmd(f'curve {rightBottomCurve} copy')
            baseCurveIDCopy=get_last_id("curve")
            offsetSign_baseCurveIDCopy=printOffsetDirectionCheck(baseCurveIDCopy,lpHpside,crossSectionNormal)
            nextStackCurves.append(baseCurveIDCopy)
            
            ### Offset camber to make gap
            #Offset is increased to create a larger clearance between HP LP shells so that the panels
            #to not self intersect during a simulation (this may not be needed)
            cubit.cmd(f'create curve offset curve {lpHpCurveDict["camberID"]} distance {camberOffsetSign*offsetSign_camberID*0.001*4} extended')
            camberOffset=get_last_id("curve")
            
            iLayer=0
            offsetDistance=1*offsetSign_baseCurveIDCopy*nextStackLayerThicknesses[0]
            offsetCurveAndCombineFragmentsIfNeeded(baseCurveIDCopy,offsetDistance)
            nextStackCurves.append(get_last_id("curve"))


            offsetDistance=1*offsetSign_baseCurveIDCopy*sum(nextStackLayerThicknesses)
            offsetCurveAndCombineFragmentsIfNeeded(baseCurveIDCopy,offsetDistance)
            topBoundingCurve=get_last_id("curve")

            keepCurve=2
            topBoundingCurve,intersectionVertex=streamlineCurveIntersections(camberOffset,topBoundingCurve,keepCurve)

            v1,_=selCurveVerts(baseCurveIDCopy)
            cubit.cmd(f'split curve {topBoundingCurve} at vertex {v1}')
            topBoundingCurve=get_last_id("curve")
            cubit.cmd(f'delete curve {get_last_id("curve")-1}')
            nextStackCurves.append(topBoundingCurve)

            offsetSign_topBoundingCurve=printOffsetDirectionCheck(topBoundingCurve,lpHpside,crossSectionNormal)

            
            offsetDistance=-1*offsetSign_topBoundingCurve*nextStackLayerThicknesses[-1]
            offsetCurveAndCombineFragmentsIfNeeded(topBoundingCurve,offsetDistance)
            lastLayerOffset=get_last_id("curve")
            nextStackCurves.insert(-1,lastLayerOffset)


        for iModeledLayers in range(nModeledLayers):
            currentStackOffset=currentStackLayerThicknesses[iModeledLayers]
            nextStackOffset=nextStackLayerThicknesses[iModeledLayers]


            offsetSign_leftBottomCurve=printOffsetDirectionCheck(leftBottomCurve,lpHpside,crossSectionNormal)
            offsetSign_rightBottomCurve=printOffsetDirectionCheck(rightBottomCurve,lpHpside,crossSectionNormal)

            #Special Treatment for TE
            if iPerimeter==0:
                #Create Left Areas Only

                materialName=currentStack.plygroups[iModeledLayers].materialid
                plyAngle=currentStack.plygroups[iModeledLayers].angle

                partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,currentStackLeftCurves[iModeledLayers+1],currentStackLeftCurves[iModeledLayers],materialName,plyAngle,partNameID,iModeledLayers,materialsUsed)
                
                if TEangle > 120:
                    lpHpCurveDict['TEadhesiveCurveList'][lpHpsideIndex].append(surfaceDict[get_last_id('surface')]['curves'][-1])
                        
                if iStation== lastRoundStation-1:
                    v1=surfaceDict[get_last_id("surface")]['verts'][0]
                    cubit.cmd(f'vertex {v1} rename "linear"')
                    v1=surfaceDict[get_last_id("surface")]['verts'][-1]
                    cubit.cmd(f'vertex {v1} rename "linear"')

                leftBottomCurve=currentStackRightCurves[iModeledLayers]
                leftTopCurve=currentStackRightCurves[iModeledLayers+1]
                [bottomLeftVertexCurveLeft,bottomRightVertexCurveLeft]=selCurveVerts(currentStackRightCurves[iModeledLayers])
                [topLeftVertexCurveLeft,topRightVertexCurveLeft]=selCurveVerts(currentStackRightCurves[iModeledLayers+1])

                #Create Left Areas Only
                rightBottomCurve=nextStackCurves[iModeledLayers]
                rightTopCurve=nextStackCurves[iModeledLayers+1]
                [bottomLeftVertexCurveRight,bottomRightVertexCurveRight]=selCurveVerts(nextStackCurves[iModeledLayers])
                [topLeftVertexCurveRight,topRightVertexCurveRight]=selCurveVerts(nextStackCurves[iModeledLayers+1])
            else:
                cubit.cmd(f'create curve offset curve {leftBottomCurve} distance {offsetSign_leftBottomCurve*currentStackOffset} extended')
                leftTopCurve=get_last_id("curve")
                [topLeftVertexCurveLeft,topRightVertexCurveLeft]=selCurveVerts(get_last_id("curve"))

                offsetCurveAndCombineFragmentsIfNeeded(rightBottomCurve,offsetSign_rightBottomCurve*nextStackOffset)
                lastOffsetCurve=get_last_id("curve")
                rightTopCurve=get_last_id("curve")
                [topLeftVertexCurveRight,topRightVertexCurveRight]=selCurveVerts(get_last_id("curve"))


            if iPerimeter==lastPerimeter:
                curveStartOrEnd='end'
                print(f'iPerimeter {iPerimeter}')
                lastOffsetCurve=extendCurvePastCurveAndTrim(lastOffsetCurve,curveStartOrEnd,lpHpCurveDict['LEadhesiveCurveIDs'][lpHpsideIndex][1])
                rightTopCurve=lastOffsetCurve
                [bottomLeftVertexCurveRight,bottomRightVertexCurveRight]=selCurveVerts(rightBottomCurve)
                [topLeftVertexCurveRight,topRightVertexCurveRight]=selCurveVerts(rightTopCurve)



            cubit.cmd(f'create curve vertex {topRightVertexCurveLeft} {topLeftVertexCurveRight}')
            transitionTopCurve=get_last_id("curve")
            cubit.cmd(f'create curve vertex {bottomLeftVertexCurveLeft} {topLeftVertexCurveLeft}')
            cubit.cmd(f'create curve vertex {bottomRightVertexCurveLeft} {topRightVertexCurveLeft}')
            nCurvesFinal=get_last_id("curve")
            leftSideCurves = [i for i in range(nCurvesFinal-1,nCurvesFinal+1)]
            cubit.cmd(f'create curve vertex {bottomLeftVertexCurveRight} {topLeftVertexCurveRight}')
            cubit.cmd(f'create curve vertex {bottomRightVertexCurveRight} {topRightVertexCurveRight}')
            nCurvesFinal=get_last_id("curve")
            rightSideCurves = [i for i in range(nCurvesFinal-1,nCurvesFinal+1)]

            if iPerimeter==lastPerimeter: 
                if isFlatback:
                    cubit.cmd(f'split curve {lpHpCurveDict["LEadhesiveCurveIDs"][lpHpsideIndex][1]} at vertex {topRightVertexCurveRight}')
                    rightSideCurves[-1]=get_last_id("curve")-1
                    lpHpCurveDict['LEadhesiveCurveIDs'][lpHpsideIndex][1]=get_last_id("curve")
                lpHpCurveDict['LEadhesiveCurveList'][lpHpsideIndex].append(rightSideCurves[-1])


            ### Create Surfaces ###
            #Sufaces for currentStack
            materialName=currentStack.plygroups[iModeledLayers].materialid
            plyAngle=currentStack.plygroups[iModeledLayers].angle
            partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,leftTopCurve,leftBottomCurve,materialName,plyAngle,partNameID,iModeledLayers,materialsUsed)
            currentStackSurfaceList.append(get_last_id("surface"))

            #Surfaces for transitionStack
            materialName=transitionStack.plygroups[iModeledLayers].materialid
            plyAngle=transitionStack.plygroups[iModeledLayers].angle
            partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,transitionTopCurve,transitionBottomCurve,materialName,plyAngle,partNameID,iModeledLayers,materialsUsed)
            transitionStackSurfaceList.append(get_last_id("surface"))
                
            #Surfaces for nextStack
            materialName=nextStack.plygroups[iModeledLayers].materialid
            plyAngle=nextStack.plygroups[iModeledLayers].angle
            partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,rightTopCurve,rightBottomCurve,materialName,plyAngle,partNameID,iModeledLayers,materialsUsed)
            nextStackSurfaceList.append(get_last_id("surface"))

                
            ### Reset ###
            #Reset curves
            leftBottomCurve=leftTopCurve
            transitionBottomCurve=transitionTopCurve
            rightBottomCurve=rightTopCurve

            #Reset vertices
            bottomLeftVertexCurveLeft=topLeftVertexCurveLeft
            bottomRightVertexCurveLeft=topRightVertexCurveLeft
            bottomLeftVertexCurveRight=topLeftVertexCurveRight
            bottomRightVertexCurveRight=topRightVertexCurveRight


        #Build spar caps
        if iPerimeter==1:
            lpHpCurveDict['webInterfaceCurves'][lpHpsideIndex]=[rightTopCurve]
            for ic, currentCurveID in enumerate(lpHpCurveDict['sparCapBaseCurves'][lpHpsideIndex]):
                bottomCurve=currentCurveID
                offSetSign=printOffsetDirectionCheck(bottomCurve,lpHpside,crossSectionNormal)
                
                for it, thickness in enumerate(nextStackLayerThicknesses):
                    cubit.cmd(f'create curve offset curve {bottomCurve} distance {offSetSign*thickness} extended')
                    topCurve=get_last_id("curve")

                    materialName=nextStack.plygroups[it].materialid
                    plyAngle=nextStack.plygroups[it].angle
                    partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,topCurve,bottomCurve,materialName,plyAngle,partNameID,it,materialsUsed)
                    nextStackSurfaceList.append(get_last_id("surface"))     

                    if it==2 and ic!=3:
                        lpHpCurveDict['webInterfaceCurves'][lpHpsideIndex].append(topCurve)
                    bottomCurve=topCurve
        elif iPerimeter==2:
            lpHpCurveDict['webInterfaceCurves'][lpHpsideIndex].append(leftTopCurve)
    return partNameID,lpHpCurveDict
####################################################
####################################################
####################################################
####################################################
####################################################


def createSimplistSurfaceForTEorLEadhesive(iStation,surfaceDict,partName,adhesiveCurveList,adhesiveMatID,partNameID,nModeledLayers,materialsUsed):
    for iCurve in range(len(adhesiveCurveList[0])):
        plyAngle=0 #Ply angle is always zero since adhesive is always assumed as isotropic
        partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName,adhesiveCurveList[1][iCurve],adhesiveCurveList[0][iCurve],adhesiveMatID,plyAngle,partNameID,nModeledLayers+1,materialsUsed)
        
    return partNameID
    
def printSineCurveBetweenTwoVerts(vBot,vTop,amplitude,direction):
    nSineCurveSamplePoints=7
    cubit.cmd(f'create curve vertex {vBot} {vTop}')
    
    idCurve=get_last_id("curve")

    if round(amplitude,3)>0:
        nStart=get_last_id("vertex")+1
        cubit.cmd(f'create vertex on curve {get_last_id("curve")} segment {nSineCurveSamplePoints-1}')
        nEnd=get_last_id("vertex")
        sineCurveSamplePoints=np.linspace(0,pi,nSineCurveSamplePoints)
        vertexOffsets=amplitude*np.sin(sineCurveSamplePoints)
        vertexList=list(range(nStart,nEnd+1))
        for iVert, vertexOffset in enumerate(vertexOffsets[1:-1]): #skip first and last point since those are considered fixed and the offset is zero anyway
            cubit.cmd(f'move vertex {vertexList[iVert]} {direction} {vertexOffset}')
        vertexList.insert(0,vBot)
        vertexList.append(vTop)
        cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
        cubit.cmd(f'delete curve {idCurve}')
    return get_last_id("curve")

def makeCrossSectionLayerAreas_web(surfaceDict,iStation,aftWebStack,foreWebStack,webInterfaceCurves,crosssectionParams,partNameID,crossSectionNormal,nModeledLayers,materialsUsed):
    aftWebOverwrapThickness=(aftWebStack.layerThicknesses()[0]+aftWebStack.layerThicknesses()[-1])/1000
    foreWebOverwrapThickness=(foreWebStack.layerThicknesses()[0]+foreWebStack.layerThicknesses()[-1])/1000
    partName='web'
    ### First create the first two layers. The first layer is the adhesive. The second layer is the web overwrap layer    
    for iCurveList,curveList in enumerate(webInterfaceCurves):
        nBaseCurvesWeb=len(curveList)
        if iCurveList == 0:
            lpHpside='HP'
        else:
            lpHpside='LP'
        for iCurve,bottomCurve in enumerate(curveList):
            
            offSetSign=printOffsetDirectionCheck(bottomCurve,lpHpside,crossSectionNormal) 
            

            if iCurve<nBaseCurvesWeb/2:
                layerThicknesses=[crosssectionParams['Web_aft_adhesive'][iStation],aftWebOverwrapThickness]
            else:
                layerThicknesses=[crosssectionParams['Web_fore_adhesive'][iStation],foreWebOverwrapThickness,]
            


            for it, thickness in enumerate(layerThicknesses):
                cubit.cmd(f'create curve offset curve {bottomCurve} distance {offSetSign*thickness} extended')
                topCurve=get_last_id("curve")

                if it==0:
                    materialName=crosssectionParams['adhesiveMatID']
                    plyAngle=0 
                    
                else:
                    if iCurve<nBaseCurvesWeb/2:
                        materialName=aftWebStack.plygroups[0].materialid
                        plyAngle=aftWebStack.plygroups[0].angle
                    else:
                        materialName=foreWebStack.plygroups[0].materialid
                        plyAngle=foreWebStack.plygroups[0].angle

                partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName, topCurve,bottomCurve,materialName,plyAngle,partNameID,nModeledLayers+it,materialsUsed)

                bottomCurve=topCurve

            #update web interface curves for vertical web retions
            webInterfaceCurves[iCurveList][iCurve]=topCurve

    ### Create vertical web regions  

    #remove curves that are not going to be part of the vertical web
    for iCurveList,curveList in enumerate(webInterfaceCurves):
        curveList.pop(3)
        curveList.pop(3)


    nBaseCurvesWeb=len(webInterfaceCurves[0])
    for iCurve in range(nBaseCurvesWeb):

        vHP,_=selCurveVerts(webInterfaceCurves[0][iCurve])
        vLP,_=selCurveVerts(webInterfaceCurves[1][iCurve])
        topCurve=printSineCurveBetweenTwoVerts(vHP,vLP,crosssectionParams['maxWebImperfectionDistance'][iStation],'x')
        _,vHP=selCurveVerts(webInterfaceCurves[0][iCurve])
        _,vLP=selCurveVerts(webInterfaceCurves[1][iCurve])
        bottomCurve=printSineCurveBetweenTwoVerts(vHP,vLP,crosssectionParams['maxWebImperfectionDistance'][iStation],'x')


        if iCurve<nBaseCurvesWeb/2:
            materialName=aftWebStack.plygroups[iCurve].materialid
            plyAngle=aftWebStack.plygroups[iCurve].angle
        else:
            materialName=foreWebStack.plygroups[iCurve-int(nBaseCurvesWeb/2)].materialid
            plyAngle=foreWebStack.plygroups[iCurve-int(nBaseCurvesWeb/2)].angle
        partNameID,materialsUsed=makeCrossSectionSurface(surfaceDict,iStation,partName, topCurve,bottomCurve,materialName,plyAngle,partNameID,nModeledLayers+it+2+iCurve,materialsUsed)
    return partNameID, (vHP,vLP)

def writeCubitCrossSection(surfaceDict,iStation,iStationGeometry,blade,hasWebs,aftWebStack,foreWebStack,iLE,crosssectionParams,geometryScaling,thicknessScaling,isFlatback,lastRoundStation,materialsUsed):
   
    with open('cubitBlade.log', 'a') as logFile:
        logFile.write(f'Working on Station: {iStation}\n')

    partNameID=0
    crossSectionNormal=getCrossSectionNormalVector(np.array([blade.keypoints[2,:,iStationGeometry],blade.keypoints[3,:,iStationGeometry],blade.keypoints[7,:,iStationGeometry]]))


    #### Step one create outer mold line
    xyz = getBladeGeometryForStation(blade,iStationGeometry) * geometryScaling
    
    #Start indexing from 1 (not 0) to ignore first point: because first point is not on the LP or HP surface but rather is the midpoint at the TE
    splinePoints=xyz[1:iLE,:]
    writeSplineFromCoordinatePoints(cubit,splinePoints)
    hpKeyCurve=get_last_id("curve")

    xyz=np.flip(xyz,0)
    splinePoints=xyz[1:iLE,:]
    writeSplineFromCoordinatePoints(cubit,splinePoints)
    lpKeyCurve=get_last_id("curve")

    # curveID=hpKeyCurve
    # npts=150
    # maxCurvature=0
    # maxCurveFraction=0.25 #Only check curvatue in the first 25% of curve length
    # for iPoint in range(npts):
    #     position = cubit.curve(curveID).position_from_fraction(iPoint/npts*maxCurveFraction)
    #     curvature=math.sqrt(sum(i**2 for i in cubit.curve(curveID).curvature(position)))
    #     print(f'curvature {curvature}')
    #     if curvature > maxCurvature:
    #         maxCurvature=curvature
    #         maxCurvatureCoords=position
    # print(f'maxCurvature {maxCurvature}')
    # create_vertex(position[0],position[1],position[2])

    # create_vertex(maxCurvatureCoords[0],maxCurvatureCoords[1],maxCurvatureCoords[2])

       


    flatBack_vBot,_=selCurveVerts(hpKeyCurve)
    flatBack_vTop,_=selCurveVerts(lpKeyCurve)

    firstPoint = xyz[-2,:]
    secondPont = xyz[1,:]
    flatBackLength = np.linalg.norm(secondPont - firstPoint)


    flatBackCurve=cubit.create_curve(cubit.vertex(flatBack_vBot),cubit.vertex(flatBack_vTop))
    flatBackCurveID=flatBackCurve.id()

    #### Extend flatback ###
    curveStartOrEnd='start'
    extensionLength=100*blade.ichord[iStationGeometry]*cubit.curve(flatBackCurveID).length()
    flatBackCurveID=extendCurveAtVertexToLength(flatBackCurveID,extensionLength,curveStartOrEnd)

    curveStartOrEnd='end'
    extensionLength=0.5*cubit.curve(flatBackCurveID).length()
    flatBackCurveID=extendCurveAtVertexToLength(flatBackCurveID,extensionLength,curveStartOrEnd)

    # #remove beginning portion of curve due to sometimes high curvature (bad geometry fix)
    # hpKeyCurve=removeBadTEgeometry(blade,iStation,hpKeyCurve,flatBackCurveID)
    # lpKeyCurve=removeBadTEgeometry(blade,iStation,lpKeyCurve,flatBackCurveID)


    curveFraction=0
    TEangle = getTEangle(hpKeyCurve,lpKeyCurve,curveFraction)
    print(f'station {iStation}')
    print(f'edgeLength={flatBackLength*1000}')
    print(crosssectionParams)
    print(f'athickness={crosssectionParams["TE_adhesive"][iStation]*1000}') 
    print(f'TE_adhesive_width {crosssectionParams["TE_adhesive_width"][iStation]*1000}')
    print(f'TEangle {TEangle}')
    if isFlatback:
        #Crate camber line
        offsetDistance=0
        npts=100
        spacing=blade.ichord[iStationGeometry]*0.5/npts

        cubit.cmd(f'curve {flatBackCurveID} copy')
        flatBackOffsetCurveID=get_last_id("curve")

        vertexList=[]
        #Do first point manually outside of loop
        flatBack_vBot,_=selCurveVerts(hpKeyCurve)
        flatBack_vTop,_=selCurveVerts(lpKeyCurve)
        coordsHP=np.array(cubit.vertex(flatBack_vBot).coordinates())
        coordsLP=np.array(cubit.vertex(flatBack_vTop).coordinates())
        coords=np.mean(np.vstack((coordsHP,coordsLP)),0)
        cubit.create_vertex(coords[0],coords[1],coords[2])
        vertexList.append(get_last_id("vertex"))
        for iPoint in range(npts):
            offsetDistance+=spacing
            axialDirection=crossSectionNormal
            position = cubit.curve(flatBackOffsetCurveID).position_from_fraction(1.0)
            tangentDirection=cubit.curve(flatBackOffsetCurveID).tangent(position)
        
            normalDirection=crossProd(axialDirection,tangentDirection)
            normalDirection=-1*spacing*np.array(vectNorm([normalDirection[0],normalDirection[1],normalDirection[2]]))
            cubit.cmd(f'curve {flatBackOffsetCurveID} move x {normalDirection[0]} y {normalDirection[1]} z {normalDirection[2]} nomesh')
        
            cubit.cmd(f'create vertex atintersection curve {flatBackOffsetCurveID} {hpKeyCurve}')
            cubit.cmd(f'create vertex atintersection curve {flatBackOffsetCurveID} {lpKeyCurve}')
            
            coordsHP=np.array(cubit.vertex(get_last_id("vertex")-1).coordinates())
            coordsLP=np.array(cubit.vertex(get_last_id("vertex")).coordinates())

            coords=np.mean(np.vstack((coordsHP,coordsLP)),0)
            cubit.create_vertex(coords[0],coords[1],coords[2])
            vertexList.append(get_last_id("vertex"))


        cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
        cubit.cmd(f'delete vertex {l2s(vertexList[1:-1])}')
        cubit.cmd(f'delete curve {flatBackOffsetCurveID}')

        camberID=get_last_id("curve")  
        



    # if isFlatback:
    #     #Make camber line based on geometry modified by "removeBadTEgeometry"
    #     cubit.cmd(f'curve {hpKeyCurve} copy')
    #     cubit.cmd(f'split curve {get_last_id("curve")} fraction 0.5 from start')
    #     camberBaseCurve=get_last_id("curve")-1

    #     vertexList=[]
    #     nStart=get_last_id("vertex")+1
    #     cubit.cmd(f'create vertex on curve {camberBaseCurve} segment 100')
    #     nEnd=get_last_id("vertex")
    #     vertexList+=list(range(nStart,nEnd+1))
    #     _,v1=selCurveVerts(camberBaseCurve) 
    #     vertexList.append(v1)
    #     print(f'vertexList {vertexList}')
    #     nStart=get_last_id("vertex")+1

    #     #Do first point manually outside of loop
    #     flatBack_vBot,_=selCurveVerts(hpKeyCurve)
    #     flatBack_vTop,_=selCurveVerts(lpKeyCurve)
    #     coordsHP=np.array(cubit.vertex(flatBack_vBot).coordinates())
    #     coordsLP=np.array(cubit.vertex(flatBack_vTop).coordinates())
    #     coords=np.mean(np.vstack((coordsHP,coordsLP)),0)
    #     cubit.create_vertex(coords[0],coords[1],coords[2])

    #     for vertexID in vertexList:
    #         coordsHP=np.array(cubit.vertex(vertexID).coordinates())

    #         coordsLP=np.array(cubit.curve(lpKeyCurve).closest_point(coordsHP))
    #         #print(f'lpKeyCurve {lpKeyCurve} vertexID {vertexID}')

    #         coords=np.mean(np.vstack((coordsHP,coordsLP)),0)
    #         cubit.create_vertex(coords[0],coords[1],coords[2])

    

    #     nEnd=get_last_id("vertex")
    #     vertexList=list(range(nStart,nEnd+1))



    #     cubit.cmd(f'create curve spline vertex {l2s(vertexList)}')
    #     cubit.cmd(f'delete vertex {l2s(vertexList[1:-1])}')
    #     camberID=get_last_id("curve")

    
    #     # cubit.cmd(f'save as "debug.cub" overwrite')        
    #     # foo  
    else:
        xyz=getMidLine(blade,iLE,iStationGeometry,geometryScaling)
        npts,_=xyz.shape
        npts=round(npts*0.75) #Only need write first 3/4 of camber line since LE is constructed another way

        splinePoints=xyz[0:npts,:]
        writeSplineFromCoordinatePoints(cubit,splinePoints)
        camberID=get_last_id("curve")
    


    nStacks=len(blade.stacks)

    LEHPstackThickness=sum(blade.stacks[int(nStacks/2.)-1,iStation].layerThicknesses())/1000
    LELPstackThickness=sum(blade.stacks[int(nStacks/2.),iStation].layerThicknesses())/1000

    #Define variables with HP side in index 0, LP side in index 1
    
    lpHpCurveDict={}
    lpHpCurveDict['sparCapBaseCurves']=[[],[]] 
    lpHpCurveDict['webInterfaceCurves']=[[],[]]
    lpHpCurveDict['baseCurveIDs']=[[],[]]
    lpHpCurveDict['TEadhesiveCurveList']=[[],[]]
    lpHpCurveDict['LEadhesiveCurveList']=[[],[]]

    hpKeyCurve,lpKeyCurve,lpHpCurveDict['LEadhesiveCurveIDs']=writeLEAdhesiveCurves(LEHPstackThickness,LELPstackThickness,crosssectionParams['LE_adhesive'][iStation],hpKeyCurve,lpKeyCurve,crossSectionNormal)
    
    keyCurves=splitCurveAtCoordintePoints(blade.keypoints[1:5,:,iStationGeometry],hpKeyCurve)
    web_adhesive_width=crosssectionParams["Web_adhesive_width"][iStation]
    lpHpCurveDict['baseCurveIDs'][0],lpHpCurveDict['sparCapBaseCurves'][0]=splitKeyCurves(keyCurves,aftWebStack,foreWebStack,web_adhesive_width)


    temp=np.flip(blade.keypoints[:,:,iStationGeometry],0)
    keyCurves=splitCurveAtCoordintePoints(temp[1:5,:],lpKeyCurve)
    lpHpCurveDict['baseCurveIDs'][1],lpHpCurveDict['sparCapBaseCurves'][1]=splitKeyCurves(keyCurves,aftWebStack,foreWebStack,web_adhesive_width)


###################
###################
    # xyz=getMidLine(blade,iLE,iStationGeometry,geometryScaling)
    # npts,_=xyz.shape
    # npts=round(npts*0.75) #Only need write first 3/4 of camber line since LE is constructed another way

    # splinePoints=xyz[0:npts,:]
    # writeSplineFromCoordinatePoints(cubit,splinePoints)
    # camberID=get_last_id("curve")
###################
###################


    #Extend
    curveStartOrEnd='start'
    extensionLength=0.5*cubit.curve(camberID).length()
    camberID=extendCurveAtVertexToLength(camberID,extensionLength,curveStartOrEnd)

    lpHpCurveDict['camberID']=camberID
    lpHpCurveDict['flatBackCurveID']=flatBackCurveID


    nModeledLayers=3

    lpHpside='HP'
    
    
    partNameID,lpHpCurveDict=makeCrossSectionLayerAreas_perimeter(surfaceDict,iStation,blade.stacks[1:6,iStation],crosssectionParams,thicknessScaling,lpHpside,isFlatback,TEangle,lastRoundStation,partNameID,nModeledLayers,crossSectionNormal,lpHpCurveDict,materialsUsed)


    lpHpside='LP'
    temp=blade.stacks[:,iStation]
    temp=np.flip(temp)
    partNameID,lpHpCurveDict=makeCrossSectionLayerAreas_perimeter(surfaceDict,iStation,temp[1:6],crosssectionParams,thicknessScaling,lpHpside,isFlatback,TEangle,lastRoundStation,partNameID,nModeledLayers,crossSectionNormal,lpHpCurveDict,materialsUsed)
    
    partName='shell'
    partNameID=createSimplistSurfaceForTEorLEadhesive(iStation,surfaceDict,partName,lpHpCurveDict['LEadhesiveCurveList'],crosssectionParams['adhesiveMatID'],partNameID,nModeledLayers,materialsUsed)
    
    

    partNameID=0 #Reset since outer areoshell is complete (LE adhesive is accouted for as aeroshell)
    if iStation>lastRoundStation:    
        partName='flatTEadhesive'
    else:
        partName='roundTEadhesive'
        
    partNameID=createSimplistSurfaceForTEorLEadhesive(iStation,surfaceDict,partName,lpHpCurveDict['TEadhesiveCurveList'],crosssectionParams['adhesiveMatID'],partNameID,nModeledLayers,materialsUsed)



    birdsMouthVerts=[]
    if hasWebs:
        partNameID=0 #Reset since outer areoshell is complete (LE adhesive is accouted for as aeroshell)

        partNameID,birdsMouthVerts=makeCrossSectionLayerAreas_web(surfaceDict,iStation,aftWebStack,foreWebStack,lpHpCurveDict['webInterfaceCurves'],crosssectionParams,partNameID,crossSectionNormal,nModeledLayers,materialsUsed)   
    
    


    parseString=f'with name "*station{iStation}*"'
    crossSectionalSurfaces=parse_cubit_list('surface', parseString)
    for surfaceID in crossSectionalSurfaces:
        n=get_surface_normal(surfaceID)
        if n[2]<0:
            cubit.cmd(f'reverse surface {surfaceID}')

    cubit.cmd(f'delete vertex all with Is_Free')

    return birdsMouthVerts



def writeVABSinput(surfaceDict,blade,crosssectionParams,directory,fileName, surfaceIDs,materialsUsed):
    #Write VABS inputfile
    if crosssectionParams['elementShape'].lower() == 'quad':
        expandedConnectivityString='face'
    elif crosssectionParams['elementShape'].lower() == 'tri':
        expandedConnectivityString='tri'
    else:
        raise NameError(f'Element type: {crosssectionParams["elementShape"]} not supported' ) 




    ######Write VABS input file
    nnodes=get_node_count()
    nelem=get_element_count()
    nmate=len(materialsUsed)
    nlayer=len(surfaceIDs) #One VABS layer is defined for each surface
    

    pathName=directory+'/'+fileName

    if not os.path.exists(directory):
        os.makedirs(directory)
        
    with open(pathName, 'w') as f:
        f.write(f'1 {nlayer}\n') #New format. One layer definiton for each element.
        f.write('1 0 0\n')
        f.write('0 0 0 0\n')
        f.write(f'{nnodes} {nelem} {nmate}\n')
        #Write Nodes
        for iNd in range(nnodes):
            nodeID=iNd+1
            coords=list(get_nodal_coordinates(nodeID))
            f.write(f'{nodeID} {coords[0]} {coords[1]}\n')
        f.write('\n\n')
        #Write Elements
        maxNumberOfPossibleNodes=9
        for iEl in range(nelem):
            elementID=iEl+1
            nodesIDs=cubit.get_expanded_connectivity(crosssectionParams['elementShape'],elementID)
            
            if nodesIDs[0]==0 or nodesIDs[0]==0.0:
                foo
            tempStr=str(nodesIDs) #convert tuple to string
            tempStr=tempStr.replace('(','')
            tempStr=tempStr.replace(')','')
            tempStr=tempStr.replace(',',' ')
            nZeros=maxNumberOfPossibleNodes-len(nodesIDs)
            tempStr2=str([0]*nZeros)
            tempStr2=tempStr2.replace('[','')
            tempStr2=tempStr2.replace(']','')
            tempStr2=tempStr2.replace(',',' ')
            f.write(f'{elementID} {tempStr} {tempStr2}\n') 
        #Write ply angle for all but the TE adhesive

        for iSurface, surfaceID in enumerate(surfaceIDs):
            for iEl,elementID in enumerate(get_surface_quads(surfaceID)):
                nodesIDs=cubit.get_expanded_connectivity('face',elementID)
                coords=[]
                for iNd,nodeID in enumerate(nodesIDs):
                    coords.append(list(get_nodal_coordinates(nodeID)))
                coords=np.array(coords)
                # #######For Plotting - find the larges element side length #######
                # distances=[]
                # for iNd,nodeID in enumerate(nodesIDs):
                    # for jNd,nodeIDj in enumerate(nodesIDs):
                        # distances.append(norm(vectSub(coords[iNd],coords[jNd])))
                # length=max(distances)
                # #######For Plotting - find the larges element side length #######
                coords=np.mean(coords,0)
                #coords=cubit.get_center_point(crosssectionParams['elementShape'], elementID)


    #                             minDist=inf #initialize
    #                             closestCurveID=nan #initialize
    #                             #Since there are possibly many curves due to the offset operation, see which curve is closeset to element center
    #                             for iCurve, curveID in enumerate(curves):
    #                                 temp=cubit.curve(curveID).closest_point(coords)
    #                                 distance=getDist(coords,temp)[0]
    #                                 if distance < minDist:
    #                                     minDist=distance
    #                                     closestCurveID=curveID

                curveIDforMatOri=cubit.surface(surfaceID).curves()[0]
                curveLocationForTangent=curveIDforMatOri.closest_point(coords)
                x=curveIDforMatOri.tangent(curveLocationForTangent)[0]
                y=curveIDforMatOri.tangent(curveLocationForTangent)[1]
                z=curveIDforMatOri.tangent(curveLocationForTangent)[2]
                tangentDirection=vectNorm([x,y,z])  #Unit vector of tangent
                theta1=math.atan2(tangentDirection[1],tangentDirection[0])*180/pi
                f.write(f'{elementID} {iSurface+1} {theta1}\n')
                # #######Only needed For Plotting Orientation Check#######
                # cubit.create_vertex(coords[0],coords[1],coords[2])
                # iVert1=get_last_id("vertex")
                # cubit.create_vertex(coords[0]+length*tangentDirection[0],coords[1]+length*tangentDirection[1],coords[2]+length*tangentDirection[2])
                # iVert2=get_last_id("vertex")
                # cubit.cmd(f'create curve vertex {iVert1} {iVert2}')
                # #######Only needed For Plotting Orientation Check#######
            ####Normal to curve
                # axialDirection=crossSectionNormal #There will be a slight error here for highly tapeded regions
                # normalDirection=crossProd(axialDirection,tangentDirection)
                # #######Only needed For Plotting Orientation Check#######
                # cubit.create_vertex(coords[0]+length*normalDirection[0],coords[1]+length*normalDirection[1],coords[2]+length*normalDirection[2])
                # c
                # cubit.cmd(f'create curve vertex {iVert1} {iVert2}')
                # #######Only needed For Plotting Orientation Check#######
        #Define Plies
        for iSurface, surfaceID in enumerate(surfaceIDs):
            materialID=list(materialsUsed).index(surfaceDict[surfaceID]['materialName'])+1
            plyAngle=surfaceDict[surfaceID]['plyAngle']
            f.write(f'{iSurface+1} {materialID} {plyAngle}\n') 
        #Define Materials
        for imat,matName in enumerate(materialsUsed):
            materialID=imat+1
            material=blade.materials[matName]
            f.write(f'{materialID} {1} \n')
            f.write(f'{material.ex} {material.ey} {material.ez}\n')
            f.write(f'{material.gxy} {material.gxz} {material.gyz}\n')
            f.write(f'{material.prxy} {material.prxz} {material.pryz}\n')
            f.write(f'{material.density}\n')
    print('Done writing VABS input')
    return
#Main script fuctions

def getTEangle(hpKeyCurve,lpKeyCurve,fraction):
    c1=cubit.curve(hpKeyCurve)
    c2=cubit.curve(lpKeyCurve)

    coords=list(c1.position_from_fraction(fraction))
    v1=np.array(c1.tangent(coords))
    coords=list(c2.position_from_fraction(fraction))
    v2=np.array(c2.tangent(coords))
    
    return math.degrees(math.acos(v1.dot(np.transpose(v2)) / (np.linalg.norm(v1) * np.linalg.norm(v2))))

