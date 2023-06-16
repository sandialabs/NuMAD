from pynumad.analysis.cubit.cubitUtils import *
from pynumad.analysis.cubit.solidModelUtils import *

import numpy as np
import os
import glob


def generateCubitCrossSections(blade, wt_name, settings, crosssectionParams, model2Dor3D, stationList=None, directory='.'):

    if stationList is None or len(stationList)==0:
        stationList = list(range(len(blade.ispan)))

    # Initialize variables
    surfaceDict = {}
    # Uniquly track which materiall IDs are actuall used in blade model
    materialsUsed = set()
    iLE = blade.LEindex+1
    thicknessScaling = 0.001
    geometryScaling = thicknessScaling*1000

    # Set up Cubit
    cubit.init(['cubit','-nojournal'])
    
    cubit.cmd('undo off')
    cubit.cmd('set geometry accuracy 1e-6')
    # making numerus 3D volumes is very slow with autosize on
    cubit.cmd('set default autosize off')

    # Modify blade object to accomodate actual layer thicknesses
    
    expandTEthicknesses=list(crosssectionParams['TE_adhesive']+6*crosssectionParams['minimumLayerThickness'])
    blade.expandBladeGeometryTEs(expandTEthicknesses)

    

    blade.editStacksForSolidMesh()

    hasWebs = []
    webNumber = 1
    for iStation in range(len(blade.swstacks[webNumber])):
        if not len(blade.swstacks[webNumber][iStation].plygroups) == 0:
            hasWebs.append(True)
        else:
            hasWebs.append(False)

    # WARNING - Last station never has webs. Fix later
    hasWebs.append(False)
    # WARNING - Last station never has webs. Fix later

    # Create Referece line as a spline

    refLineCoords = np.vstack(([blade.sweep, blade.prebend, blade.ispan])).transpose()
    spanwiseMatOriCurve = 1
    
    roundStations=np.argwhere(np.array(blade.TEtype)=='round')
    roundStations=list(roundStations[:,0])
    lastRoundStation = roundStations[-1]


    with open('cubitBlade.log', 'w') as logFile:
        logFile.write(f'Making cross sections for {wt_name}\n')


    pathName=directory+'/'+wt_name+'-crossSections'
    
    for iStation in stationList:

        cubit.cmd('reset ') # This is needed to restart node numbering for VABS. VABS neeeds every element and node starting from 1 to nelem/nnode should be present
        writeSplineFromCoordinatePoints(cubit, refLineCoords)
        iStationGeometry = iStation
        if iStation == len(blade.ispan)-1:  # Only do this for the last station
            blade.addInterpolatedStation(blade.ispan[-1]*0.999)
            blade.editStacksForSolidMesh()
            expandTEthicknesses.append(expandTEthicknesses[-1])
            blade.expandBladeGeometryTEs(expandTEthicknesses)


            # adjustLastStackAfterNewTipStation(iStation)

            iStationGeometry = iStation+1

        if blade.getprofileTEtype(iStationGeometry) == 'flat':
            isFlatback = True
        else:
            isFlatback = False

        iStationFirstWeb = np.argwhere(hasWebs)[0][0]
        iStationLastWeb = np.argwhere(hasWebs)[-1][0]

        if hasWebs[iStation] == True:
            webNumber = 1
            aftWebStack = blade.swstacks[webNumber][iStation]
            webNumber = 0
            foreWebStack = blade.swstacks[webNumber][iStation]
        else:
            if iStation < iStationFirstWeb:
                iWebStation = iStationFirstWeb

    #         elif iStationLastWeb == len(blade.ispan) - 1-1:
            else:
                iWebStation = iStationLastWeb
    #         else:
    #             raise Exception('assuming web ends at last station for now. ')

            webNumber = 1
            aftWebStack = blade.swstacks[webNumber][iWebStation]
            webNumber = 0
            foreWebStack = blade.swstacks[webNumber][iWebStation]


        # Only save birdsMouthVerts for the right cross-section
        if iStation == iStationFirstWeb:
            birdsMouthVerts = writeCubitCrossSection(surfaceDict, iStation, iStationGeometry, blade,
                                                     hasWebs[iStation], aftWebStack, foreWebStack, iLE, crosssectionParams, geometryScaling, thicknessScaling, isFlatback, lastRoundStation, materialsUsed)
        else:

            writeCubitCrossSection(surfaceDict, iStation, iStationGeometry, blade, hasWebs[iStation], aftWebStack, foreWebStack,
                                   iLE, crosssectionParams, geometryScaling, thicknessScaling, isFlatback, lastRoundStation, materialsUsed)
            birdsMouthVerts = []

        cubit.cmd(f'delete curve all with Is_Free except {spanwiseMatOriCurve}')

        # Chord line for rotation of cross-section for homogenization
        if model2Dor3D.lower() == '2d':
            #         #Blocks

            for imat, materialName in enumerate(materialsUsed):
                cubit.cmd(f'block {imat+1} add surface with name "*{materialName}*"')

            addColor(blade, 'surface')

            # create_vertex(blade.geometry[0,0,iStation]*geometryScaling,blade.geometry[0,1,iStation]*geometryScaling,blade.geometry[0,2,iStation]*geometryScaling)
            # TEvert=get_last_id("vertex")
            # create_vertex(blade.geometry[iLE-1,0,iStation]*geometryScaling,blade.geometry[iLE-1,1,iStation]*geometryScaling,blade.geometry[iLE-1,2,iStation]*geometryScaling)
            # LEvert=get_last_id("vertex")

            # cubit.cmd(f'create curve vertex {TEvert} {LEvert}')
            # coords=cubit.vertex(TEvert).coordinates()
            # tangent=cubit.curve(get_last_id("curve")).tangent(coords)
            # tangentDirection=vectNorm(list(tangent))  #Unit vector of tangent.
            # crossSectionRotationAngle=math.atan2(tangentDirection[1],tangentDirection[0])*180/pi

            parseString = f'with name "*Station{str(iStation)}*"'
            volumeIDs = parse_cubit_list('surface', parseString)

            # Undo initial twist
            cubit.cmd(f'rotate Surface {l2s(volumeIDs)} angle {blade.degreestwist[iStation]} about Z include_merged ')

            # Undo prebend
            if blade.prebend[iStation] != 0:
                cubit.cmd(f'move surface {l2s(volumeIDs)} y {-1*blade.prebend[iStation]} include_merged')

            # Undo sweep
            if blade.sweep[iStation] != 0:
                raise ValueError('Presweep is untested for cross-sectional meshing')

            # Mesh the cross-section
            cubit.cmd(f'curve with name "layerThickness*" interval {crosssectionParams["nelPerLayer"]}')
            #cubit.cmd(f'imprint volume {l2s(surfaceIDs)}')
            cubit.cmd(f'merge volume {l2s(volumeIDs)}')
            cubit.cmd(f'set default autosize on')

            if crosssectionParams['elementShape'].lower() == 'tri':
                cubit.cmd(f'surface {l2s(volumeIDs)} scheme tri')
            else:
                cubit.cmd(f'surface {l2s(volumeIDs)} scheme map')

            cubit.cmd(f'mesh surface {l2s(volumeIDs)}')

            fileName = wt_name+'-'+str(iStation)+'-t-0.in'

            if not os.path.exists(directory):
                os.makedirs(directory)


            if settings['export'] is not None:
                if 'g' in settings['export'].lower():
                    cubit.cmd(f'export mesh "{pathName}.g" overwrite')
                elif 'cub' in settings['export'].lower():
                    cubit.cmd(f'delete curve {spanwiseMatOriCurve}')
                    cubit.cmd(f'save as "{pathName}-{str(iStation)}.cub" overwrite')
                else:
                    raise NameError(f'Unknown model export format: {settings["export"]}')
            
            if settings['solver'] is not None:
                if  'vabs' in settings['solver'].lower():
                    writeVABSinput(surfaceDict, blade, crosssectionParams,directory,fileName, volumeIDs,materialsUsed)
           

                elif 'anba' in settings['solver'].lower():
                    raise ValueError('ANBA currently not supported')
                else:
                    raise NameError(f'Unknown beam cross-sectional solver: {settings["solver"]}')

    #Import all cross-sections into one cub file
    if settings['export'] is not None and'cub' in settings['export'].lower():
        cubit.cmd('reset ') 
        writeSplineFromCoordinatePoints(cubit, refLineCoords)

        for iStation in stationList:
            cubit.cmd(f'import cubit "{pathName}-{str(iStation)}.cub"')
        cubit.cmd(f'save as "{pathName}.cub" overwrite')

        # Remove unnecessary files to save space
        for filePath in glob.glob(f'{pathName}-*.cub'):
            os.remove(filePath)
    return cubit, blade, surfaceDict, birdsMouthVerts, iStationFirstWeb, iStationLastWeb, materialsUsed, spanwiseMatOriCurve


def generateCubitSolidModel(blade, wt_name, settings, crosssectionParams, stationList=None):

    if stationList is None or len(stationList)==0:
        stationList = list(range(len(blade.ispan)))
    elif len(stationList)==1:
        raise ValueError('Need more that one cross section to make a solid model')

    cubit, blade, surfaceDict, birdsMouthVerts, iStationFirstWeb, iStationLastWeb, materialsUsed, spanwiseMatOriCurve = generateCubitCrossSections(
        blade, wt_name, settings, crosssectionParams, '3D', stationList)

    iStationStart = stationList[0]
    iStationEnd = stationList[-1]
    ### ### ### ###
    # Make volumes along the span.
    ### ### ### ###
    meshVolList = []

    partName = 'shell'
    orderedList = getOrderedList(partName)
    if len(orderedList) > 0:
        meshVolList = makeAeroshell(
            surfaceDict, orderedList, meshVolList, iStationEnd)
#     cubit.cmd(f'save as "python2.cub" overwrite')
#     foo

    partName = 'web'
    orderedList = getOrderedList(partName)
    orderedListWeb = orderedList.copy()
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(
            surfaceDict, orderedList, meshVolList, iStationEnd)

    partName = 'roundTEadhesive'
    orderedList = getOrderedList(partName)
    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(
            surfaceDict, orderedList, meshVolList, iStationEnd)

    partName = 'flatTEadhesive'
    orderedList = getOrderedList(partName)

    if orderedList and len(orderedList[0]) > 1:
        meshVolList = makeAeroshell(
            surfaceDict, orderedList, meshVolList, iStationEnd)

    if orderedListWeb and len(orderedListWeb[0]) > 1 and crosssectionParams['amplitudeFraction'] and birdsMouthVerts:

        makeBirdsMouth(blade, birdsMouthVerts,
                       crosssectionParams['amplitudeFraction'], iStationFirstWeb, iStationLastWeb)

    cubit.cmd(f'merge volume {l2s(meshVolList)}')
    cubit.cmd(f'reset volume all')

    cubit.cmd(f'delete surface with Is_Free')
    cubit.cmd(
        f'curve with name "layerThickness*" interval {crosssectionParams["nelPerLayer"]}')
    cubit.cmd('set default autosize on')
    cubit.cmd(f'mesh volume {l2s(meshVolList)}')
    cubit.cmd(f'draw volume {l2s(meshVolList)}')

    # Blocks
    # for imat,material in enumerate(blade.materials):
    for imat, materialName in enumerate(materialsUsed):
        cubit.cmd(f'block {imat+1} add volume with name "*{materialName}*"')
        cubit.cmd(f'block {imat+1} name "{materialName}"')

    addColor(blade, 'volume')

    # Adding Nodesets
    # Root Nodeset
    parseString = f'with name "*station{iStationStart}*"'
    print(f'parseString{parseString}')
    surfaceIDs = parse_cubit_list('surface', parseString)
    cubit.cmd(f'nodeset 1 add surface {l2s(surfaceIDs)} ')
    cubit.cmd(f'nodeset 1 name "root"')

    for iLoop, iStation in enumerate(stationList[1:-1]):
        parseString = f'with name "*station{iStation}*"'
        print(f'parseString{parseString}')
        surfaceIDs = parse_cubit_list('surface', parseString)
        cubit.cmd(f'nodeset {iLoop+2} add surface {l2s(surfaceIDs)} ')
        cubit.cmd(f'nodeset {iLoop+2} name "station{iStation}"')
    if not stationList[1:-1]:
        iLoop = -1
    # Tip Nodeset
    parseString = f'with name "*station{iStationEnd}*"'
    surfaceIDs = parse_cubit_list('surface', parseString)
    cubit.cmd(f'nodeset {iLoop+3} add surface {l2s(surfaceIDs)} ')
    cubit.cmd(f'nodeset {iLoop+3} name "tip"')

    # Outer mold-line nodeset
    cubit.cmd('draw surf with name "*layer0_bottomFace*"')
    parseString = f'with name "*layer0_bottomFace*"'
    surfaceIDs = parse_cubit_list('surface', parseString)
    cubit.cmd(f'nodeset {iLoop+4} add surface {l2s(surfaceIDs)} ')
    cubit.cmd(f'nodeset {iLoop+4} name "oml"')

    # ####################################
    # ### Assign material orientations ###
    # ####################################

    # parseString = f'in volume with name "*volume*"'
    # allVolumeIDs = parse_cubit_list('volume', parseString)

    # for iVol, volumeID in enumerate(allVolumeIDs):
    #     surfIDforMatOri, sign = getMatOriSurface(volumeID, spanwiseMatOriCurve)

    #     for iEl, elementID in enumerate(get_volume_hexes(volumeID)):

    #         coords = cubit.get_center_point("hex", elementID)

    #         cubit.create_vertex(coords[0], coords[1], coords[2])
    #         iVert1 = get_last_id("vertex")
    #         if surfIDforMatOri:
    #             surfaceNormal = sign * \
    #                 np.array(get_surface_normal_at_coord(
    #                     surfIDforMatOri, coords))

    #             curveLocationForTangent = cubit.curve(
    #                 spanwiseMatOriCurve).closest_point(coords)
    #             x = cubit.curve(spanwiseMatOriCurve).tangent(
    #                 curveLocationForTangent)[0]
    #             y = cubit.curve(spanwiseMatOriCurve).tangent(
    #                 curveLocationForTangent)[1]
    #             z = cubit.curve(spanwiseMatOriCurve).tangent(
    #                 curveLocationForTangent)[2]
    #             spanwiseDirection = vectNorm([x, y, z])

    #             perimeterDirection = crossProd(
    #                 surfaceNormal, spanwiseDirection)
    #         else:
    #             perimeterDirection = [1, 0, 0]
    #             surfaceNormal = [0, 1, 0]

    #         length = 0.1
    #         cubit.create_vertex(coords[0]+length*perimeterDirection[0], coords[1] +
    #                             length*perimeterDirection[1], coords[2]+length*perimeterDirection[2])
    #         iVert2 = get_last_id("vertex")
    #         cubit.cmd(f'create curve vertex {iVert1} {iVert2}')
    if settings['export'] is not None:
        if 'g' in settings['export'].lower():
            cubit.cmd(f'export mesh "{wt_name}.g" overwrite')
        if 'cub' in settings['export'].lower():
            cubit.cmd(f'save as "{wt_name}.cub" overwrite')
