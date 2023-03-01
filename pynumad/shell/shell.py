import numpy as np
import warnings
from os import getcwd
from os.path import join

from pynumad.utils.interpolation import interpolator_wrap
from pynumad.shell.shellClasses import shellRegion, elementSet, NuMesh3D, spatialGridList2D, spatialGridList3D
from pynumad.io.ansys import writeANSYSshellModel


def shellMeshGeneral(blade, forSolid, includeAdhesive):
    """
    This method generates a finite element shell mesh for the blade, based on what is
    stored in blade.geometry, blade.keypoints, and blade.profiles.  Element sets are
    returned corresponding to blade.stacks and blade.swstacks

    Parameters
    -----------

    Returns
    -------
    """

    geomSz = blade.geometry.shape
    lenGeom = geomSz[0]
    numXsec = geomSz[2]
    XSCurvePts = np.array([],dtype=int)

    ## Determine the key curve points along the OML at each cross section
    for i in range(numXsec):
        keyPts = np.array([0])
        minDist = 1
        lePt = 0
        for j in range(lenGeom):
            prof = blade.profiles[j,:,i]
            mag = np.linalg.norm(prof)
            if (mag < minDist):
                minDist = mag
                lePt = j

        for j in range(5):
            kpCrd = blade.keypoints[j,:,i]
            minDist = blade.ichord[i]
            pti = 1
            for k in range(lePt):
                ptCrd = blade.geometry[k,:,i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if (mag < minDist):
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts,[pti]))

        keyPts = np.concatenate((keyPts,[lePt]))
        for j in range(5,10):
            kpCrd = blade.keypoints[j,:,i]
            minDist = blade.ichord[i]
            pti = 1
            for k in range(lePt,lenGeom):
                ptCrd = blade.geometry[k,:,i]
                vec = ptCrd - kpCrd
                mag = np.linalg.norm(vec)
                if (mag < minDist):
                    minDist = mag
                    pti = k
            keyPts = np.concatenate((keyPts,[pti]))

        keyPts = np.concatenate((keyPts,[lenGeom-1]))
        allPts = np.array([keyPts[0]])
        for j in range(0,len(keyPts) - 1):
            secPts = np.linspace(keyPts[j],keyPts[j+1],4)
            secPts = np.round(secPts).astype(int)
            allPts = np.concatenate((allPts,secPts[1:4]))

        XSCurvePts = np.vstack((XSCurvePts,allPts)) if XSCurvePts.size else allPts
    rws,cls = XSCurvePts.shape
    ## Create longitudinal splines down the blade through each of the key X-section points

    splineX = blade.geometry[XSCurvePts[0,:],0,0]
    splineY = blade.geometry[XSCurvePts[0,:],1,0]
    splineZ = blade.geometry[XSCurvePts[0,:],2,0]
    for i in range(1,rws):
        Xrow = blade.geometry[XSCurvePts[i,:],0,i]
        splineX = np.vstack((splineX,Xrow.T))
        Yrow = blade.geometry[XSCurvePts[i,:],1,i]
        splineY = np.vstack((splineY,Yrow.T))
        Zrow = blade.geometry[XSCurvePts[i,:],2,i]
        splineZ = np.vstack((splineZ,Zrow.T))
    
    spParam = np.transpose(np.linspace(0,1,rws))
    nSpi = rws + 2 * (rws - 1)
    spParami = np.transpose(np.linspace(0,1,nSpi))
    splineXi = interpolator_wrap(spParam,splineX[:,0],spParami,'pchip')
    splineYi = interpolator_wrap(spParam,splineY[:,0],spParami,'pchip')
    splineZi = interpolator_wrap(spParam,splineZ[:,0],spParami,'pchip')
    for i in range(1,cls):
        splineXi = np.vstack([splineXi,interpolator_wrap(spParam,splineX[:,i],spParami,'pchip')])
        splineYi = np.vstack([splineYi,interpolator_wrap(spParam,splineY[:,i],spParami,'pchip')])
        splineZi = np.vstack([splineZi,interpolator_wrap(spParam,splineZ[:,i],spParami,'pchip')])
    splineXi = splineXi.T
    splineYi = splineYi.T
    splineZi = splineZi.T
    ## Determine the first spanwise section that needs adhesive
    if (includeAdhesive == 1):
        stPt = 0
        frstXS = 0
        while (frstXS == 0 and stPt < splineXi.shape[0]):

            v1x = splineXi[stPt,6] - splineXi[stPt,4]
            v1y = splineYi[stPt,6] - splineYi[stPt,4]
            v1z = splineZi[stPt,6] - splineZi[stPt,4]
            v2x = splineXi[stPt,30] - splineXi[stPt,32]
            v2y = splineYi[stPt,30] - splineYi[stPt,32]
            v2z = splineZi[stPt,30] - splineZi[stPt,32]
            mag1 = np.sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
            mag2 = np.sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
            dp = (1 / (mag1 * mag2)) * (v1x * v2x + v1y * v2y + v1z * v2z)
            if (dp > 0.7071):
                frstXS = stPt
            stPt = stPt + 3

        if (frstXS == 0):
            frstXS = splineXi.shape[0]
    else:
        frstXS = splineXi.shape[0]
    
    ## Generate the mesh using the splines as surface guides
    nodes = np.array([])
    elements = np.array([])
    ## Outer shell sections
    outerShellElSets = np.array([])
    stPt = 0
    for i in range(rws - 1):
        if (stPt < frstXS):
            setCol = np.array([])
            stSec = 0
            endSec = 11
            stSp = 0
        else:
            setCol = np.array([elementSet(blade.stacks[0,i].name,blade.stacks[0,i].plygroups,[])])
            stSec = 1
            endSec = 10
            stSp = 3
        for j in range(stSec,endSec+1):
            shellKp = np.array([
                [splineXi[stPt,stSp],splineYi[stPt,stSp],splineZi[stPt,stSp]],
                [splineXi[stPt,stSp + 3],splineYi[stPt,stSp + 3],splineZi[stPt,stSp + 3]],
                [splineXi[stPt + 3,stSp + 3],splineYi[stPt + 3,stSp + 3],splineZi[stPt + 3,stSp + 3]],
                [splineXi[stPt + 3,stSp],splineYi[stPt + 3,stSp],splineZi[stPt + 3,stSp]],
                [splineXi[stPt,stSp + 1],splineYi[stPt,stSp + 1],splineZi[stPt,stSp + 1]],
                [splineXi[stPt,stSp + 2],splineYi[stPt,stSp + 2],splineZi[stPt,stSp + 2]],
                [splineXi[stPt + 1,stSp + 3],splineYi[stPt + 1,stSp + 3],splineZi[stPt + 1,stSp + 3]],
                [splineXi[stPt + 2,stSp + 3],splineYi[stPt + 2,stSp + 3],splineZi[stPt + 2,stSp + 3]],
                [splineXi[stPt + 3,stSp + 2],splineYi[stPt + 3,stSp + 2],splineZi[stPt + 3,stSp + 2]],
                [splineXi[stPt + 3,stSp + 1],splineYi[stPt + 3,stSp + 1],splineZi[stPt + 3,stSp + 1]],
                [splineXi[stPt + 2,stSp],splineYi[stPt + 2,stSp],splineZi[stPt + 2,stSp]],
                [splineXi[stPt + 1,stSp],splineYi[stPt + 1,stSp],splineZi[stPt + 1,stSp]],
                [splineXi[stPt + 1,stSp + 1],splineYi[stPt + 1,stSp + 1],splineZi[stPt + 1,stSp + 1]],
                [splineXi[stPt + 1,stSp + 2],splineYi[stPt + 1,stSp + 2],splineZi[stPt + 1,stSp + 2]],
                [splineXi[stPt + 2,stSp + 2],splineYi[stPt + 2,stSp + 2],splineZi[stPt + 2,stSp + 2]],
                [splineXi[stPt + 2,stSp + 1],splineYi[stPt + 2,stSp + 1],splineZi[stPt + 2,stSp + 1]]
            ])
            vec = shellKp[1,:] - shellKp[0,:]
            mag = np.linalg.norm(vec)
            nEl = np.array([],dtype=int)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[2,:] - shellKp[1,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[3,:] - shellKp[2,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[0,:] - shellKp[3,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            shell = shellRegion('quad16',shellKp,nEl)
            regNodes,regElements = shell.createShellMesh('quad','structured')
            numNds = nodes.shape[0]
            numEls = elements.shape[0]
            for k in range(regElements.shape[0]):
                for m in range(regElements.shape[1]):
                    if (regElements[k,m] < 0):
                        regElements[k,m] = - numNds
            regElements = regElements + numNds

            nodes = np.vstack([nodes,regNodes]) if nodes.size else regNodes
            elements = np.vstack([elements,regElements]) if elements.size else regElements
            elList = np.array(np.arange(numEls + 1,elements.shape[0]+1))
            newSet = elementSet(blade.stacks[j,i].name,blade.stacks[j,i].plygroups,elList)
            setCol = np.concatenate([setCol,[newSet]])
            stSp = stSp + 3
        if (stPt >= frstXS):
            newSet = elementSet(blade.stacks[11,i].name,blade.stacks[11,i].plygroups,[])
            setCol = np.concatenate([setCol,[newSet]])
        outerShellElSets = np.concatenate([outerShellElSets,setCol])
        stPt = stPt + 3
    
    ## Shift the appropriate splines if the mesh is for a solid model seed
    if (forSolid == 1):
        caseIndex = np.array([
            [9,27,3],
            [12,24,3],
            [24,12,8],
            [27,9,8]])
        # 5,33,2; ...
        # 6,32,2; ...
        # 7,31,2; ...
        # 33,5,11; ...
        # 32,6,11; ...
        # 31,7,11];
        for i in range(caseIndex.shape[0]):
            spl = caseIndexc[i,0]
            tgtSp = caseIndex[i,1]
            sec = caseIndex[i,2]
            stPt = 0
            for j in range(rws-1):
                totalThick = 0
                for k in range(3):
                    tpp = 0.001 * blade.stacks[sec,j].plygroups[k].thickness
                    npls = blade.stacks[sec,j].plygroups[k].nPlies
                    totalThick = totalThick + tpp * npls
                for k in range(3):
                    vx = splineXi[stPt,tgtSp] - splineXi[stPt,spl]
                    vy = splineYi[stPt,tgtSp] - splineYi[stPt,spl]
                    vz = splineZi[stPt,tgtSp] - splineZi[stPt,spl]
                    magInv = 1 / np.sqrt(vx * vx + vy * vy + vz * vz)
                    ux = magInv * vx
                    uy = magInv * vy
                    uz = magInv * vz
                    splineXi[stPt,spl] = splineXi[stPt,spl] + totalThick * ux
                    splineYi[stPt,spl] = splineYi[stPt,spl] + totalThick * uy
                    splineZi[stPt,spl] = splineZi[stPt,spl] + totalThick * uz
                    stPt = stPt + 1
    
    ## Shear web sections
    stPt = 0
    web1Sets = np.array([])
    web2Sets = np.array([])
    for i in range(rws-1):
        if blade.swstacks[0][i].plygroups:
            shellKp = np.zeros((16,3))
            shellKp[0,:] = np.array([splineXi[stPt,12],splineYi[stPt,12],splineZi[stPt,12]])
            shellKp[1,:] = np.array([splineXi[stPt,24],splineYi[stPt,24],splineZi[stPt,24]])
            shellKp[2,:] = np.array([splineXi[stPt + 3,24],splineYi[stPt + 3,24],splineZi[stPt + 3,24]])
            shellKp[3,:] = np.array([splineXi[stPt + 3,12],splineYi[stPt + 3,12],splineZi[stPt + 3,12]])
            shellKp[6,:] = np.array([splineXi[stPt + 1,24],splineYi[stPt + 1,24],splineZi[stPt + 1,24]])
            shellKp[7,:] = np.array([splineXi[stPt + 2,24],splineYi[stPt + 2,24],splineZi[stPt + 2,24]])
            shellKp[10,:] = np.array([splineXi[stPt + 2,12],splineYi[stPt + 2,12],splineZi[stPt + 2,12]])
            shellKp[11,:] = np.array([splineXi[stPt + 1,12],splineYi[stPt + 1,12],splineZi[stPt + 1,12]])
            shellKp[4,:] = 0.6666 * shellKp[0,:] + 0.3333 * shellKp[1,:]
            shellKp[5,:] = 0.3333 * shellKp[0,:] + 0.6666 * shellKp[1,:]
            shellKp[8,:] = 0.6666 * shellKp[2,:] + 0.3333 * shellKp[3,:]
            shellKp[9,:] = 0.3333 * shellKp[2,:] + 0.6666 * shellKp[3,:]
            shellKp[12,:] = 0.6666 * shellKp[11,:] + 0.3333 * shellKp[6,:]
            shellKp[13,:] = 0.3333 * shellKp[11,:] + 0.6666 * shellKp[6,:]
            shellKp[14,:] = 0.6666 * shellKp[7,:] + 0.3333 * shellKp[10,:]
            shellKp[15,:] = 0.3333 * shellKp[7,:] + 0.6666 * shellKp[10,:]

            vec = shellKp[1,:] - shellKp[0,:]
            mag = np.linalg.norm(vec)

            nEl = np.array([],dtype=int)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[2,:] - shellKp[1,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[3,:] - shellKp[2,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[0,:] - shellKp[3,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])

            shell = shellRegion('quad16',shellKp,nEl)
            regNodes,regElements = shell.createShellMesh('quad','structured')
            numNds = nodes.shape[0]
            numEls = elements.shape[0]
            for k in range(regElements.shape[0]):
                for m in range(regElements.shape[1]):
                    if (regElements[k,m] < 1):
                        regElements[k,m] = - numNds
            regElements = regElements + numNds
            if nodes.size ==0:
                nodes = regNodes
            else:
                nodes = np.vstack([nodes,regNodes])

            if elements.size ==0:
                elements = regElements
            else:
                elements = np.vstack([elements,regElements])
            elList = np.array([np.arange(numEls,elements.shape[0])])
            newSet = elementSet(blade.swstacks[0][i].name,blade.swstacks[0][i].plygroups,elList)
            web1Sets = np.concatenate([web1Sets,[newSet]])
        else:
            newSet = elementSet(blade.swstacks[0][i].name,blade.swstacks[0][i].plygroups,[])
            web1Sets = np.concatenate([web1Sets,[newSet]])
        if blade.swstacks[1][i].plygroups:
            shellKp = np.zeros((16,3))
            shellKp[0,:] = np.array([splineXi[stPt,27],splineYi[stPt,27],splineZi[stPt,27]])
            shellKp[1,:] = np.array([splineXi[stPt,9],splineYi[stPt,9],splineZi[stPt,9]])
            shellKp[2,:] = np.array([splineXi[stPt + 3,9],splineYi[stPt + 3,9],splineZi[stPt + 3,9]])
            shellKp[3,:] = np.array([splineXi[stPt + 3,27],splineYi[stPt + 3,27],splineZi[stPt + 3,27]])
            shellKp[6,:] = np.array([splineXi[stPt + 1,9],splineYi[stPt + 1,9],splineZi[stPt + 1,9]])
            shellKp[7,:] = np.array([splineXi[stPt + 2,9],splineYi[stPt + 2,9],splineZi[stPt + 2,9]])
            shellKp[10,:] = np.array([splineXi[stPt + 2,27],splineYi[stPt + 2,27],splineZi[stPt + 2,27]])
            shellKp[11,:] = np.array([splineXi[stPt + 1,27],splineYi[stPt + 1,27],splineZi[stPt + 1,27]])
            shellKp[4,:] = 0.6666 * shellKp[0,:] + 0.3333 * shellKp[1,:]
            shellKp[5,:] = 0.3333 * shellKp[0,:] + 0.6666 * shellKp[1,:]
            shellKp[8,:] = 0.6666 * shellKp[2,:] + 0.3333 * shellKp[3,:]
            shellKp[9,:] = 0.3333 * shellKp[2,:] + 0.6666 * shellKp[3,:]
            shellKp[12,:] = 0.6666 * shellKp[11,:] + 0.3333 * shellKp[6,:]
            shellKp[13,:] = 0.3333 * shellKp[11,:] + 0.6666 * shellKp[6,:]
            shellKp[14,:] = 0.6666 * shellKp[7,:] + 0.3333 * shellKp[10,:]
            shellKp[15,:] = 0.3333 * shellKp[7,:] + 0.6666 * shellKp[10,:]

            vec = shellKp[1,:] - shellKp[0,:]
            mag = np.linalg.norm(vec)

            nEl = np.array([],dtype=int)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[2,:] - shellKp[1,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[3,:] - shellKp[2,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])
            vec = shellKp[0,:] - shellKp[3,:]
            mag = np.linalg.norm(vec)
            nEl = np.concatenate([nEl,[np.ceil(mag / blade.mesh).astype(int)]])

            shell = shellRegion('quad16',shellKp,nEl)
            regNodes,regElements = shell.createShellMesh('quad','structured')
            numNds = nodes.shape[0]
            numEls = elements.shape[0]
            for k in range(regElements.shape[0]):
                for m in range(regElements.shape[1]):
                    if (regElements[k,m] < 1):
                        regElements[k,m] = - numNds
            regElements = regElements + numNds
            if nodes.size ==0:
                nodes = regNodes
            else:
                nodes = np.vstack([nodes,regNodes])

            if elements.size ==0:
                elements = regElements
            else:
                elements = np.vstack([elements,regElements])
            elList = np.array([range(numEls+1,elements.shape[0])])
            newSet = elementSet(blade.swstacks[0][i].name,blade.swstacks[1][i].plygroups,elList)
            web2Sets = np.concatenate([web2Sets,[newSet]])
        else:
            newSet = elementSet(blade.swstacks[0][i].name,blade.swstacks[1][i].plygroups,[])
            web2Sets = np.concatenate([web2Sets,[newSet]])
        stPt = stPt + 3
    
    # plotShellMesh(nodes,elements);
    # keyboard
    shearWebElSets = [web1Sets, web2Sets]
    ## Eliminate duplicate nodes from the global list and update element connectivity
    minX = np.amin(nodes[:,0]) - blade.mesh
    maxX = np.amax(nodes[:,0]) + blade.mesh
    minY = np.amin(nodes[:,1]) - blade.mesh
    maxY = np.amax(nodes[:,1]) + blade.mesh
    minZ = np.amin(nodes[:,2]) - blade.mesh
    maxZ = np.amax(nodes[:,2]) + blade.mesh
    gS = 2 * blade.mesh
    nodeGL = spatialGridList3D(minX,maxX,minY,maxY,minZ,maxZ,gS,gS,gS)
    numNds = len(nodes)
    for i in range(numNds):
        nodeGL = nodeGL.addEntry(i,nodes[i,:])
    
    nodeElim = np.zeros((numNds),dtype=int)
    newLabel = np.zeros((numNds),dtype=int)
    newNodes = np.array([])
    lab = 0
    for i in range(numNds):
        if (nodeElim[i] == 0):
            lab = lab + 1
            newLabel[i] = lab
            try:
                newNodes = np.vstack([newNodes,nodes[i,:]])
            except ValueError:
                newNodes = nodes[i,:]
            nearNds = nodeGL.findInRadius(nodes[i,:],blade.mesh)
            for j in range(len(nearNds)):
                k = nearNds[j]
                if (k > i and nodeElim[k] == 0):
                    k = nearNds[j]
                    vec = nodes[i,:] - nodes[k,:]
                    dist = np.linalg.norm(vec)
                    if (dist < blade.mesh * 1e-10):
                        nodeElim[k] = i
    
    nodes = newNodes
    numEl,elNds = elements.shape
    for i in range(numEl):
        for j in range(elNds):
            if (elements[i,j] != 0):
                orig = elements[i,j]
                k = nodeElim[orig]
                if (k != 0):
                    elements[i,j] = newLabel[k]
                else:
                    elements[i,j] = newLabel[orig]
    
    ## Generate mesh for trailing edge adhesive if requested
    if (includeAdhesive == 1):
        stPt = frstXS
        v1x = splineXi[stPt,6] - splineXi[stPt,4]
        v1y = splineYi[stPt,6] - splineYi[stPt,4]
        v1z = splineZi[stPt,6] - splineZi[stPt,4]
        mag1 = np.sqrt(v1x * v1x + v1y * v1y + v1z * v1z)
        v2x = splineXi[stPt,30] - splineXi[stPt,32]
        v2y = splineYi[stPt,30] - splineYi[stPt,32]
        v2z = splineZi[stPt,30] - splineZi[stPt,32]
        mag2 = np.sqrt(v2x * v2x + v2y * v2y + v2z * v2z)
        v3x = splineXi[stPt,6] - splineXi[stPt,30]
        v3y = splineYi[stPt,6] - splineYi[stPt,30]
        v3z = splineZi[stPt,6] - splineZi[stPt,30]
        mag3 = np.sqrt(v3x * v3x + v3y * v3y + v3z * v3z)
        v4x = splineXi[stPt,4] - splineXi[stPt,32]
        v4y = splineYi[stPt,4] - splineYi[stPt,32]
        v4z = splineZi[stPt,4] - splineZi[stPt,32]
        mag4 = np.sqrt(v4x * v4x + v4y * v4y + v4z * v4z)
        nE1 = np.ceil(mag1 / blade.mesh).astype(int)
        nE2 = np.ceil(mag3 / blade.mesh).astype(int)
        nE3 = np.ceil(mag2 / blade.mesh).astype(int)
        nE4 = np.ceil(mag4 / blade.mesh).astype(int)
        nEl = np.array([nE1,nE2,nE3,nE4])
        gdLayer = 0
        sweepElements = []
        while (stPt < splineXi.shape[0]):

            #                     shellKp = zeros(16,3);
#                     shellKp[0,:] = [splineXi(stPt,4),splineYi(stPt,4),splineZi(stPt,4)];
#                     shellKp[1,:] = [splineXi(stPt,7),splineYi(stPt,7),splineZi(stPt,7)];
#                     shellKp[2,:] = [splineXi(stPt,31),splineYi(stPt,31),splineZi(stPt,31)];
#                     shellKp[3,:] = [splineXi(stPt,34),splineYi(stPt,34),splineZi(stPt,34)];
#                     shellKp(5,:) = [splineXi(stPt,5),splineYi(stPt,5),splineZi(stPt,5)];
#                     shellKp(6,:) = [splineXi(stPt,6),splineYi(stPt,6),splineZi(stPt,6)];
#                     shellKp[6,:] = 0.6666*shellKp[1,:] + 0.3333*shellKp[2,:];
#                     shellKp[7,:] = 0.3333*shellKp[1,:] + 0.6666*shellKp[2,:];
#                     shellKp(9,:) = [splineXi(stPt,32),splineYi(stPt,32),splineZi(stPt,32)];
#                     shellKp(10,:) = [splineXi(stPt,33),splineYi(stPt,33),splineZi(stPt,33)];
#                     shellKp[10,:] = 0.3333*shellKp[0,:] + 0.6666*shellKp[3,:];
#                     shellKp[11,:] = 0.6666*shellKp[0,:] + 0.3333*shellKp[3,:];
#                     shellKp(13,:) = 0.6666*shellKp(5,:) + 0.3333*shellKp(10,:);
#                     shellKp(14,:) = 0.6666*shellKp(6,:) + 0.3333*shellKp(9,:);
#                     shellKp(15,:) = 0.3333*shellKp(6,:) + 0.6666*shellKp(9,:);
#                     shellKp(16,:) = 0.3333*shellKp(5,:) + 0.6666*shellKp(10,:);
#                     shell = shellRegion('quad16',shellKp,nEl);
#                     [bNodes,bFaces] = shell.createShellMesh('quad','free');
            shellKp = np.zeros((9,3))
            shellKp[0,:] = np.array([splineXi[stPt,4],splineYi[stPt,4],splineZi[stPt,4]])
            shellKp[1,:] = np.array([splineXi[stPt,6],splineYi[stPt,6],splineZi[stPt,6]])
            shellKp[2,:] = np.array([splineXi[stPt,30],splineYi[stPt,30],splineZi[stPt,30]])
            shellKp[3,:] = np.array([splineXi[stPt,32],splineYi[stPt,32],splineZi[stPt,32]])
            shellKp[4,:] = np.array([splineXi[stPt,5],splineYi[stPt,5],splineZi[stPt,5]])
            shellKp[5,:] = 0.5 * shellKp[1,:] + 0.5 * shellKp[2,:]
            shellKp[6,:] = np.array([splineXi[stPt,31],splineYi[stPt,31],splineZi[stPt,31]])
            shellKp[7,:] = 0.5 * shellKp[0,:] + 0.5 * shellKp[3,:]
            shellKp[8,:] = 0.5 * shellKp[4,:] + 0.5 * shellKp[6,:]
            shell = shellRegion('quad9',shellKp,nEl)
            bNodes,bFaces = shell.createShellMesh('quad','free')
            #                     plotShellMesh(bNodes,bFaces);
#                     keyboard
            if (stPt == frstXS):
                adhesMesh = NuMesh3D(bNodes,bFaces)
            else:
                gdLayer = gdLayer + 1
                guideNds[gdLayer] = bNodes
                layerSwEl = np.ceil((splineZi(stPt,4) - splineZi(stPt - 3,4)) / blade.mesh).astype(int)
                sweepElements = np.array([sweepElements,layerSwEl])
            stPt = stPt + 3

        #                 sweepElements = ceil((splineZi(end,4) - splineZi(frstXS,4))/blade.mesh);
        adhesNds,adhesEls = adhesMesh.createSweptMesh('to_dest_nodes',[],[],sweepElements,[],guideNds)
    else:
        adhesNds = []
        adhesEls = []
    
    return nodes,elements,outerShellElSets,shearWebElSets,adhesNds,adhesEls
  
    
def generateShellModel(blade, feaCode, includeAdhesive, varargin): 
    # This method generates a shell FEA model in one of the supported FEA codes; w/ or w/o adhesieve
    
    if str(feaCode.lower()) == str('ansys'):
        global ansysPath
        # define ANSYS model settings (can be options in generateFEA)
        config = {}
        config["BoundaryCondition"] = 'cantilevered'
        config["elementType"] = '181'
        config["MultipleLayerBehavior"] = 'multiply'
        config["dbgen"] = 1
        config["dbname"] = 'master'
        # Generate a mesh using shell elements
        APDLname = 'buildAnsysShell.src'
        ansys_product = 'ANSYS'
        blade.paths.job = getcwd()
        filename = join(blade.paths.job,APDLname)
        if len(varargin)==0:
            forSolid = 0
            meshData.nodes,meshData.elements,meshData.outerShellElSets,meshData.shearWebElSets,meshData.adhesNds,meshData.adhesEls = blade.shellMeshGeneral(forSolid,includeAdhesive)
        else:
            meshData = varargin[0]
        writeANSYSshellModel(blade,filename,meshData,config,includeAdhesive)
        if config["dbgen"]:
            if len(ansysPath)==0:
                errordlg('Path to ANSYS not specified. Aborting.','Operation Not Permitted')
                return meshData
            try:
                #tcl: exec "$ANSYS_path" -b -p $AnsysProductVariable -I shell7.src -o output.txt
                ansys_call = sprintf('"%s" -b -p %s -I %s -o output.txt',ansysPath,ansys_product,APDLname)
                status,result = dos(ansys_call)
                if status==0:
                    # dos command completed successfully; log written to output.txt
                    if 1:
                        print('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.')
                    else:
                        helpdlg('ANSYS batch run to generate database (.db) has completed. See "output.txt" for any warnings.','ANSYS Call Completed')
                if status==7:
                    # an error has occured which is stored in output.txt
                    if 1:
                        print('Could not complete ANSYS call. See "output.txt" for details.')
                    else:
                        warndlg('Could not complete ANSYS call. See "output.txt" for details.','Error: ANSYS Call')
            finally:
                pass
    else:
        raise Exception('FEA code "%s" not supported.',feaCode)
    
    return meshData
    
    
def getSolidMesh(blade, layerNumEls): 
    ## Edit stacks to be usable for 3D solid mesh
    blade.editStacksForSolidMesh()
    ## Create shell mesh as seed
    shNodes,shElements,outerShellElSets,shearWebElSets,adhesNds,adhesEls = blade.shellMeshGeneral(1,1)
    ## Initialize 3D solid mesh from the shell mesh
    bladeMesh = NuMesh3D(shNodes,shElements)
    ## Calculate unit normal vectors for all nodes
    numNds = shNodes.shape[1-1]
    nodeNorms = np.zeros((numNds,3))
    for i in np.arange(1,shElements.shape[1-1]+1).reshape(-1):
        n1 = shElements[i,0]
        n2 = shElements[i,1]
        n3 = shElements[i,2]
        n4 = shElements[i,3]
        if (n4 == 0):
            v1 = shNodes[n3,:] - shNodes[n1,:]
            v2 = shNodes[n2,:] - shNodes[n1,:]
        else:
            v1 = shNodes[n4,:] - shNodes[n2,:]
            v2 = shNodes[n3,:] - shNodes[n1,:]
        v3x = v1[1] * v2[2] - v1[2] * v2[1]
        v3y = v1[2] * v2[0] - v1[0] * v2[2]
        v3z = v1[0] * v2[1] - v1[1] * v2[0]
        v3 = np.array([v3x,v3y,v3z])
        mag = np.sqrt(v3 * v3.T)
        uNorm = (1 / mag) * v3
        for j in range(4):
            nj = shElements[i,j]
            if (nj != 0):
                nodeNorms[nj,:] = nodeNorms[nj,:] + uNorm
    
    for i in range(numNds):
        mag = np.sqrt(nodeNorms[i,:] * np.transpose(nodeNorms[i,:]))
        nodeNorms[i,:] = (1 / mag) * nodeNorms[i,:]
    
    ## Extrude shell mesh into solid mesh
    if (len(layerNumEls)==0):
        layerNumEls = np.array([1,1,1])
    
    for i in np.arange(1,len(layerNumEls)+1).reshape(-1):
        nodeDist = np.zeros((numNds,1))
        nodeHitCt = np.zeros((numNds,1))
        numSec,numStat = outerShellElSets.shape
        for j in np.arange(1,numSec+1).reshape(-1):
            for k in np.arange(1,numStat+1).reshape(-1):
                eSet = outerShellElSets(j,k)
                projDist = 0
                for m in np.arange(1,i+1).reshape(-1):
                    tpp = 0.001 * eSet.plygroups(m).thickness
                    npls = eSet.plygroups(m).nPlies
                    projDist = projDist + tpp * npls
                eLst = eSet.elementList
                for m in range(len(eLst)):
                    el = eLst(m)
                    for p in range(4):
                        nd = shElements(el,p)
                        if (nd != 0):
                            if (j == 3 or j == 8):
                                nodeDist[nd] = projDist
                                nodeHitCt[nd] = - 1
                            else:
                                if (nodeHitCt(nd) != - 1):
                                    nodeDist[nd] = nodeDist(nd) + projDist
                                    nodeHitCt[nd] = nodeHitCt(nd) + 1
        for j in range(2):
            setLst = shearWebElSets[j]
            for k in range(len(setLst)):
                eSet = setLst(k)
                if (not len(eSet.elementList)==0 ):
                    projDist = 0
                    for m in range(i):
                        tpp = 0.001 * eSet.plygroups(m).thickness
                        npls = eSet.plygroups(m).nPlies
                        projDist = projDist + tpp * npls
                    eLst = eSet.elementList
                    for m in range(len(eLst)):
                        el = eLst(m)
                        for p in range(4):
                            nd = shElements(el,p)
                            if (nd != 0):
                                nodeDist[nd] = nodeDist[nd] + projDist
                                nodeHitCt[nd] = nodeHitCt[nd] + 1
        newLayer = shNodes
        for j in np.arange(numNds):
            if (nodeHitCt(j) != - 1):
                nodeDist[j] = (1 / nodeHitCt(j)) * nodeDist(j)
            newLayer[j,:] = newLayer[j,:] + nodeDist[j] * nodeNorms[j,:]
        guideNds[i] = newLayer
    
    nodes,elements = bladeMesh.createSweptMesh('to_dest_nodes',[],[],layerNumEls,[],guideNds)
    numShEls = shElements.shape[0]
    totLayers = np.sum(layerNumEls, 'all'-1)
    for i in range(numSec):
        for j in range(numStat):
            eSet = outerShellElSets[i,j]
            elst = eSet.elementList
            newELst = elst
            if (not len(elst)==0 ):
                for k in range(totLayers-1):
                    newELst = np.array([[newELst],[(elst + k * numShEls)]])
            outerShellElSets[i,j].elementList = newELst
    
    for i in range(2):
        setLst = shearWebElSets[i]
        for j in range(len(setLst)):
            eSet = setLst(j)
            elst = eSet.elementList
            newELst = elst
            if (not len(elst)==0 ):
                for k in range(totLayers-1):
                    newELst = np.array([[newELst],[(elst + k * numShEls)]])
            shearWebElSets[i][j].elementList = newELst
    
    numSolidNds = nodes.shape[0]
    for i in range(adhesEls.shape[0]):
        for j in range(8):
            k = adhesEls[i,j]
            if (k != 0):
                adhesEls[i,j] = k + numSolidNds
    
    nodes = np.array([nodes,adhesNds])
    elements = np.array([[elements],[adhesEls]])
    lastEl = elements.shape[0]
    stAdEl = lastEl - adhesEls.shape[0] + 1
    eList = np.array([range(stAdEl,lastEl)])
    adhesiveElSet = elementSet('adhesive',[],eList)
    return nodes,elements,outerShellElSets,shearWebElSets,adhesiveElSet
    
    return nodes,elements,outerShellElSets,shearWebElSets,adhesNds,adhesEls

