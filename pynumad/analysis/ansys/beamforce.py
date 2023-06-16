import numpy as np
import matplotlib.pyplot as plt
import numpy.matlib
from os.path import join
import numpy as np

from pynumad.objects.Subobjects import SkinArea
from pynumad.analysis.ansys.read import *
from pynumad.utils.fatigue import *
from pynumad.utils.interpolation import *
from pynumad.analysis.ansys.write import *


def beamForceToAnsysShell(nodeData, loads, maptype = 'map3D_fxM0'):
    """
    #beamForceToAnsysShell:  Maps AeroDyn forces to an ANSYS blade FE model
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    beamForceToAnsysShell(maptype,nodeData,forcesfile,outfile)
    
    Note: you may omit any or all of the filename arguments and the script
        will provide a file selection box.

    Parameters
    ----------
    nodeData  = node numbers and coordinates for each node
    loads = variables needed to apply nodal forces
        loads["rBlade"]: spanwise location of forces (center of aero element)
        loads["Fxb"]: forces aligned with the X axis of the ANSYS blade model
        loads["Fyb"]: forces aligned with the Y axis of the ANSYS blade model
        loads["Fzb"]: forces aligned with the Z axis of the ANSYS blade model
        loads["Mxb"]: moments about x axis
        loads["Myb"]: moments about y axis
        loads["Mzb"]: moments about z axis
        loads["Alpha"]: CURRENTLY UNUSED
        loads["prebend"]:
        loads["presweep"]:
    maptype = 'map2D_fxM0' Maintains Fx, Fy, and Mz
            fx forces produce zero Mz moment
            'map3D_fxM0' (default) Maintains Fx, Fy, Mz, and root moments
            fx forces produce zero Mz moment
    
    Returns
    -------
    forcemap

    forcesums
    """

    # Transform loads from the FAST blade coordinate system the
    # ANYS Blade Coordinate System
    beta=[[0, -1, 0],   #Direction cosine matrix for a 90 degree clockwise rotation
          [1, 0, 0],  #about the z axis.
          [0, 0, 1]]
    for i in range(len(loads["rBlade"])):
        F = beta @ np.array([loads["Fxb"][i],loads["Fyb"][i],loads["Fzb"][i]])
        M = beta @ np.array([loads["Mxb"][i],loads["Myb"][i],loads["Mzb"][i]])
        #Overwrite loads in the FAST CSYS with the ANSYS CSYS
        loads["Fxb"][i] = F[0]
        loads["Fyb"][i] = F[1]
        loads["Fzb"][i] = F[2]
        loads["Mxb"][i] = M[0]
        loads["Myb"][i] = M[1]
        loads["Mzb"][i] = M[2]
    
    if 'map2D_fxM0' == maptype:
        forcemap = map2D_fxM0(nodeData,loads)
    elif 'map3D_fxM0' == maptype:
        forcemap = map3D_fxM0(nodeData,loads)
    
    forcesums = check_sums(nodeData,loads,forcemap)
    
    # write the ansys commands which apply the calculated forces

    return forcemap, forcesums

    
def beamForceToAnsysShellFollower(nodeData, loads, maptype = 'map3D_fxM0'): 
    #AD2ANSYS  Maps AeroDyn forces to an ANSYS blade FE model
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #   ad2ansys(maptype,nodeData,forcesfile,outfile)
    
    #   Note: you may omit any or all of the filename arguments and the script
    #         will provide a file selection box.
    
    #     maptype = 'map2D_fxM0' Maintains Fx, Fy, and Mz
    #                            fx forces produce zero Mz moment
    #               'map3D_fxM0' (default) Maintains Fx, Fy, Mz, and root moments
    #                            fx forces produce zero Mz moment
    
    #     nodeData  = node numbers and coordinates for each node
    
    #     forcesfile = name of file containing the load definition with
    #     columns:
    #       Z  - spanwise location of forces (center of aero element)
    #       Fx - forces aligned with the X axis of the ANSYS blade model
    #       Fy - forces aligned with the Y axis of the ANSYS blade model
    #       M  - moments about Z axis
    #       Alpha - CURRENTLY UNUSED
    #       x_off -
    #       y_off -
    
    #     outfile = name of the output file to be /INPUT in ANSYS
    
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Transform loads from the FAST blade coordinate system the
    # ANYS Blade Coordinate System
    
    beta=[[0, -1, 0], #Direction cosine matrix for a 90 degree clockwise rotation
          [1, 0, 0],  #about the z axis.
          [0, 0, 1]]
    for i in range(len(loads["rBlade"])):
        F = beta * np.array([[loads["Fxb"](i)],[loads["Fyb"](i)],[loads["Fzb"](i)]])
        M = beta * np.array([[loads["Mxb"](i)],[loads["Myb"](i)],[loads["Mzb"](i)]])
        #Overwrite loads in the FAST CSYS with the ANSYS CSYS
        loads["Fxb"][i] = F[0]
        loads["Fyb"][i] = F[1]
        loads["Fzb"][i] = F[2]
        loads["Mxb"][i] = M[0]
        loads["Myb"][i] = M[1]
        loads["Mzb"][i] = M[2]
    
    #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    # set the default mapping algorithm
    if not ('maptype' is not None)  or len(maptype)==0:
        maptype = 'map3D_fxM0'
    
    #Look through nodeData, remove root nodes
    nodeData_new = []
    for l in range(len(nodeData[:,2])):
        if nodeData[l,2] == 0:
            continue    
            #         nodeData_new(l,:) = [];
        else:
            nodeData_new.append(nodeData[l,:])
    
    nodeData = np.array(nodeData_new)
    if 'map2D_fxM0' == maptype:
        forcemap = map2D_fxM0(nodeData,loads)
    else:
        if 'map3D_fxM0' == maptype:
            forcemap = map3D_fxM0(nodeData,loads)
    
    forcesums = check_sums(nodeData, loads, forcemap)
    
    # write the ansys commands which apply the calculated forces
    # writeforcefile(outfile,forcemap,forcesums,maptype)
    
    return forcemap, forcesums
    
    
def check_sums(nodeData, loads, forcemap):
    forcesums = dict()
    forcesums["Z"] = np.asarray(loads["rBlade"])
    forcesums["Fx"] = np.zeros((len(loads["Fxb"]), 2))
    forcesums["Fx"][:,0] = np.asarray(loads["Fxb"])
    forcesums["Fy"] = np.zeros((len(loads["Fyb"]), 2))
    forcesums["Fy"][:,0] = np.asarray(loads["Fyb"])
    forcesums["M"] = np.zeros((len(loads["Mzb"]), 2))
    forcesums["M"][:,0] = np.asarray(loads["Mzb"])
    ## EMA original:
    # forcesums["RootMx"][:,0] = loads["rBlade"] .* loads["Fyb"];
    # forcesums["RootMy"][:,0] = loads["rBlade"] .* loads["Fxb"];
    ## changed to:
    forcesums["RootMx"] = np.zeros((len(loads["rBlade"]), 2))
    forcesums["RootMx"][:,0] = np.multiply(- np.asarray(loads["rBlade"]),np.asarray(loads["Fyb"])) + np.multiply(loads["prebend"],loads["Fzb"])
    forcesums["RootMy"] = np.zeros((len(loads["rBlade"]), 2))
    forcesums["RootMy"][:,0] = np.multiply(np.asarray(loads["rBlade"]), np.asarray(loads["Fxb"])) - np.multiply(loads["presweep"],loads["Fzb"])
    ## END
    
    for bk in range(len(loads["Fxb"])):
        i = np.where(forcemap["bin"] == bk)
        x = nodeData[i,1] - loads["presweep"][bk]
        y = nodeData[i,2] - loads["prebend"][bk]
        z = nodeData[i,3]
        forcesums["Fx"][bk,1] = sum(forcemap["fx"][i])
        forcesums["Fy"][bk,1] = sum(forcemap["fy"][i])
        forcesums["M"][bk,1] = sum(np.multiply(- y.reshape(-1),forcemap["fx"][i]) + np.multiply(x.reshape(-1),forcemap["fy"][i]))
        ## EMA original:
        #     forcesums["RootMx"](bk,2) = sum(z.*forcemap["fy"][i]);
        #     forcesums["RootMy"](bk,2) = sum(z.*forcemap["fx"][i]);
        ## changed to:
        x = nodeData[i,0]
        y = nodeData[i,1]
        forcesums["RootMx"][bk,0] = sum(np.multiply(- z.reshape(-1),forcemap["fy"][i]) + np.multiply(y.reshape(-1),forcemap["fz"][i]))
        forcesums["RootMy"][bk,0] = sum(np.multiply(z.reshape(-1),forcemap["fx"][i]) - np.multiply(x.reshape(-1),forcemap["fz"][i]))
        ## END
    
    return forcesums
    
 
def map2D_fxM0(nodeData, loads): 
    # Map forces such that the equivalent Fx Fy and M are maintained for each
# section.  It is assumed that only the fy forces contribute to M and fx
# are balanced to produce no moment.
    mesh = dict()
    mesh["n"] = nodeData[:,0]
    mesh["x"] = nodeData[:,1]
    mesh["y"] = nodeData[:,2]
    mesh["z"] = nodeData[:,3]
    # divide nodes into spanwise groups
    bin = np.zeros((mesh["n"].shape,mesh["n"].shape))
    
    for nk in range(len(mesh["n"])):
        halfDZ = loads["rBlade"][1]
        for bk in np.arange(1,np.asarray(loads["rBlade"]).size+1).reshape(-1):
            OBedge = loads["rBlade"][bk] + halfDZ
            if (mesh["z"](nk) <= OBedge) or (bk == np.asarray(loads["rBlade"]).size):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"](bk + 1) - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    forcemap["fy"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    for bk in range(len(loads["Fxb"])):
        i = np.where(bin == bk)
        N = len(i)
        x = mesh["x"][i] - loads["presweep"][bk]
        y = mesh["y"][i] - loads["prebend"][bk]
        mx = np.mean(x)
        my = np.mean(y)
        mxx = np.mean(x ** 2)
        myy = np.mean(y ** 2)
        #mxy = np.mean(x.*y);
        A = np.array([
            [my,1,0,0],
            [0,0,mx,1],
            [0,0,mxx,mx],
            [- myy,- my,0,0]])
        F = 1 / N * np.array([loads["Fxb"][bk],loads["Fyb"][bk],loads["Mzb"][bk],0])
        ab = np.linalg.solve(A,F)
        forcemap["fx"][i] = ab[0] * y + ab[1]
        forcemap["fy"][i] = ab[2] * x + ab[3]
    
    return forcemap
    
    
def map3D_fxM0(nodeData, loads): 
    # Map forces such that the equivalent Fx Fy and M are maintained for each
    # section.  It is assumed that only the fy forces contribute to M and fx
    # are balanced to produce no moment.
    mesh = dict()
    mesh["n"] = nodeData[:,0]
    mesh["x"] = nodeData[:,1]
    mesh["y"] = nodeData[:,2]
    mesh["z"] = nodeData[:,3]
    # aero.Z = forces{1};
    # aero.Fx = forces{2};
    # aero.Fy = forces{3};
    # aero.M = forces{4};
    # aero.alpha = forces{5};
    # aero.xoff = forces{6};
    # aero.yoff = forces{7};
    
    # divide nodes into spanwise groups
    bin = np.zeros(mesh["n"].shape)
    
    for nk in range(len(mesh["n"])):
        halfDZ = loads["rBlade"][0]
        for bk in range(len(loads["rBlade"])):
            OBedge = loads["rBlade"][bk] + halfDZ
            if (mesh["z"][nk] <= OBedge) or (bk == len(loads["rBlade"])-1):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"][bk + 1] - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros(mesh["n"].shape)
    forcemap["fy"] = np.zeros(mesh["n"].shape)
    ## EMA added:
    forcemap["fz"] = np.zeros(mesh["n"].shape)
    ## END
    for bk in range(len(loads["Fxb"])):
        i = np.where(bin == bk)[0]
        N = len(i)
        x = mesh["x"][i] - loads["presweep"][bk]
        y = mesh["y"][i] - loads["prebend"][bk]
        z = mesh["z"][i]
        mx = np.mean(x)
        my = np.mean(y)
        mz = np.mean(z)
        mxx = np.mean(x ** 2)
        myy = np.mean(y ** 2)
        mzz = np.mean(z ** 2)
        #mxy = np.mean(x.*y);
        mzx = np.mean(np.multiply(z,x))
        mzy = np.mean(np.multiply(z,y))
        ## Original
        #     A = [mz,   my,   1,   0,   0,   0;
        #         0,    0,    0,   mz,  mx,  1;
        #         0,    0,    0,   mzx, mxx, mx;
        #         -mzy, -myy, -my,  0,   0,   0;
        #         0,    0,    0,   mzz, mzx, mz;
        #         mzz,  mzy,  mz,  0,   0,   0];
        #     F = 1/N * [aero.Fx[bk]; aero.Fy[bk]; aero.M[bk]; 0; ...
        #         aero.Z[bk]*aero.Fy[bk]; aero.Z[bk]*aero.Fx[bk]];
        #     abc = A\F;
        ## Changed to:
        ## Add/modify equations to include forces in the Z-direction.
        A = np.array([[mz,my,1,0,0,0,0],
                      [0,0,0,mz,mx,1,0],
                      [0,0,0,0,0,0,1],
                      [0,0,0,mzx,mxx,mx,0],
                      [- mzy,- myy,- my,0,0,0,0],
                      [0,0,0,mzz,mzx,mz,- my],
                      [mzz,mzy,mz,0,0,0,- mx]])
        F = 1 / N * np.array([
            loads["Fxb"][bk],
            loads["Fyb"][bk],
            loads["Fzb"][bk],
            loads["Mzb"][bk],
            0,
            loads["rBlade"][bk] * loads["Fyb"][bk],
            loads["rBlade"][bk] * loads["Fxb"][bk]
            ])
        abc = np.linalg.solve(A,F)
        ## END
        forcemap["fx"][i] = abc[0] * z + abc[1] * y + abc[2]
        forcemap["fy"][i] = abc[3] * z + abc[4] * x + abc[5]
        ## EMA added:
        forcemap["fz"][i] = abc[6]
        ## END
    
    return forcemap
    
  