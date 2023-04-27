import numpy as np
import matplotlib.pyplot as plt
import numpy.matlib
from os.path import join
import numpy as np

from pynumad.objects.Subobjects import SkinArea
from pynumad.analysis.ansys.read import *
from pynumad.utils.fatigue import *
from pynumad.utils.interpolation import *

    
def beamForceToAnsysShell(maptype, nodeData, loads, outfile): 
    #beamForceToAnsysShell:  Maps AeroDyn forces to an ANSYS blade FE model
    # **********************************************************************
    # *                   Part of the SNL NuMAD Toolbox                    *
    # * Developed by Sandia National Laboratories Wind Energy Technologies *
    # *             See license.txt for disclaimer information             *
    # **********************************************************************
    #   beamForceToAnsysShell(maptype,nodeData,forcesfile,outfile)
    
    #   Note: you may omit any or all of the filename arguments and the script
    #         will provide a file selection box.
    
    #     maptype = 'map2D_fxM0' Maintains Fx, Fy, and Mz
    #                            fx forces produce zero Mz moment
    #               'map3D_fxM0' (default) Maintains Fx, Fy, Mz, and root moments
    #                            fx forces produce zero Mz moment
    
    #     nodeData  = node numbers and coordinates for each node
    
    #     loads = variables needed to apply nodal forces
    #          loads["rBlade"]: spanwise location of forces (center of aero element)
    #          loads["Fxb"]: forces aligned with the X axis of the ANSYS blade model
    #          loads["Fyb"]: forces aligned with the Y axis of the ANSYS blade model
    #          loads["Fzb"]: forces aligned with the Z axis of the ANSYS blade model
    #          loads["Mxb"]: moments about x axis
    #          loads["Myb"]: moments about y axis
    #          loads["Mzb"]: moments about z axis
    #          loads["Alpha"]: CURRENTLY UNUSED
    #          loads["prebend"]:
    #          loads["presweep"]:
    #     outfile = name of the output file to be /INPUT in ANSYS
    
    #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Transform loads from the FAST blade coordinate system the
    # ANYS Blade Coordinate System
    
    beta=[[0, -1, 0],   #Direction cosine matrix for a 90 degree clockwise rotation
          [1, 0, 0],  #about the z axis.
          [0, 0, 1]]
    for i in range(len(loads["rBlade"])):
        F = beta * np.array([[loads["Fxb"][i]],[loads["Fyb"][i]],[loads["Fzb"][i]]])
        M = beta * np.array([[loads["Mxb"][i]],[loads["Myb"][i]],[loads["Mzb"][i]]])
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
    
    if 'map2D_fxM0' == maptype:
        forcemap = map2D_fxM0(nodeData,loads)
    else:
        if 'map3D_fxM0' == maptype:
            forcemap = map3D_fxM0(nodeData,loads)
    
    forcesums = check_sums(nodeData,loads,forcemap)
    # open file selector if 'outfile' not given
    if not ('outfile' is not None)  or len(outfile)==0:
        outfile = join(nlistpath,'forces.src')
        fn,pn = uiputfile(np.array(['*.*','All files (*.*)']),'Save Output As',outfile)
        if fn==0:
            # the user canceled file selection
            return
        outfile = join(pn,fn)
    
    # write the ansys commands which apply the calculated forces
    writeforcefile(outfile,forcemap,forcesums,maptype)
    # plot the results
    if 1:
        plt.figure(1007)
        subplot(2,2,np.array([1,3]))
        ds = np.ceil(nodeData.shape[0] / 1000.0)
        k = np.arange(1,nodeData.shape[0]+ds,ds)
        quiver3(nodeData[k,1],nodeData[k,2],nodeData[k,3],forcemap["fx"][k],forcemap["fy"][k],np.zeros((forcemap["fx"][k].shape,forcemap["fx"][k].shape)))
        plt.title('Force vectors visual check (downsampled)')
        plt.axis('equal')
        subplot(2,2,2)
        __,k = __builtint__.sorted(nodeData[:,3])
        plt.plot(nodeData[k,3],forcemap["fx"][k],'.')
        plt.xlabel('z')
        plt.ylabel('fx')
        subplot(2,2,4)
        plt.plot(nodeData[k,3],forcemap["fy"][k],'.')
        plt.xlabel('z')
        plt.ylabel('fy')
    
    return
    
    
def read_forces(filename): 
    # Open the file
    fid = open(filename)
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    header = fid.readline()
    
    filecontents = np.transpose(fid.read())
    fid.close()
    # 'Z (m)	Fx (N)	Fy (N)	M (N-m)	Alpha	x_off	y_off'
    forces = textscan(filecontents,np.matlib.repmat('%f',1,7))
    return forces
    
    
def map2D_fxM0(nodeData, loads):
    """
    Map forces such that the equivalent Fx Fy and M are maintained for each
    section.  It is assumed that only the fy forces contribute to M and fx
    are balanced to produce no moment.
    """
    mesh = dict()
    mesh["n"] = nodeData[:,0]
    mesh["x"] = nodeData[:,1]
    mesh["y"] = nodeData[:,2]
    mesh["z"] = nodeData[:,3]
    # divide nodes into spanwise groups
    bin = np.zeros((mesh["n"].shape,mesh["n"].shape))
    
    for nk in range(len(mesh["n"])):
        halfDZ = loads["rBlade"][1]
        for bk in range(loads["rBlade"].size):
            OBedge = loads["rBlade"][bk] + halfDZ
            if (mesh["z"][nk] <= OBedge) or (bk == loads["rBlade"].size):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"][bk + 1] - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    forcemap["fy"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    for bk in np.arange(1,len(loads["Fxb"])+1).reshape(-1):
        i = np.where(bin == bk)
        N = len(i)
        x = mesh["x"][i] - loads["presweep"][bk]
        y = mesh["y"][i] - loads["prebend"][bk]
        mx = np.mean[x]
        my = np.mean[y]
        mxx = np.mean(x ** 2)
        myy = np.mean(y ** 2)
        #mxy = np.mean(x.*y);
        A = np.array([[my,1,0,0],[0,0,mx,1],[0,0,mxx,mx],[- myy,- my,0,0]])
        F = 1 / N * np.array([[loads["Fxb"][bk]],[loads["Fyb"][bk]],[loads["Mzb"][bk]],[0]])
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
    # divide nodes into spanwise groups
    bin = np.zeros((mesh["n"].shape,mesh["n"].shape))
    
    for nk in range(len(mesh["n"])):
        halfDZ = loads["rBlade"][1]
        for bk in range(np.asarray(loads["rBlade"]).size):
            OBedge = loads["rBlade"][bk] + halfDZ
            if (mesh["z"][nk] <= OBedge) or (bk == np.asarray(loads["rBlade"]).size):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"][bk+1] - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    forcemap["fy"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    ## EMA added:
    forcemap["fz"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    ## END
    for bk in np.arange(1,len(loads["Fxb"])+1).reshape(-1):
        i = np.where(bin == bk)
        N = len(i)
        x = mesh["x"][i] - loads["presweep"][bk]
        y = mesh["y"][i] - loads["prebend"][bk]
        z = mesh["z"][i]
        mx = np.mean[x]
        my = np.mean[y]
        mz = np.mean[z]
        mxx = np.mean(x**2)
        myy = np.mean(y**2)
        mzz = np.mean(z**2)
        #mxy = np.mean(x.*y);
        mzx = np.mean(np.multiply(z,x))
        mzy = np.mean(np.multiply(z,y))
        ## EMA original:
        #     A = [mz,   my,   1,   0,   0,   0;
        #        0,    0,    0,   mz,  mx,  1;
        #        0,    0,    0,   mzx, mxx, mx;
        #        -mzy, -myy, -my,  0,   0,   0;
        #        0,    0,    0,   mzz, mzx, mz;
        #        mzz,  mzy,  mz,  0,   0,   0];
        #     F = 1/N * [loads["Fxb"][bk]; loads["Fyb"][bk]; loads["Mzb"][bk]; 0; ...
        #        loads["rBlade"][bk]*loads["Fyb"][bk]; loads["rBlade"][bk]*loads["Fxb"][bk]];
        #     abc = A\F;
        ## Changed to:
        ## Add/modify equations to include forces in the Z-direction.
        A = np.array([[mz,my,1,0,0,0,0],[0,0,0,mz,mx,1,0],[0,0,0,0,0,0,1],[0,0,0,mzx,mxx,mx,0],[- mzy,- myy,- my,0,0,0,0],[0,0,0,mzz,mzx,mz,- my],[mzz,mzy,mz,0,0,0,- mx]])
        F = 1 / N * np.array([[loads["Fxb"][bk]],[loads["Fyb"][bk]],[loads["Fzb"][bk]],[loads["Mzb"][bk]],[0],[loads["rBlade"][bk] * loads["Fyb"][bk]],[loads["rBlade"][bk] * loads["Fxb"][bk]]])
        abc = np.linalg.solve(A,F)
        ## END
        forcemap["fx"][i] = abc[0] * z + abc[1] * y + abc[2]
        forcemap["fy"][i] = abc[3] * z + abc[4] * x + abc[5]
        ## EMA added:
        forcemap["fz"][i] = abc[6]
        ## END
    
    return forcemap
    
    
def check_sums(nodeData, loads, forcemap):
    forcesums = dict()
    forcesums["Z"] = loads["rBlade"]
    forcesums["Fx"][:,0] = loads["Fxb"]
    forcesums["Fy"][:,0] = loads["Fyb"]
    forcesums["M"][:,0] = loads["Mzb"]
    ## EMA original:
    # forcesums["RootMx"][:,0] = loads["rBlade"] .* loads["Fyb"];
    # forcesums["RootMy"][:,0] = loads["rBlade"] .* loads["Fxb"];
    ## changed to:
    forcesums["RootMx"][:,0] = np.multiply(- loads["rBlade"],loads["Fyb"]) + np.multiply(loads["prebend"],loads["Fzb"])
    forcesums["RootMy"][:,0] = np.multiply(loads["rBlade"],loads["Fxb"]) - np.multiply(loads["presweep"],loads["Fzb"])
    ## END
    
    for bk in range(len(loads["Fxb"])):
        i = np.where(forcemap["bin"] == bk)
        x = nodeData[i,1] - loads["presweep"][bk]
        y = nodeData[i,2] - loads["prebend"][bk]
        z = nodeData[i,3]
        forcesums["Fx"][bk,1] = sum(forcemap["fx"][i])
        forcesums["Fy"][bk,1] = sum(forcemap["fy"][i])
        forcesums["M"][bk,1] = sum(np.multiply(- y,forcemap["fx"][i]) + np.multiply(x,forcemap["fy"][i]))
        ## EMA original:
        #     forcesums["RootMx"](bk,2) = sum(z.*forcemap["fy"][i]);
        #     forcesums["RootMy"](bk,2) = sum(z.*forcemap["fx"][i]);
        ## changed to:
        x = nodeData[i,0]
        y = nodeData[i,1]
        forcesums["RootMx"][bk,0] = sum(np.multiply(- z,forcemap["fy"][i]) + np.multiply(y,forcemap["fz"][i]))
        forcesums["RootMy"][bk,0] = sum(np.multiply(z,forcemap["fx"][i]) - np.multiply(x,forcemap["fz"][i]))
        ## END
    
    return forcesums
    
    
def writeforcefile(filename, forcemap, forcesums, maptype): 
    fid = open(filename,'wt')
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    if ('forcesums' is not None):
        fid.write('!========== FORCE MAPPING SUMMARY ==========' % ())
        fid.write('\n!maptype = "%s"' % (maptype))
        s = print('\n!                Z =%s      TOTAL     ',('  %14.6e',forcesums["Z"]))
        fid.write('%s' % (s))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input          Fx =' % ())
        fid.write('  %14.6e' % (forcesums["Fx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,0])))
        fid.write('\n!Output    sum(fx) =' % ())
        fid.write('  %14.6e' % (forcesums["Fx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input          Fy =' % ())
        fid.write('  %14.6e' % (forcesums["Fy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,0])))
        fid.write('\n!Output    sum(fy) =' % ())
        fid.write('  %14.6e' % (forcesums["Fy"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input           M =' % ())
        fid.write('  %14.6e' % (forcesums["M"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,0])))
        fid.write('\n!sum(-y*fx + x*fy) =' % ())
        fid.write('  %14.6e' % (forcesums["M"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input        Z*Fy =' % ())
        fid.write('  %14.6e' % (forcesums["RootMx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,0])))
        fid.write('\n!Output  sum(z*fy) =' % ())
        fid.write('  %14.6e' % (forcesums["RootMx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input        Z*Fx =' % ())
        fid.write('  %14.6e' % (forcesums["RootMy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMy"][:,0])))
        fid.write('\n!Output  sum(z*fx) =' % ())
        fid.write('  %14.6e' % (forcesums["RootMy"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["RootMy"][:,1])))
        fid.write('\n\n' % ())
    
    fid.write('finish\n/prep7\n\n' % ())
    for nk in range(len(forcemap["n"])):
        if forcemap["fx"][nk]:
            fid.write('f,%d,fx,%g\n' % (forcemap["n"][nk],forcemap["fx"][nk]))
        if forcemap["fy"][nk]:
            fid.write('f,%d,fy,%g\n' % (forcemap["n"][nk],forcemap["fy"][nk]))
        ## EMA added:
        if forcemap["fz"][nk]:
            fid.write('f,%d,fz,%g\n' % (forcemap["n"][nk],forcemap["fz"][nk]))
        ## END
    
    fid.close()
    return

    
def beamForceToAnsysShellFollower(maptype, nodeData, loads, outfile): 
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
    # open file selector if 'outfile' not given
    if not ('outfile' is not None)  or len(outfile)==0:
        outfile = join(nlistpath,'forces.src')
        fn,pn = uiputfile(np.array(['*.*','All files (*.*)']),'Save Output As',outfile)
        if fn==0:
            # the user canceled file selection
            return
        outfile = join(pn,fn)
    
    # write the ansys commands which apply the calculated forces
    writeforcefile(outfile,forcemap,forcesums,maptype)
    # plot the results
    # if 1:
    #     plt.figure(1007)
    #     subplot(2,2,np.array([1,3]))
    #     ds = np.ceil(nodeData.shape[1-1] / 1000.0)
    #     k = np.arange(1,nodeData.shape[1-1]+ds,ds)
    #     quiver3(nodeData[k,1],nodeData[k,2],nodeData[k,3],forcemap["fx"][k],forcemap["fy"][k],np.zeros((forcemap["fx"][k].shape,forcemap["fx"][k].shape)))
    #     plt.title('Force vectors visual check (downsampled)')
    #     plt.axis('equal')
    #     subplot(2,2,2)
    #     __,k = __builtint__.sorted(nodeData[:,3])
    #     plt.plot(nodeData[k,3],forcemap["fx"][k],'.')
    #     plt.xlabel('z')
    #     plt.ylabel('fx')
    #     subplot(2,2,4)
    #     plt.plot(nodeData[k,3],forcemap["fy"][k],'.')
    #     plt.xlabel('z')
    #     plt.ylabel('fy')
    
    return
    
    
def read_forces(filename): 
    # Open the file
    fid = open(filename)
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    header = fid.readline()
    
    filecontents = np.transpose(fid.read())
    fid.close()
    # 'Z (m)	Fx (N)	Fy (N)	M (N-m)	Alpha	x_off	y_off'
    forces = textscan(filecontents,np.matlib.repmat('%f',1,7))
    return forces
    
    
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
            OBedge = loads["rBlade"](bk) + halfDZ
            if (mesh["z"](nk) <= OBedge) or (bk == np.asarray(loads["rBlade"]).size):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"](bk + 1) - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    forcemap["fy"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    for bk in np.arange(1,len(loads["Fxb"])+1).reshape(-1):
        i = np.where(bin == bk)
        N = len(i)
        x = mesh["x"](i) - loads["presweep"](bk)
        y = mesh["y"](i) - loads["prebend"](bk)
        mx = np.mean(x)
        my = np.mean(y)
        mxx = np.mean(x ** 2)
        myy = np.mean(y ** 2)
        #mxy = np.mean(x.*y);
        A = np.array([[my,1,0,0],[0,0,mx,1],[0,0,mxx,mx],[- myy,- my,0,0]])
        F = 1 / N * np.array([[loads["Fxb"](bk)],[loads["Fyb"](bk)],[loads["Mzb"](bk)],[0]])
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
    bin = np.zeros((mesh["n"].shape,mesh["n"].shape))
    
    for nk in np.arange(1,len(mesh["n"])+1).reshape(-1):
        halfDZ = loads["rBlade"][0]
        for bk in np.arange(1,np.asarray(loads["rBlade"]).size+1).reshape(-1):
            OBedge = loads["rBlade"](bk) + halfDZ
            if (mesh["z"](nk) <= OBedge) or (bk == np.asarray(loads["rBlade"]).size):
                bin[nk] = bk
                break
            halfDZ = loads["rBlade"](bk + 1) - OBedge
    
    forcemap = dict()
    forcemap["n"] = nodeData[:,0]
    forcemap["bin"] = bin
    forcemap["fx"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    forcemap["fy"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    ## EMA added:
    forcemap["fz"] = np.zeros((mesh["n"].shape,mesh["n"].shape))
    ## END
    for bk in np.arange(1,len(loads["Fxb"])+1).reshape(-1):
        i = np.where(bin == bk)
        N = len(i)
        x = mesh["x"](i) - loads["presweep"](bk)
        y = mesh["y"](i) - loads["prebend"](bk)
        z = mesh["z"](i)
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
        #     F = 1/N * [aero.Fx(bk); aero.Fy(bk); aero.M(bk); 0; ...
        #         aero.Z(bk)*aero.Fy(bk); aero.Z(bk)*aero.Fx(bk)];
        #     abc = A\F;
        ## Changed to:
        ## Add/modify equations to include forces in the Z-direction.
        A = np.array([[mz,my,1,0,0,0,0],[0,0,0,mz,mx,1,0],[0,0,0,0,0,0,1],[0,0,0,mzx,mxx,mx,0],[- mzy,- myy,- my,0,0,0,0],[0,0,0,mzz,mzx,mz,- my],[mzz,mzy,mz,0,0,0,- mx]])
        F = 1 / N * np.array([[loads["Fxb"](bk)],[loads["Fyb"](bk)],[loads["Fzb"](bk)],[loads["Mzb"](bk)],[0],[loads["rBlade"](bk) * loads["Fyb"](bk)],[loads["rBlade"](bk) * loads["Fxb"](bk)]])
        abc = np.linalg.solve(A,F)
        ## END
        forcemap["fx"][i] = abc[0] * z + abc[1] * y + abc[2]
        forcemap["fy"][i] = abc[3] * z + abc[4] * x + abc[5]
        ## EMA added:
        forcemap["fz"][i] = abc[6]
        ## END
    
    return forcemap
    
    
def check_sums(nodeData = None,loads = None,forcemap = None):
    forcesums = dict()
    forcesums["Z"] = loads["rBlade"]
    forcesums["Fx"][:,1] = loads["Fxb"]
    forcesums["Fy"][:,1] = loads["Fyb"]
    forcesums["M"][:,1] = loads["Mzb"]
    ## EMA original:
# forcesums["RootMx"][:,0] = loads["rBlade"] .* loads["Fyb"];
# forcesums["RootMy"][:,0] = loads["rBlade"] .* loads["Fxb"];
## changed to:
    forcesums["RootMx"][:,1] = np.multiply(- loads["rBlade"],loads["Fyb"]) + np.multiply(loads["prebend"],loads["Fzb"])
    forcesums["RootMy"][:,1] = np.multiply(loads["rBlade"],loads["Fxb"]) - np.multiply(loads["presweep"],loads["Fzb"])
    ## END
    
    for bk in np.arange(1,len(loads["Fxb"])+1).reshape(-1):
        i = np.where(forcemap["bin"] == bk)
        x = nodeData[i,1] - loads["presweep"](bk)
        y = nodeData[i,2] - loads["prebend"](bk)
        z = nodeData[i,3]
        forcesums["Fx"][bk,2] = sum(forcemap["fx"](i))
        forcesums["Fy"][bk,2] = sum(forcemap["fy"](i))
        forcesums["M"][bk,2] = sum(np.multiply(- y,forcemap["fx"](i)) + np.multiply(x,forcemap["fy"](i)))
        ## EMA original:
#     forcesums["RootMx"](bk,2) = sum(z.*forcemap["fy"](i));
#     forcesums["RootMy"](bk,2) = sum(z.*forcemap["fx"](i));
## changed to:
        x = nodeData[i,1]
        y = nodeData[i,2]
        forcesums["RootMx"][bk,2] = sum(np.multiply(- z,forcemap["fy"](i)) + np.multiply(y,forcemap["fz"](i)))
        forcesums["RootMy"][bk,2] = sum(np.multiply(z,forcemap["fx"](i)) - np.multiply(x,forcemap["fz"](i)))
        ## END
    
    return forcesums
    
    
def writeforcefile(filename, forcemap, forcesums, maptype): 
    fid = open(filename,'wt')
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    if ('forcesums' is not None):
        fid.write('!========== FORCE MAPPING SUMMARY ==========' % ())
        fid.write('\n!maptype = "%s"' % (maptype))
        s = print('\n!                Z =%s      TOTAL     ',('  %14.6e',forcesums["Z"]))
        fid.write('%s' % (s))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input          Fx =' % ())
        fid.write('  %14.6e' % (forcesums["Fx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,0])))
        fid.write('\n!Output    sum(fx) =' % ())
        fid.write('  %14.6e' % (forcesums["Fx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input          Fy =' % ())
        fid.write('  %14.6e' % (forcesums["Fy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,0])))
        fid.write('\n!Output    sum(fy) =' % ())
        fid.write('  %14.6e' % (forcesums["Fy"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input           M =' % ())
        fid.write('  %14.6e' % (forcesums["M"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,0])))
        fid.write('\n!sum(-y*fx + x*fy) =' % ())
        fid.write('  %14.6e' % (forcesums["M"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input        Z*Fy =' % ())
        fid.write('  %14.6e' % (forcesums["RootMx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,0])))
        fid.write('\n!Output  sum(z*fy) =' % ())
        fid.write('  %14.6e' % (forcesums["RootMx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,1])))
        fid.write('\n!%s' % (np.matlib.repmat('-',1,len(s))))
        fid.write('\n!Input        Z*Fx =' % ())
        fid.write('  %14.6e' % (forcesums["RootMy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMy"][:,0])))
        fid.write('\n!Output  sum(z*fx) =' % ())
        fid.write('  %14.6e' % (forcesums["RootMy"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["RootMy"][:,1])))
        fid.write('\n\n' % ())
    
    fid.write('finish\n/prep7\n\n' % ())
    fid.write('csys,0\n' % ())
    fid.write('allsel\n' % ())
    fid.write('et,33,FOLLW201\n' % ())
    fid.write('KEYOPT,33,2,0\n' % ())
    fid.write('TYPE, 33\n\n' % ())
    
    for nk in range(len(forcemap["n"])):
        F = np.array([forcemap["fx"](nk),forcemap["fy"](nk),forcemap["fz"](nk)])
        vlen = np.linalg.norm(F,2)
        FX = forcemap["fx"](nk) / vlen
        FY = forcemap["fy"](nk) / vlen
        FZ = forcemap["fz"](nk) / vlen
        fid.write('R,%d,%g,%g,%g,0,0,0 \n' % (forcemap["n"](nk),FX,FY,FZ))
        fid.write('REAL, %d\n' % (forcemap["n"](nk)))
        fid.write('E,%d \n' % (forcemap["n"](nk)))
        #Select follower element without knowing its element number
        fid.write('NSEL, S,NODE,,%d\n' % (forcemap["n"](nk)))
        fid.write('ESLN, S, 1\n' % ())
        #Apply Load
        if vlen:
            fid.write('sfe,all,1,pres,,%g\n' % (vlen))
        fid.write('allsel\n\n' % ())
        #     if forcemap["fx"](nk)  # if not zero
#         fprintf(fid,'f,#d,fx,#g\n',forcemap["n"](nk),forcemap["fx"](nk));
#     end
#     if forcemap["fy"](nk)  # if not zero
#         fprintf(fid,'f,#d,fy,#g\n',forcemap["n"](nk),forcemap["fy"](nk));
#     end
    
    fid.close()
    return

    
def getMatrialLayerInfoWithOutGUI(blade): 
    #Temparary workaround to extract data needed for:
    
    app = getApp('numad.nmd','ansys',blade.mesh)
    #From write_shell7.m
    TotalStations = app.station.size
    TotalShearwebs = app.shearweb.size
    skin_areas = []
    for kStation in range(skin_areas.size):
        stationIB = app.station(kStation)
        stationOB = app.station(kStation + 1)
        for kdp in range(stationIB.size - 1):
            if np.array(['single','double']) == stationIB.dptype[kdp]:
                # single and double are equivalent on the area inboard edge
                # start and end are current and next DP
                skin_areas[kStation].startIB[-1] = kdp
                skin_areas[kStation].endIB[-1] = kdp + 1
                skin_areas[kStation].Material[-1] = stationIB.sm[kdp]
            else:
                if np.array(['flare','hourglass']) == stationIB.dptype[kdp]:
                    # flare and hourglass are equivalent on the area inboard edge
                    # start and end of first area is current DP
                    # start and end of next area is current and next DP
                    skin_areas[kStation].startIB[end() + [np.arange[1,2+1]]] = np.array([kdp,kdp])
                    skin_areas[kStation].endIB[end() + [np.arange[1,2+1]]] = np.array([kdp,kdp + 1])
                    skin_areas[kStation].Material[end() + 1] = stationIB.dpmaterial[kdp]
                    skin_areas[kStation].Material[end() + 1] = stationIB.sm[kdp]
        for kdp in np.arange(1,np.asarray(stationOB.dp).size - 1+1).reshape(-1):
            if np.array(['single','flare']) == stationOB.dptype[kdp]:
                # single and flare are equivalent on the area outboard edge
                # start and end are current and next DP
                skin_areas[kStation].startOB[end() + 1] = kdp
                skin_areas[kStation].endOB[end() + 1] = kdp + 1
            else:
                if np.array(['double','hourglass']) == stationOB.dptype[kdp]:
                    # double and hourglass are equivalent on the area outboard edge
                    # start and end of first area is current DP
                    # start and end of next area is current and next DP
                    skin_areas[kStation].startOB[end() + [np.arange[1,2+1]]] = np.array([kdp,kdp])
                    skin_areas[kStation].endOB[end() + [np.arange[1,2+1]]] = np.array([kdp,kdp + 1])
    
    #tcl: Determine which composite materials are used in the model
    compsInModel = np.array([])
    #tcl: search shear web materials
    for k in np.arange(1,TotalShearwebs+1).reshape(-1):
        compsInModel = np.array([[compsInModel],[np.array([app.shearweb(k).Material])]])
    
    #tcl: search skin materials
    for k in np.arange(1,np.asarray(skin_areas).size+1).reshape(-1):
        compsInModel = np.array([[compsInModel],[transpose(skin_areas(k).Material)]])
    
    compsInModel = np.unique(compsInModel)
    # load the material database and create searchable list
    if not len(app.settings.job_name)==0 :
        # job name exists, use the local material database
        app.matdb_path = fullfile(app.settings.job_path,'MatDBsi.txt')
    else:
        # if no job name, use the master material database
        # jcb: we shouldn't arrive here because a file save is required first
        app.matdb_path = fullfile(app.numadpath,'MatDBsi.txt')
    
    app.matdb = readMatDB(app.matdb_path)
    for k in np.arange(1,np.asarray(app.matdb).size+1).reshape(-1):
        app.matlist[k] = app.matdb(k).name
        mattype[k] = app.matdb(k).type
    
    app.isotropic = str('isotropic') == str(mattype)
    app.orthotropic = str('orthotropic') == str(mattype)
    app.composite = str('composite') == str(mattype)
    # Determine which isotropic and orthotropic materials are used in the model
    isoorthoInModel = np.array([])
    for kcomp in np.arange(1,np.asarray(compsInModel).size+1).reshape(-1):
        n = str(compsInModel(kcomp)) == str(app.matlist)
        if not np.any(n) :
            errordlg(sprintf('Material "%s" not found in database.',compsInModel[kcomp]),'Error')
            raise Exception('Material "%s" not found in database.',compsInModel[kcomp])
        mat = app.matdb(n)
        layerNames = []
        for klay in range(mat.layer.size):
            layerNames.append(np.array([mat.layer(klay).layerName]))
        isoorthoInModel = np.unique(np.array([[isoorthoInModel],[layerNames]]))
    
    return isoorthoInModel,compsInModel,skin_areas,app


def txt2mat(filename):
    """
    Reads text file containing lines of space separated numbers and
    converts to a matrix of floats.

    Parameters
    ----------
    filename : str

    Returns
    -------
    mat : numpy array
    """
    with open(filename) as fid:
        lines = fid.readlines()
    rows = []
    for line in lines:
        line = line.strip()
        nums = line.split(' ')
        nums = list(filter(''.__ne__,nums))
        nums = list(map(float,nums))
        row = np.array(nums)
        rows.append(row)
    mat = np.array(rows)
    return mat


def postprocessANSYSfatigue(blade, meshData, wt, rccdata, IEC, loadsTable, config): 
    if np.any('all' in config.analysisFlags.fatigue.lower()): # NOTE probably need to workshop this -kb
        nSegments = 1
    else:
        nSegments = np.asarray(config.analysisFlags.fatigue).size
    
    # Order of the segment names in segmentNamesReference
    # is very important. config.analysisFlags.fatigue can
    # be any order
    segmentNamesReference = ['HP_TE_FLAT','HP_TE_ReINF','HP_TE_PANEL','HP_SPAR','HP_LE_PANEL','HP_LE','LP_LE','LP_LE_PANEL','LP_SPAR','LP_TE_PANEL','LP_TE_REINF','LP_TE_FLAT']
    nsegmentNamesReference = len(segmentNamesReference)
    markovSize = 16
    designVar = np.array([])
    
    Yr = IEC.designLife
    #fst=readFastMain(['IEC_' IEC.fstfn '.fst']);
    #simtime=IEC.numSeeds*(fst.SimCtrl.TMax-IEC.delay); # simulated and rainflow counted time, seconds
    simtime = IEC.numSeeds * (IEC.SimTime - IEC.delay)
    nSpace = 90 / loadsTable[2].theta
    
    nDirections = len(loadsTable)
    
    for kTheta in np.arange(1,nDirections / 2+1).reshape(-1):
        #since blade movements along single direction constitues
        #two directions (e.g positive flap deflections and negative
        #ones are two directions; both wich make
        #the flap cycles)
        loadsTableTheta = loadsTable[kTheta]
        theta = loadsTableTheta.theta
        loadsTableThetaPlus90 = loadsTable[kTheta + nSpace]
        nGage = loadsTableTheta.input.rGagesize
        gageNumber = np.transpose((np.arange(nGage)))
        criticalElement,fatigueDamage,criticalLayerNo,criticalMatNo = deal(np.zeros((nGage,1)))
        criticalMat = cell(nGage,1)
        rGage = loadsTableTheta.input.rGage
        MrTheta = loadsTableTheta.input.Mrb
        MrThetaPlus90 = loadsTableThetaPlus90.input.Mrb
        plotFatigue = []
        fileNameTheta = 'plateStrains-all-'+str(kTheta)+'.txt'
        fileNameThetaPlus90 = 'plateStrains-all-'+str(kTheta + nSpace)+'.txt'
        print(fileNameTheta)
        print(fileNameThetaPlus90)
        #Used for reading element stresses
        pat = 'ELEM\s*ZCENT\s*EPS11\s*EPS22\s*EPS12\s*KAPA11\s*KAPA22\s*KAPA12\s*GAMMA13\s*GAMMA23'
        NCOLS = 10
        plateStrainsTheta = readANSYSElementTable(fileNameTheta,pat,NCOLS)
        plateStrainsThetaPlus90 = readANSYSElementTable(fileNameThetaPlus90,pat,NCOLS)
        for i in range(nSegments):
            iSegment = np.where(segmentNamesReference == config.analysisFlags.fatigue[i])
            if np.any('all' in config.analysisFlags.fatigue.lower()):
                title = 'All segments'
            else:
                if not 'webs' == config.analysisFlags.fatigue[i] :
                    title = config.analysisFlags.fatigue[i]
                    __,nSpanRegions = meshData.outerShellElSets.shape
                    elementList = []
                    for iSpan in range(nSpanRegions):
                        elementList = [elementList,meshData.outerShellElSets[iSegment,iSpan].elementList]
                else:
                    title = 'Webs'
                    __,nWebs = meshData.shearWebElSets.shape
                    elementList = []
                    for iWeb in range(nWebs):
                        __,nSpanRegions = meshData.shearWebElSets[iWeb].shape
                        for iSpan in range(nSpanRegions):
                            elementList = np.array([elementList,meshData.shearWebElSets[iWeb][iSpan].elementList])
                plateStrainsThetaSet = plateStrainsTheta[elementList,:]
                plateStrainsThetaPlus90Set = plateStrainsThetaPlus90[elementList,:]
            for chSpan in range(nGage):
                direction = str(theta)
                Ltheta = getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction)
                if theta + 90 < 180:
                    direction = str(theta + 90)
                else:
                    direction = str(theta - 90)
                LthetaPlus90 = getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction)
                Mtheta = interpolator_wrap(rGage,MrTheta,rGage[chSpan])
                MthetaPlus90 = interpolator_wrap(rGage,MrThetaPlus90,rGage[chSpan])
                zwidth = 0.75
                # at a blade gage location.
                z1 = rGage[chSpan] - zwidth / 2
                z2 = rGage[chSpan] + zwidth / 2
                binnedElements = np.intersect(np.where(plateStrainsThetaSet[:,1] < z2),np.where(plateStrainsThetaSet(:,2) > z1))
                fdData,plotFatigueChSpan = calcFatigue(blade,meshData,IEC,Ltheta,LthetaPlus90,Mtheta,MthetaPlus90,binnedElements,plateStrainsThetaSet,plateStrainsThetaPlus90Set,iSegment)
                plotFatigue = np.array([[plotFatigue],[plotFatigueChSpan]])
                criticalElement[chSpan] = fdData[0]
                fatigueDamage[chSpan] = fdData[1]
                criticalLayerNo[chSpan] = fdData[4]
                criticalMatNo[chSpan] = fdData[7]
                criticalMat[chSpan] = blade.materials(fdData(8)).name
            #         plotFatigueFileName=['plotFatigue-' str(kTheta)];
#         writePlotFatigue(plotFatigueFileName,plotFatigue)
            print('\n\n\n ************************ Segment No-%i: %s ************************\n' % (i,title))
            table(gageNumber,criticalElement,fatigueDamage,criticalLayerNo,criticalMatNo,criticalMat)
            #         designVar{end+1}=max(fatigueDamage);
            designVar[kTheta].fatigueDamage[i,:] = fatigueDamage
            designVar[kTheta].criticalElement[i,:] = criticalElement
            designVar[kTheta].criticalLayerNo[i,:] = criticalLayerNo
            designVar[kTheta].criticalMatNo[i,:] = criticalMatNo
    
    #delete stresses-*-*.txt;
    return designVar


def writePlotFatigue(fname, plotFatigue): 
    #Write fatigue damage for each element. ANSYS requires elements to be
    #sorted
    n = len(plotFatigue[:,0])
    fid = open(np.array([fname,'.txt']),'w+')
    fid.write('Element fatigueDamage\n' % ())
    plotFatigue = sortrows(plotFatigue,1)
    for i in range(len(plotFatigue[:,0])):
        fid.write('%8i  %6.5E\n' % (plotFatigue[i,0],plotFatigue[i,1]))
    
    fid.close()
    #Write plot commands
    fid = open(fname+'.mac'),'w+')
    fid.write('/post1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('plnsol,u,sum\n' % ())
    fid.write('etab,test,u,X\n' % ())
    fid.write('*get,max_e,elem,0,count\n' % ())
    fid.write('*dim,d_res,array,max_e,3\n' % ())
    fid.write('!Column 1 = Element Number\n' % ())
    fid.write("!Column 2 = Where I'm putting result data\n" % ())
    fid.write('*vget,d_res(1,1),elem,,elist\n' % ())
    fid.write('*vfill,d_res(1,2),ramp,0,0\n' % ())
    fid.write('*dim,d_results,array,%i,2   !Need to specify the same size array as the data being read in\n' % (n))
    fid.write('*vread,d_results(1,1),%s,txt,,jik,2,%i,,1\n' % (fname,n))
    fid.write('(F8.0,E13.5)\n' % ())
    fid.write('*get,d_temp,parm,d_results,dim,x\n' % ())
    fid.write('j=1\n' % ())
    fid.write('i=1\n' % ())
    fid.write('d_run=1\n' % ())
    fid.write('*dowhile,d_run\n' % ())
    fid.write('*if,d_res(i,1),EQ,d_results(j,1),THEN\n' % ())
    fid.write('d_res(i,2)=d_results(j,2)\n' % ())
    fid.write('j=j+1\n' % ())
    fid.write('*endif\n' % ())
    fid.write('*if,j,GT,d_temp,THEN\n' % ())
    fid.write('d_run=0\n' % ())
    fid.write('*endif\n' % ())
    fid.write('i=i+1\n' % ())
    fid.write('*enddo\n' % ())
    #     fprintf(fid,'j=1\n');
    #     fprintf(fid,'*do,i,1,max_e\n');
    #     fprintf(fid,'*if,j,LT,#i,THEN\n',n+1);
    #     fprintf(fid,'*if,d_res(i,1),EQ,d_results(j,1),THEN\n');
    #     fprintf(fid,'d_res(i,2)=d_results(j,2)\n');
    #     fprintf(fid,'j=j+1\n');
    #     fprintf(fid,'*endif\n');
    #     fprintf(fid,'*endif\n');
    #     fprintf(fid,'*enddo\n');
    
    fid.write('allsel,all\n' % ())
    fid.write('*vput,d_res(1,2),elem,1,etab,test\n' % ())
    fid.write('Pretab\n' % ())
    fid.write('pletab,test\n' % ())
    
    return

    
def getWindSpeedDistribution(avgws): 
    scipy.io.loadmat('rccdata.mat','rccdata')
    # determine the wind speeds that are saved in rccdata
    for w in range(rccdata.shape[1]):
        ws[w] = rccdata[1,w].windspeed
    
    # check to make sure wind speeds are spaced evenly
    if std(np.diff(ws)) != 0:
        print(ws)
        raise Exception('Your windspeeds are not spaced evenly')
    
    # define bin edges based on wind speeds in rccdata
    binwidth = ws[1] - ws[0]
    binedges = np.array([ws[0] - binwidth / 2,ws + binwidth / 2])
    # define wind bins
    windbins = np.zeros((len(binedges),2))
    for jj in range(len(binedges)-1):
        windbins[jj,:] = np.array([binedges[jj],binedges[jj + 1]])
    
    # Calculate weights for each bin according to Rayleigh distribution
    sig = avgws / np.sqrt(np.pi / 2)
    # find PDF and CDF of Rayleigh distribution
    #     pdf=windbins.*exp(-windbins.^2/(2*sig^2))/sig^2;
    cdf = 1 - np.exp(- windbins ** 2 / (2 * sig ** 2))
    # calculate weights
    wt = np.diff(cdf,1,2)
    # show sum of weights (should be close to one, and not greater than one)
    print(' ')
    print('Sum of the Rayleigh weights is ', sum(wt))
    print(' ')
    # check to make sure that the number of weights is equal to the number of
    # simulated wind speeds in rccdata
    if len(wt) != rccdata.shape[1]:
        raise Exception('The number of Rayleigh weights is not equal to the number of wind speeds contained in the rccdata.mat file')
    
    return wt,rccdata


    
def getLoadFactorsForElementsWithSameSection(LF, ansysSecNumber, avgFaceStress,app, mat, coreMatName): 
    #This is a recursive function for LF. It appends to the list of LF for each
    #element that has a positive LF.  EC
    m,__ = avgFaceStress.shape
    for i in np.arange(1,m+1).reshape(-1):
        elno = avgFaceStress(i,1)
        S11a = avgFaceStress(i,2)
        S22a = avgFaceStress(i,3)
        S12a = avgFaceStress(i,7)
        #Ignoring other stresses for the time being.
        #if elno==305
        #disp('press pause')
        #pause(10)
        lf,phicr = checkWrinkle(np.array([[S11a],[S22a],[S12a]]),mat,app,coreMatName)
        if lf >= 0:
            LF = np.array([[LF],[ansysSecNumber,elno,lf,phicr]])
        #end
    
    
    def checkWrinkle(S_alphaBeta, mat, app, coreMatName): 
        # For a single finite element, given the average in-plane stresses
        # in a face-sheet of that element, compute the load factor for that
        # element. EC
        
        # lf    - scalar load factor for the element
        # phicr - an angle, degrees. The direction of wrinkling
        # S11a,S22a,S12a - respective average face sheet stress
        # mat - material object
        # app - blade data
        
        #Locate the face sheet
        cellMat = np.array([])
        for i in range(len(mat.layer)):
            cellMat = np.array([[cellMat],[np.array([mat.layer[i].layerName])]])
        
        kbalsa = np.find(str(coreMatName) == str(cellMat))
        iLayer = np.arange(1,(kbalsa - 1)+1)
        
        #ilayer=(kbalsa+1):numel(cellMat)); #Number of distinct materials in the bottom face
        matCore = app.matdb(np.find(str(coreMatName) == str(app.matlist)))
        if str(matCore.type) == str('orthotropic'):
            Ec = matCore.ez
        else:
            if str(matCore.type) == str('isotropic'):
                Ec = matCore.ex
                Gc = Ec / (2 * (1 + matCore.nuxy))
            else:
                print('Material  "%s"  not found in database.',matCore.type)
                raise Exception('Material type "%s" not found in database.',matCore.type)
        
        #ilayer=(kbalsa+1):numel(cellMat)); #Number of distinct materials in the bottom face
        dangle = 2
        
        N = 180 / dangle + 1
        
        angle = 0
        invLF = np.zeros((N,1))
        #Apt=zeros(3,3);
        #Bpt=zeros(3,3);
        Dpt = np.zeros((3,3))
        #Find total height of facesheet
        h = 0
        for klay in range(iLayer.size):
            h = h + mat.layer(iLayer[klay]).thicknessA * mat.layer(iLayer[klay]).quantity
        
        for kang in range(N):
            if str(matCore.type) == str('orthotropic'):
                Gc = 1 / (np.sin(np.pi/180*angle) ** 2 * (1 / matCore.gyz) + np.cos(np.pi/180*angle) ** 2 * (1 / matCore.gxz))
            #Asssuming all ply angles are zero
            # R_sig=[cosd(angle)^2, sind(angle)^2, -2*sind(angle)*cosd(angle);  #In-plate clockwise rotation of: angle
            #        sind(angle)^2, cosd(angle)^2, 2*sind(angle)*cosd(angle)
            #        sind(angle)*cosd(angle), -sind(angle)*cosd(angle), cosd(angle)^2-sind(angle)^2];
            z1 = - h / 2
            for klay in range(iLayer.size):
                z2 = z1 + mat.layer(iLayer(klay)).thicknessA * mat.layer(iLayer(klay)).quantity
                matklay = app.matdb(np.find(str(mat.layer(iLayer(klay)).layerName) == str(app.matlist)))
                # Bulid Plane Stress reduced compliance matrix for each
                # layer
                #fprintf('z1 = #f z2 = #f mat = #f   #s\n',z1,z2,matListnumber,mat.layer(klay).layerName)
                #Entries common to either isotropic or orthotropic entries
                Se = np.zeros((3,3))
                Se[1,1] = 1 / matklay.ex
                Se[1,3] = 0
                Se[2,3] = 0
                Se[3,1] = 0
                Se[3,2] = 0
                if str(matklay.type) == str('orthotropic'):
                    Se[1,2] = - matklay.prxy / matklay.ex
                    Se[2,1] = - matklay.prxy / matklay.ex
                    Se[2,2] = 1 / matklay.ey
                    Se[3,3] = 1 / matklay.gxy
                else:
                    if str(matklay.type) == str('isotropic'):
                        Se[1,2] = - matklay.nuxy / matklay.ex
                        Se[2,1] = - matklay.nuxy / matklay.ex
                        Se[2,2] = 1 / matklay.ex
                        Se[3,3] = 2 * (1 + matklay.nuxy) / matklay.ex
                    else:
                        print('Material  "%s"  not found in database.',matklay.type)
                        raise Exception('Material type "%s" not found in database.',matklay.type)
                #Apt=Apt+R_sig*inv(Se)*R_sig'*(z2-z1);
                #Bpt=Bpt+1/2*R_sig*inv(Se)*R_sig'*(z2^2-z1^2);
                Dpt = Dpt + 1 / 3 * R_sig * inv(Se) * np.transpose(R_sig) * (z2 ** 3 - z1 ** 3)
                z1 = z2
            Pcr = - 3 / 2 * (2 * Dpt(1,1) * Ec * Gc) ** (1 / 3)
            Pphi = (R_sig(1,:) * S_alphaBeta) * h
            invLF[kang] = Pphi / Pcr
            angle = angle + dangle
        
        invlf,phicr_i = np.amax(invLF)
        lf = 1 / invlf
        phicr = (phicr_i - 1) * dangle
        #         if lf>1e6
        #             figure(2)
        #             plot(angle, Pcr,'k')
        #             xlabel('Angle, \phi [deg]')
        #             ylabel('Load [N/m]')
        #             hold on;
        
        #             plot(angle, Pphi,'r')
        #             plot(angle, lf*Pphi,'b')
        #             legend('P_c_r','P_\phi',strcat(num2str(lf),'P_\phi'))
        #             hold off;
        #             figure(3)
        #             plot(angle,invLF,'r')
        #             hold on;
        #             plot(phicr,invlf,'*')
        #             xlabel('Angle, \phi [deg]')
        #             ylabel('1/\lambda [ ]')
        #             hold off;
        #             fprintf('\n #8.2g  #i | #8.2g #8.2g | #8.2g #8.2g #8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
        #             pause(10)
        #         else
        #             #fprintf('\n #8.2g  #i   #8.2g #8.2g   #8.2g #8.2g #8.2g\n',lf,ansysSecNumber,Pcr(phicr_i),Pphi(phicr_i),S11a,S22a,S12a)
        #         end
        return lf,phicr
    
    return LF