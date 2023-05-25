import numpy as np
from pynumad.utils.interpolation import interpolator_wrap
import os
def readVABShomogenization(fileName):


    #Stiffness
    beam_stiff=np.zeros((6,6))
    with open(fileName) as f:
        lines=f.readlines()
        
    for lineNumber,line in enumerate(lines):
        if 'Timoshenko Stiffness Matrix' in line:
            break
    lineStart=lineNumber+3
    lineEnd=lineStart+6
    ct=0
    for iLine in range(lineStart,lineEnd):
        dataList=[float(i) for i in lines[iLine].split() if i.strip()]
        beam_stiff[ct,:]=dataList
        ct+=1

    
    #mass
    beam_inertia=np.zeros((6,6))
    for lineNumber,line in enumerate(lines):
        if 'Mass Matrix' in line:
            break
    lineStart=lineNumber+3
    lineEnd=lineStart+6
    ct=0
    for iLine in range(lineStart,lineEnd):
        dataList=[float(i) for i in lines[iLine].split() if i.strip()]
        beam_inertia[ct,:]=dataList
        ct+=1
  
    return beam_stiff, beam_inertia

    
def transformMatrixToBeamDyn(beam_stiff,beam_inertia):
    beamDynData={}

    B = np.array([[0, 0, 1], [0, -1, 0], [1, 0, 0]])  # NEW transformation matrix
    T = np.dot(np.identity(3), np.linalg.inv(B))
    
    nStations, _,_=np.shape(beam_stiff)

    for iStation in range(nStations):
        beam_stiff[iStation,:,:]=trsf_sixbysix(beam_stiff[iStation,:,:], T)
        beam_inertia[iStation,:,:]=trsf_sixbysix(beam_inertia[iStation,:,:], T)
   
    return(beam_stiff,beam_inertia)

def trsf_sixbysix(M, T):
    """
    Transform six-by-six compliance/stiffness matrix. 
    change of reference frame in engineering (or Voigt) notation.
    
    Parameters
    ----------
    M : np.ndarray
        6x6 Siffness or Mass Matrix
    T : np.ndarray
        Transformation Matrix
        
    Returns
    ----------
    res : np.ndarray
        Transformed 6x6 matrix
    """

    TS_1 = np.dot(np.dot(T.T, M[0:3, 0:3]), T)
    TS_2 = np.dot(np.dot(T.T, M[3:6, 0:3]), T)
    TS_3 = np.dot(np.dot(T.T, M[0:3, 3:6]), T)
    TS_4 = np.dot(np.dot(T.T, M[3:6, 3:6]), T)

    tmp_1 = np.vstack((TS_1, TS_2))
    tmp_2 = np.vstack((TS_3, TS_4))
    res = np.hstack((tmp_1, tmp_2))
    return res

# --- Write BeamDyn file with blade reference line locations ---#
def write_beamdyn_axis(directory, wt_name, blade,radial_stations):

    n_pts = 50
    grid = np.linspace(0, 1, n_pts)

    kp_xr=interpolator_wrap(radial_stations,blade.prebend,grid,'pchip', axis=1)
    kp_yr=interpolator_wrap(radial_stations,blade.sweep,grid,'pchip', axis=1)
    kp_zr=interpolator_wrap(radial_stations,blade.ispan,grid,'pchip', axis=1)
    twist_interp=interpolator_wrap(radial_stations,blade.idegreestwist,grid,'pchip', axis=1)


    data = np.vstack((kp_xr, kp_yr, kp_zr, twist_interp)).T

    if not os.path.exists(directory):
        os.makedirs(directory)

    axisFileName=wt_name + '_BeamDyn.dat'
    
    file = open(directory +'/'+ axisFileName, 'w')
    file.write('--------- BEAMDYN with OpenFAST INPUT FILE -------------------------------------------\n')
    file.write('%s blade\n' % (wt_name))
    file.write('---------------------- SIMULATION CONTROL --------------------------------------\n')
    file.write('True          Echo            - Echo input data to "<RootName>.ech" (flag)\n')
    file.write('True          QuasiStaticInit - Use quasistatic pre-conditioning with centripetal accelerations in initialization (flag) [dynamic solve only]\n')
    file.write(' 0            rhoinf          - Numerical damping parameter for generalized-alpha integrator\n')
    file.write(' 2            quadrature      - Quadrature method: 1=Gaussian; 2=Trapezoidal (switch)\n')
    file.write(' 2            refine          - Refinement factor for trapezoidal quadrature (-). DEFAULT = 1 [used only when quadrature=2]\n')
    file.write('"DEFAULT"     n_fact          - Factorization frequency (-). DEFAULT = 5\n')
    file.write('"DEFAULT"     DTBeam          - Time step size (s).\n')
    file.write(' 50           load_retries    - Number of factored load retries before quitting the aimulation\n')
    file.write('"DEFAULT"     NRMax           - Max number of iterations in Newton-Ralphson algorithm (-). DEFAULT = 10\n')
    file.write('"DEFAULT"     stop_tol        - Tolerance for stopping criterion (-)\n')
    file.write('"DEFAULT"     tngt_stf_fd     - Flag to use finite differenced tangent stiffness matrix (-)\n')
    file.write('"DEFAULT"     tngt_stf_comp   - Flag to compare analytical finite differenced tangent stiffness matrix  (-)\n')
    file.write('"DEFAULT"     tngt_stf_pert   - perturbation size for finite differencing (-)\n')
    file.write('"DEFAULT"     tngt_stf_difftol- Maximum allowable relative difference between analytical and fd tangent stiffness (-)\n')
    file.write('True          RotStates       - Orient states in the rotating frame during linearization? (flag) [used only when linearizing]\n')
    file.write('---------------------- GEOMETRY PARAMETER --------------------------------------\n')
    file.write('          1   member_total    - Total number of members (-)\n')
    file.write('         %u   kp_total        - Total number of key points (-) [must be at least 3]\n' % (n_pts))
    file.write('     1     %u                 - Member number; Number of key points in this member\n' % (n_pts))
    file.write('\t\t kp_xr \t\t\t kp_yr \t\t\t kp_zr \t\t initial_twist\n')
    file.write('\t\t  (m)  \t\t\t  (m)  \t\t\t  (m)  \t\t   (deg)\n')


    for i in range(n_pts):
        file.write('\t %.5e \t %.5e \t %.5e \t %.5e \n' % (data[i, 0], data[i, 1], data[i, 2], data[i, 3]))

    file.write('---------------------- MESH PARAMETER ------------------------------------------\n')
    file.write('          10   order_elem     - Order of interpolation (basis) function (-)\n')
    file.write('---------------------- MATERIAL PARAMETER --------------------------------------\n')
    file.write('"%s"    BldFile - Name of file containing properties for blade (quoted string)\n' % (wt_name + '_BeamDyn_Blade.dat'))
    file.write('---------------------- PITCH ACTUATOR PARAMETERS -------------------------------\n')
    file.write('False         UsePitchAct - Whether a pitch actuator should be used (flag)\n')
    file.write('        200   PitchJ      - Pitch actuator inertia (kg-m^2) [used only when UsePitchAct is true]\n')
    file.write('      2E+07   PitchK      - Pitch actuator stiffness (kg-m^2/s^2) [used only when UsePitchAct is true]\n')
    file.write('     500000   PitchC      - Pitch actuator damping (kg-m^2/s) [used only when UsePitchAct is true]\n')
    file.write('---------------------- OUTPUTS -------------------------------------------------\n')
    file.write('False          SumPrint       - Print summary data to "<RootName>.sum" (flag)\n')
    file.write('"ES10.3E2"    OutFmt         - Format used for text tabular output, excluding the time channel.\n')
    file.write('          1   NNodeOuts      - Number of nodes to output to file [0 - 9] (-)\n')
    file.write('3, 6, 9, 12, 15, 18, 21, 24, 27    OutNd          - Nodes whose values will be output  (-)\n')
    file.write('          OutList            - The next line(s) contains a list of output parameters. See OutListParameters.xlsx for a listing of available output channels, (-)\n')

    coordinate={}
    coordinate['F']='l'
    coordinate['M']='l'
    coordinate['RD']='r'
    coordinate['TD']='r'
    
    channelList=['F','M','RD','TD']
    for iNode in range(9):
        for load in channelList:
            for dir in ['x','y','z']:
                file.write(f'"N{iNode+1}{load}{dir}{coordinate[load]}"\n')

    
    #Root
    coordinate={}
    coordinate['F']='r'
    coordinate['M']='r'
    
    channelList=['F','M']
    for iNode in ['Root']:
        for load in channelList:
            for dir in ['x','y','z']:
                file.write(f'"{iNode}{load}{dir}{coordinate[load]}"\n')

    #Tip
    coordinate={}
    coordinate['RD']='r'
    coordinate['TD']='r'
    
    channelList=['RD','TD']
    for iNode in ['Tip']:
        for load in channelList:
            for dir in ['x','y','z']:
                file.write(f'"{iNode}{load}{dir}{coordinate[load]}"\n')
                

    file.write('END of input file (the word "END" must appear in the first 3 columns of this last OutList line)\n')
    file.write('---------------------------------------------------------------------------------------\n')

    file.close()

    print('Finished writing BeamDyn File')

    return axisFileName

# --- Write BeamDyn_Blade file with blade properties ---#
def write_beamdyn_prop(folder, wt_name, radial_stations, beam_stiff, beam_inertia, mu):
    n_pts = len(radial_stations)

    if not os.path.exists(folder):
        os.makedirs(folder)
        
    propFileName= wt_name + '_BeamDyn_Blade.dat'
    
    
    file = open(folder +'/'+propFileName, 'w')
    file.write(' ------- BEAMDYN V1.00.* INDIVIDUAL BLADE INPUT FILE --------------------------\n')
    file.write(' Test Format 1\n')
    file.write(' ---------------------- BLADE PARAMETERS --------------------------------------\n')
    file.write('%u   station_total    - Number of blade input stations (-)\n' % (n_pts))
    file.write(' 1   damp_type        - Damping type: 0: no damping; 1: damped\n')
    file.write('  ---------------------- DAMPING COEFFICIENT------------------------------------\n')
    file.write('   mu1        mu2        mu3        mu4        mu5        mu6\n')
    file.write('   (-)        (-)        (-)        (-)        (-)        (-)\n')
    file.write('\t %.5e \t %.5e \t %.5e \t %.5e \t %.5e \t %.5e\n' % (mu[0], mu[1], mu[2], mu[3], mu[4], mu[5])) 
    file.write(' ---------------------- DISTRIBUTED PROPERTIES---------------------------------\n')
    
    for i in range(n_pts):
        file.write('\t %.6f \n' % (radial_stations[i]))
        # write stiffness matrices
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_stiff[i, j, 0], beam_stiff[i, j, 1], beam_stiff[i, j, 2], beam_stiff[i, j, 3], beam_stiff[i, j, 4],
            beam_stiff[i, j, 5]))
        file.write('\n')

        # write inertia properties
        for j in range(6):
            file.write('\t %.16e \t %.16e \t %.16e \t %.16e \t %.16e \t %.16e\n' % (
            beam_inertia[i, j, 0], beam_inertia[i, j, 1], beam_inertia[i, j, 2], beam_inertia[i, j, 3],
            beam_inertia[i, j, 4], beam_inertia[i, j, 5]))
        file.write('\n')
        # ToDO: check correct translation of stiffness and mass matrices from VABS and anbax !!!
    file.close()

    print('Finished writing BeamDyn_Blade File')

    return propFileName


def writeBeamDynStandAlone(fileNames,disrLoads,tipLoads,directory='.'):

    if not os.path.exists(directory):
        os.makedirs(directory)

    from pynumad.utils.misc_utils import copy_and_replace

    templateFileName='beamDynStandAlone.template'
    
    analysisFileName=fileNames[0]+'_driver.inp'

    pathName=directory+'/'+analysisFileName



    
    copy_and_replace(templateFileName, pathName,
        {
            'DISTRLOAD1' : str(disrLoads[0]),
            'DISTRLOAD2' : str(disrLoads[1]),
            'DISTRLOAD3' : str(disrLoads[2]),
            'DISTRLOAD4' : str(disrLoads[3]),
            'DISTRLOAD5' : str(disrLoads[4]),
            'DISTRLOAD6' : str(disrLoads[5]),
            'TIPLOAD1' : str(tipLoads[0]),
            'TIPLOAD2' : str(tipLoads[1]),
            'TIPLOAD3' : str(tipLoads[2]),
            'TIPLOAD4' : str(tipLoads[3]),
            'TIPLOAD5' : str(tipLoads[4]),
            'TIPLOAD6' : str(tipLoads[5]),
            'AXIS FILE NAME': fileNames[0],
        })
    return analysisFileName

def runBeamDynStandAlone(beamDynDriverFileName,log,directory='.'):
    from pynumad import path_data
    import subprocess
    try:
        this_cmd = path_data['openFast']+'beamdyn_driver '+directory+'/'+beamDynDriverFileName
        log.info(f' running: {this_cmd}')
        subprocess.run(this_cmd, shell=True, check=True, capture_output=True)

        # with open(filePath+'.ech', 'r') as f:
        #     lines = f.readlines()
        # #log the last line of .ech file:
        # log.error(f'****************************\n{lines[-1]}\n******************************')

    except subprocess.CalledProcessError as e:
        log.error(f'Error running {this_cmd}: {e}')


        

