#import logging
import subprocess
#import os
import glob
import numpy as np


def writeBeamModel(wt_name,settings,blade,mu,log,directory='.'):
    import pynumad.analysis.beamUtils as beamUtils

#     #Runs VABS or OpenSG to homogenize
#     #Makes beamDyn or GEBT files



    radial_stations=blade.ispan/blade.ispan[-1]
    nStations=len(radial_stations)
    # #Run input files

    if 'vabs' in settings['solver'].lower():

        log.info(f'\n\n\nRunning VABS homogenization.')
        
        fileCount=0
        #First remove any lck files
        pattern=directory+'/'+wt_name+'*.in'
        if len(glob.glob(pattern))==0:
            raise RuntimeError(f'Could not find files with pattern: {pattern}. Beam model generation failed')
        MAXnLicenceTries=100
        for filePath in glob.glob(directory+'/'+wt_name+'*.in'):
            fileCount+=1
            try:
                this_cmd = 'VABS ' +filePath
                log.info(f' running: {this_cmd}')

                licenseAvailable=False
                nLicenceTries=0
                while not licenseAvailable and nLicenceTries <=MAXnLicenceTries-1:
                    subprocess.run(this_cmd, shell=True, check=True, capture_output=True)

                    with open(filePath+'.ech', 'r') as f:
                        lines = f.readlines()
                    #log the last line of .ech file:
                    
                    if 'Congratulations! No errors' in lines[-1]:
                        log.info(f'****************************\n{lines[-1]}\n******************************')
                        licenseAvailable=True
                        nLicenceTries=0
                    elif 'license' in lines[-1].lower():
                        nLicenceTries+=1
                        log.info(f'****************************\nnLicenceTries: {nLicenceTries}, {lines[-1]}\n******************************')

                    else:
                        log.error(f'****************************\n{lines[-1]}\n******************************')
                        raise Exception(f'Cross-sectional homogenization for file {filePath} failed due to: \n {lines[-1]} \n Beam model creation failed.') 
                if nLicenceTries ==MAXnLicenceTries:
                        string=f'License failed to be obtained after {MAXnLicenceTries} tries. Beam model creation failed.'
                        log.error(string)
                        raise Exception(string) 

            except subprocess.CalledProcessError as e:
                log.error(f'Error running {this_cmd}: {e}')
        
        # if fileCount != nStations:
        #     raise Exception('Error. Not enough VABS input files:')

    elif 'anba' in settings['solver'].lower():
        raise ValueError('ANBA currently not supported')



### Read inputs
    extension='K'

    blade.ispan=np.delete(blade.ispan,-2)  #TEMP Need to delete the added station near tip
    blade.idegreestwist=np.delete(blade.idegreestwist,-2)  #TEMP Need to delete the added station near tip

    radial_stations=blade.ispan/blade.ispan[-1]
    beam_stiff = np.zeros([len(radial_stations), 6, 6])
    beam_inertia = np.zeros([len(radial_stations), 6, 6])



    for fileName in glob.glob(directory+'/'+wt_name+"*." +extension):
        iStation=int(fileName.split('-')[-3].split('.')[0])
        print(f'fileName {fileName} iStation {iStation}')

        beam_stiff[iStation,:,:],beam_inertia[iStation,:,:]=beamUtils.readVABShomogenization(fileName)
    

    if 'beamdyn' in settings['beamSolver'].lower():
        beam_stiff,beam_inertia=beamUtils.transformMatrixToBeamDyn(beam_stiff,beam_inertia)
        axisFileName=beamUtils.write_beamdyn_axis(directory, wt_name, blade,radial_stations)
        propFileName=beamUtils.write_beamdyn_prop(directory, wt_name, radial_stations, beam_stiff, beam_inertia, mu)
    return [axisFileName,propFileName]