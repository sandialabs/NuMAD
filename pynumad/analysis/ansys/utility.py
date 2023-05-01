import numpy as np
import matplotlib.pyplot as plt
import numpy.matlib
from os.path import join
import numpy as np

from pynumad.objects.Subobjects import SkinArea
from pynumad.utils.fatigue import *
from pynumad.utils.interpolation import *


def getMatrialLayerInfoWithOutGUI(blade): 
    #Temparary workaround to extract data needed for:    
    #From write_shell7.m
    TotalStations = blade.station.size
    TotalShearwebs = blade.shearweb.size
    skin_areas = []
    for kStation in range(skin_areas.size):
        stationIB = blade.station(kStation)
        stationOB = blade.station(kStation + 1)
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
        for kdp in range(stationOB.dp.size-1):
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
    for k in range(TotalShearwebs):
        compsInModel = np.array([[compsInModel],[np.array([blade.shearweb[k].Material])]])
    
    #tcl: search skin materials
    for k in range(skin_areas.size):
        compsInModel = np.array([[compsInModel],[transpose(skin_areas[k].Material)]])
    
    compsInModel = np.unique(compsInModel)
    # load the material database and create searchable list
    if not len(blade.settings.job_name)==0 :
        # job name exists, use the local material database
        blade.matdb_path = fullfile(blade.settings.job_path,'MatDBsi.txt')
    else:
        # if no job name, use the master material database
        # jcb: we shouldn't arrive here because a file save is required first
        blade.matdb_path = fullfile(blade.numadpath,'MatDBsi.txt')
    
    blade.matdb = readMatDB(blade.matdb_path)
    for k in range(blade.matdb.size):
        blade.matlist[k] = blade.matdb[k].name
        mattype[k] = blade.matdb[k].type
    
    blade.isotropic = str('isotropic') == str(mattype)
    blade.orthotropic = str('orthotropic') == str(mattype)
    blade.composite = str('composite') == str(mattype)
    # Determine which isotropic and orthotropic materials are used in the model
    isoorthoInModel = np.array([])
    for kcomp in range(compsInModel.size):
        n = str(compsInModel(kcomp)) == str(blade.matlist)
        if not np.any(n) :
            raise Exception('Material "%s" not found in database.',compsInModel[kcomp])
        mat = blade.matdb[n]
        layerNames = []
        for klay in range(mat.layer.size):
            layerNames.bladeend(np.array([mat.layer(klay).layerName]))
        isoorthoInModel = np.unique(np.array([[isoorthoInModel],[layerNames]]))
    
    return isoorthoInModel,compsInModel,skin_areas,blade


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
    
    for kTheta in range(nDirections / 2+1):
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
                binnedElements = np.intersect(np.where(plateStrainsThetaSet[:,1] < z2),np.where(plateStrainsThetaSet[:,1] > z1))
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
            Pcr = - 3 / 2 * (2 * Dpt[0,0] * Ec * Gc) ** (1 / 3)
            Pphi = (R_sig[0,:] * S_alphaBeta) * h
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


