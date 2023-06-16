import numpy as np
import warnings
    
def calcFatigue(blade = None,
                meshData = None,
                IEC = None,
                Ltheta = None,
                LthetaPlus90 = None,
                Mtheta = None,
                MthetaPlus90 = None,
                binnedElements = None,
                plateStrainsTheta = None,
                plateStrainsThetaPlus90 = None,
                iSegment = None): 
    #blade.materials read in to obtain fatigue exponents
    #Initialize damage variables, get section data, and get element strains
    
    numElem = np.asarray(binnedElements).size
    if numElem == 0:
        if iSegment != 13:
            raise Exception('"zwidth" is smaller than element size. Decrease element size or increase zwdith')
        else:
            fatigueDamage = np.zeros((1,10))
            plotFatigue = []
            return fatigueDamage,plotFatigue
    
    analysis = 'LU'
    env = 'noEnv'
    temp = 'noTemp'
    mfg = 'basicFlaw'
    calca = 'test'
    loada = '12_directions'
    SFs = np.round(setIEC5(analysis,env,temp,mfg,calca,loada),2)
    #     analysis='LF';
#     loada='2_directions';
#     calcb ='measuredSlope';
#     loadb = 'Markov';
#     SFf=setIEC5(analysis,env,temp,mfg,calca,loada,calcb,loadb);
    
    SFf = SFs
    warnings.warn('SFf = SFs for BAR project since we are doing more than 2 load directions but less than 12')
    fatigueDamage = np.zeros((numElem,10))
    plotFatigue = []
    for i in np.arange(1,numElem+1).reshape(-1):
        elemFD = 0
        elemFDlayer = 0
        elemFDmat = 0
        elemFDflap = 0
        elemFDedge = 0
        elNo = plateStrainsTheta(binnedElements(i),1)
        coordSys = 'local'
        localFieldsTheta = extractFieldsThruThickness(plateStrainsTheta,meshData,blade.materials,blade.stacks,blade.swstacks,elNo,coordSys)
        localFieldsThetaPlus90 = extractFieldsThruThickness(plateStrainsThetaPlus90,meshData,blade.materials,blade.stacks,blade.swstacks,elNo,coordSys)
        npts = np.asarray(getattr(localFieldsTheta,(np.array(['element',num2str(elNo)]))).x3).size
        nptsPerLayer = 2
        layer = __builtint__.sorted(np.array([np.arange(1,npts / nptsPerLayer+1),np.arange(1,npts / nptsPerLayer+1)]))
        for ix3 in np.arange(1,npts+1).reshape(-1):
            matNumber = getattr(localFieldsTheta,(np.array(['element',num2str(elNo)]))).matNumber(ix3)
            if not len(blade.materials(matNumber).m)==0 :
                #Determine mean and amplitude stress
                MthetaFactor = getattr(localFieldsTheta,(np.array(['element',num2str(elNo)]))).sig11(ix3) / Mtheta
                MthetaPlus90Factor = getattr(localFieldsThetaPlus90,(np.array(['element',num2str(elNo)]))).sig11(ix3) / MthetaPlus90
                LthetaWithFactor = Ltheta
                LthetaPlus90WithFactor = LthetaPlus90
                LthetaWithFactor[1,:] = LthetaWithFactor[0,:] * MthetaFactor
                LthetaWithFactor[:,1] = np.abs(LthetaWithFactor[:,0] * MthetaFactor)
                LthetaPlus90WithFactor[1,:] = LthetaPlus90WithFactor[0,:] * MthetaPlus90Factor
                LthetaPlus90WithFactor[:,1] = np.abs(LthetaPlus90WithFactor[:,0] * MthetaPlus90Factor)
                m = blade.materials(matNumber).m
                XTEN = blade.materials(matNumber).uts(1) * SFs
                XCMP = blade.materials(matNumber).ucs(1) * SFs
                # Determine the maximum number of cycles for failure based
                # on available fatigue failure criterion or from data
                #Calculate fatigue damage value, layer, and material for flap and edge cycles
                if isfinite(MthetaFactor):
                    if 'Shifted Goodman' == IEC.fatigueCriterion:
                        layerFDtheta = shiftedGoodman(LthetaWithFactor,XTEN,XCMP,m,SFs,SFf)
                else:
                    layerFDtheta = 0
                if isfinite(MthetaPlus90Factor):
                    if 'Shifted Goodman' == IEC.fatigueCriterion:
                        layerFDthetaPlus90 = shiftedGoodman(LthetaPlus90WithFactor,XTEN,XCMP,m,SFs,SFf)
                else:
                    layerFDthetaPlus90 = 0
                layerFD = layerFDtheta + layerFDthetaPlus90
                if layerFD > elemFD:
                    elemFD = layerFD
                    elemFDlayer = layer(ix3)
                    elemFDmat = matNumber
                if layerFDtheta >= elemFDflap:
                    elemFDflap = layerFDtheta
                    elemFDflapLayer = layer(ix3)
                    elemFDflapMat = matNumber
                if layerFDthetaPlus90 >= elemFDedge:
                    elemFDedge = layerFDthetaPlus90
                    elemFDedgeLayer = layer(ix3)
                    elemFDedgeMat = matNumber
                FDvalue = np.array([elemFD,elemFDflap,elemFDedge])
                FDlayer = np.array([elemFDlayer,elemFDflapLayer,elemFDedgeLayer])
                FDmat = np.array([elemFDmat,elemFDflapMat,elemFDedgeMat])
                fatigueDamage[i,:] = np.array([elNo,FDvalue,FDlayer,FDmat])
        #         plotFatigue=[plotFatigue;elements(binnedElements(i),1) FDvalue(1)];
    
    #Output results in comma-delimited format
    #table(fatigue_damage(:,1),fatigue_damage(:,2),fatigue_damage(:,5),fatigue_damage(:,8))
    
    __,imax = np.amax(fatigueDamage[:,1])
    fatigueDamage = fatigueDamage[imax,:]
    
    
    #         #########################
    #         for j=1:nLayers #Loop through layers
    #             #Get material values
    #             matNumber = sections.layers{sections.secID==sec_num}(j,3);
    
    #             if ~isempty(blade.materials(matNumber).m)
    #                 #Determine mean and amplitude stress
    #                 MthetaFactor=stressesTheta(binnedElements(i),j+2)/Mtheta;
    #                 MthetaPlus90Factor=stressesThetaPlus90(binnedElements(i),j+2)/MthetaPlus90;
    
    #                 LthetaWithFactor=Ltheta; #Initialize for every layer
    #                 LthetaPlus90WithFactor=LthetaPlus90;
    
    #                 LthetaWithFactor(1,:)=LthetaWithFactor(1,:)*MthetaFactor; #Means
    #                 LthetaWithFactor(:,1)=abs(LthetaWithFactor(:,1)*MthetaFactor); #Amplitudes
    
    #                 LthetaPlus90WithFactor(1,:)=LthetaPlus90WithFactor(1,:)*MthetaPlus90Factor; #Means
    #                 LthetaPlus90WithFactor(:,1)=abs(LthetaPlus90WithFactor(:,1)*MthetaPlus90Factor); #Amplitudes
    
    #                 m=blade.materials(matNumber).m;
    #                 XTEN=materials(matNumber).XTEN*SFs; #Unfactor the factored resistance
    #                 XCMP=materials(matNumber).XCMP*SFs; #Unfactor the factored resistance
    
    
    #                 # Determine the maximum number of cycles for failure based
    #                 # on available fatigue failure criterion or from data
    
    #                 #Calculate fatigue damage value, layer, and material for flap and edge cycles
    #                 if isfinite(MthetaFactor)
    #                     switch IEC.fatigueCriterion
    #                         case 'Shifted Goodman'  #one case for now
    #                             layerFDtheta=shiftedGoodman(LthetaWithFactor,XTEN,XCMP,m,SFs,SFf);
    #                     end
    #                 else
    #                     layerFDtheta=0;
    #                 end
    
    #                 if isfinite(MthetaPlus90Factor)
    #                     switch IEC.fatigueCriterion
    #                         case 'Shifted Goodman' #one case for now
    #                             layerFDthetaPlus90=shiftedGoodman(LthetaPlus90WithFactor,XTEN,XCMP,m,SFs,SFf);
    #                     end
    #                 else
    #                     layerFDthetaPlus90=0;
    #                 end
    
    
    
    #                 layerFD=layerFDtheta+layerFDthetaPlus90;
    # #                 if elNo==1950
    # #                     keyboard
    # #                 end
    
    #                 if layerFD>elemFD
    #                     elemFD=layerFD;
    #                     elemFDlayer=j;
    #                     elemFDmat=matNumber;
    #                 end
    #                 if layerFDtheta>=elemFDflap
    #                     elemFDflap=layerFDtheta;
    #                     elemFDflapLayer=j;
    #                     elemFDflapMat=matNumber;
    #                 end
    #                 if layerFDthetaPlus90>=elemFDedge
    #                     elemFDedge=layerFDthetaPlus90;
    #                     elemFDedgeLayer=j;
    #                     elemFDedgeMat=matNumber;
    #                 end
    
    #                 FDvalue=[elemFD elemFDflap elemFDedge];
    #                 FDlayer=[elemFDlayer elemFDflapLayer elemFDedgeLayer];
    #                 FDmat=[elemFDmat elemFDflapMat elemFDedgeMat];
    #                 fatigueDamage(i,:)=[elNo FDvalue FDlayer FDmat];
    #            end
    #         end
    # #         plotFatigue=[plotFatigue;elements(binnedElements(i),1) FDvalue(1)];
    
    #     end
    #     #Output results in comma-delimited format
    #     #table(fatigue_damage(:,1),fatigue_damage(:,2),fatigue_damage(:,5),fatigue_damage(:,8))
    
    #     [~,imax]=max(fatigueDamage(:,2));
    #     fatigueDamage = fatigueDamage(imax,:);
    return fatigueDamage,plotFatigue

    
def getMomentMarkov(rccdata = None,
                    wt = None,
                    Yr = None,
                    simtime = None,
                    markovSize = None,
                    chSpan = None,
                    direction = None): 
    if chSpan == 1:
        baseStr = 'Root'
    else:
        baseStr = np.array(['Spn',int2str(chSpan - 1)])
    
    #Search through first windspeed column of rccdata for the first appearance of
    #baseStr + M + direction e.g. RootMyb1
    ct = 1
    
    # of baseStr + M + direction e.g. RootMyb1.
    while not (contains(rccdata[ct,1].label,baseStr) and contains(rccdata[ct,1].label,'M') and contains(rccdata[ct,1].label,direction)) :

        ct = ct + 1

    
    means = []
    
    ampl = []
    
    cycles = []
    for w in np.arange(1,rccdata.shape[2-1]+1).reshape(-1):
        # put data from rccdata structure for this channel and wind speed into a temporary variable, data
        data = rccdata[ct,w]
        # Make sure that fatigue data are only summed accross
        # windspeeds for the same channel.
        if w == 1 or str(data.label) == str(rccdata[ct,w - 1].label):
            means = np.array([[means],[data.means]])
            ampl = np.array([[ampl],[data.amplitudes]])
            cycles = np.array([[cycles],[data.cycles / simtime * (60 * 60 * 24 * 365.24 * Yr) * wt(w)]])
        else:
            raise Exception('Data channel name from current wind speed does not match the previous wind speed. Fatigue cycles cannot be summed.')
    
    Ni,EDGESi,BINi = histcounts(ampl,markovSize)
    
    Nj,EDGESj,BINj = histcounts(means,markovSize)
    
    #BIN, bin number assignment for each element in ampl or means
    
    markov = np.zeros((markovSize + 1,markovSize + 1))
    markov[np.arange[2,end()+1],1] = 0.5 * np.transpose((EDGESi(np.arange(1,end() - 1+1)) + EDGESi(np.arange(2,end()+1))))
    
    markov[1,np.arange[2,end()+1]] = 0.5 * (EDGESj(np.arange(1,end() - 1+1)) + EDGESj(np.arange(2,end()+1)))
    #indChecki=zeros(1,markovSize);
    for j in range(markovSize):
        ampIndecies = np.find(BINj == j)
        #indCheck=0;
        for i in range(markovSize):
            meanIndecies = np.find(BINi == i)
            Indecies = np.intersect(ampIndecies,meanIndecies)
            markov[i + 1,j + 1] = sum(cycles(Indecies))
    
    #surf(markov(2:end,2:end))
    return markov
    
    
def shiftedGoodman(markov, XTEN, XCMP, m, SFs, SFf): 
    num_range = markov.shape[0] - 1
    num_mean = markov.shape[1] - 1
    damage = 0
    FSloads = 1.0
    
    for i in range(num_range):
        for j in range(num_mean):
            sa = markov[i + 1,0]
            sm = markov[0,j + 1]
            n = markov(i + 1,j + 1)
            N = ((XTEN + np.abs(XCMP) - np.abs(2 * sm * SFs * FSloads - XTEN + np.abs(XCMP))) / (2 * sa * (SFf) * FSloads)) ** m
            damage = damage + n / N
    
    return damage