import numpy as np
import os
import warnings
import subprocess
from datetime import datetime
import numpy as np
import logging
import warnings
import os
from pynumad.utils.interpolation import calcGenLinePP
from pynumad.analysis.ansys.beamforce import *

from pynumad.analysis.ansys.read import readANSYSoutputs
from pynumad.analysis.ansys.utility import txt2mat, getLoadFactorsForElementsWithSameSection,\
    getMatrialLayerInfoWithOutGUI

def writeAnsysDeflections(blade, config, iLoad, fid, deflectionFilename): 
    # Outer AeroShell
    nStationLayups,nStations = blade.stacks.shape
    maxSectionNumber = int(str(nStations)+str(nStationLayups))
    
    # The following two lines help make unique IDs for web sections
    # based on the highes section already defined for aeroshell
    orderOfMagnitude = int(np.floor(np.log10(maxSectionNumber)))
    webSectionIDstart = np.ceil(maxSectionNumber / 10 ** orderOfMagnitude) * 10 ** orderOfMagnitude
    fid.write('/POST1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('RSYS,0\n' % ())
    
    fid.write('seltol,0.05\n' % ())
    for i in range(blade.ispan.size):
        fid.write('*CFOPEN, %s,out\n' % (deflectionFilename+'-'+str(i)))
        fid.write('ESEL,S,SEC,,1,%i   \n' % (webSectionIDstart))
        #fprintf(fid,'ESEL,S,SEC,,1,999   \n');    #Selects aero shell only
        fid.write('nsle,S,   \n' % ())
        fid.write('nsel,r,loc,z,%f  \n' % (blade.ispan[i]))
        #fprintf(fid,'nsll,s,,\n');
        if i == blade.ispan.size:
            fid.write('nsel,u,node,,z_master_node_number\n' % ())
        #fprintf(fid,'nplot\n');
        fid.write('*GET, NsectionNodes, NODE,0,COUNT   !Get the number of nodes in the set\n' % ())
        fid.write('*GET, node_num, NODE,0,NUM,MIN        !Get the smallest number node in the set\n' % ())
        fid.write('*DO, i, 1, NsectionNodes                 !loop through all nodes in cross section\n' % ())
        fid.write('*GET, xpos, NODE,node_num,loc,X\n' % ())
        fid.write('*GET, ypos, NODE,node_num,loc,Y\n' % ())
        fid.write('*GET, zpos, NODE,node_num,loc,Z\n' % ())
        fid.write('*GET, u1, NODE,node_num,U,X\n' % ())
        fid.write('*GET, u2, NODE,node_num,U,Y\n' % ())
        fid.write('*GET, u3, NODE,node_num,U,Z\n' % ())
        fid.write(' *VWRITE,node_num,xpos,ypos,zpos,u1,u2,u3\n' % ())
        fid.write('(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n' % ())
        fid.write('node_num=NDNEXT(node_num)             !Get the next higher node number in the set\n' % ())
        fid.write('*ENDDO\n' % ())
        fid.write('*CFCLOS\n' % ())
        fid.write('\n \n \n' % ())
    
    fid.write('finish\n' % ())
    return


def writeAnsysFagerberWrinkling(app, SkinAreas, compsInModel, coreMatName): 
    #limitingElementData - [ansysSecNumber elno lf phicr]
    TotalStations = app.station.size
    TotalShearwebs = app.shearweb.size
    #################   Main loop #1: loop around aero shell.   #################
    LF = []
    for kStation in range(TotalStations-1):
        for kArea in range(SkinAreas[kStation].startIB.size):
            #See if the section contatins Balsa/core material name (i.e find
            #the sandwhich panels)
            n = str(SkinAreas[kStation].Material[kArea]) == str(app.matlist)
            mat = app.matdb(n)
            if coreMatName in mat.layer.layerName:
                ansysSecNumber = np.where(str(SkinAreas[kStation].Material[kArea]) == str(compsInModel))
                file = 'section-'+str(ansysSecNumber)+'-faceAvgStresses.txt'
                avgFaceStress = txt2mat(file)
                os.delete(file)
                LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName)
    
    #################   Main loop #2: loop along web.   #################
    for kShearweb in range(TotalShearwebs):
        n = str(app.shearweb(kShearweb).Material) == str(app.matlist)
        mat = app.matdb(n)
        if coreMatName in mat.layer.layerName:
            ansysSecNumber = np.where(str(app.shearweb(kShearweb).Material) == str(compsInModel))
            file = 'section-'+str(ansysSecNumber + 1000)+'-faceAvgStresses.txt'
            avgFaceStress = txt2mat(file)
            os.delete(file)
            LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber + 1000,avgFaceStress,app,mat,coreMatName)
    
    minLF,index = np.amin(LF[:,2])
    limitingElementData = LF[index,:]
    print('\n\n The minimum wrinkling LF is: %f, wrinkle angle: %.2f°' % (minLF,LF[index,3]))
    print('\n and occurs in section number %i, element number %i\n, ' % (LF[index,0],LF[index,1]))
    # [maxLF,index]=max(LF(:,3));
    # fprintf('\n\n The maximum LF is: #f, wrinkle  angle: #.2f°' ,maxLF, LF(index,4))
    # fprintf('\n and occurs in section number #i, element number #i, ',LF(index,1),LF(index,2))
    
    return limitingElementData


def writeAnsysFatigue(fid, iLoad): 
    ###################Outputs for fatigue analysis in MATLAB#################
    fid.write('! BEGIN FATIGUE SCRIPT\n' % ())
    fid.write('allsel\n' % ())
    fid.write('/prep7\n' % ())
    fid.write('esel,all\n' % ())
    fid.write('allsel\n' % ())
    fid.write('/prep7\n' % ())
    fid.write('esel,all\n' % ())
    fid.write('esel,u,type,,21  \n' % ())
    fid.write('/POST1\n' % ())
    fid.write('set,LAST\n' % ())
    fid.write('RSYS,SOLU\n' % ())
    
    ### Element strains and curvatures ###
    fid.write('ALLSEL\n' % ())
    fid.write('ETABLE, zcent,CENT,Z\n' % ())
    fid.write('ETABLE, eps11,SMISC,9 \n' % ())
    fid.write('ETABLE, eps22,SMISC,10 \n' % ())
    fid.write('ETABLE, eps12,SMISC,11 \n' % ())
    fid.write('ETABLE, kapa11,SMISC,12 \n' % ())
    fid.write('ETABLE, kapa22,SMISC,13 \n' % ())
    fid.write('ETABLE, kapa12,SMISC,14 \n' % ())
    fid.write('ETABLE, gamma13,SMISC,15 \n' % ())
    fid.write('ETABLE, gamma23,SMISC,16 \n' % ())
    fid.write('/output,plateStrains-all-%s,txt\n' % (str(iLoad)))
    fid.write('PRETAB,zcent,eps11,eps22,eps12,kapa11,kapa22,kapa12,gamma12,gamma13,gamma23\n' % ())
    fid.write('ETABLE,ERAS\n\n' % ())
    fid.write('finish\n' % ())
    fid.write('! END FATIGUE OUTPUT SCRIPT\n' % ())
    return

 
def writeAnsysGetFaceStresses(blade, fid, coreMatName): 
    isoorthoInModel,compsInModel,SkinAreas,app = getMatrialLayerInfoWithOutGUI(blade)
    #fid=fopen('getFaceStresses.mac','w+');
    TotalStations = blade.ispan.size
    for kStation in range(TotalStations-1):
        #kPanel=find(~cellfun('isempty',strfind([SkinAreas(kStation).Material],'PANEL'))); #Array that stores the kArea index that contains 'PANEL' in the name
        #for i=1:numel(kPanel)
        for kArea in range(SkinAreas[kStation].startIB.size):
            #See if the section contatins Balsa/core material name (i.e find
            #the sandwhich panels)
            n = str(SkinAreas[kStation].Material[kArea]) == str(app.matlist)
            mat = app.matdb[n]
            if coreMatName in mat.layer.layerName:
                ansysSecNumber = np.where(str(SkinAreas[kStation].Material[kArea]) == str(compsInModel))
                writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)
    
    fid.write('\n' % ())
    fid.write('\n' % ())
    fid.write('!*************** WEB ***************\n' % ())
    fid.write('\n' % ())
    fid.write('\n' % ())
    fid.close()
    #Web
    TotalShearwebs = np.asarray(app.shearweb).size
    for kShearweb in range(TotalShearwebs):
        n = str(app.shearweb[kShearweb].Material) == str(app.matlist)
        mat = app.matdb[n]
        if coreMatName in mat.layer.layerName:
            ansysSecNumber = np.where(str(app.shearweb[kShearweb].Material) == str(compsInModel))
            ansysSecNumber = ansysSecNumber + 1000
            writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)
    
    fid.write('FINISH\n' % ())
    fid.write('allsel\n' % ())
    return app,SkinAreas,compsInModel


def writeANSYSinputFile(fid, mat, ansysSecNumber, coreMatName): 
    #####Find the face sheet####
    cellMat = np.array([])
    for i in range(len(mat.layer)):
        cellMat = np.array([[cellMat],[np.array([mat.layer(i).layerName])]])
    
    kbalsa = np.where(str(coreMatName) == str(cellMat))
    iLayer = np.arange(0,(kbalsa - 1))
    
    # Find the number of layers in the face
    qty = 0
    
    for i in range(iLayer.size):
        qty = qty + mat.layer(iLayer[i]).quantity
    
    #Loop through the top facesheet layers
    
    fid.write('!*************** ansysSecNumber = %i ***************\n' % (ansysSecNumber))
    fid.write('/POST1\n' % ())
    fid.write('*DEL,iel\n' % ())
    fid.write('*DEL,enum\n' % ())
    fid.write('*DEL,nelTemp\n' % ())
    fid.write('RSYS, SOLU\n' % ())
    fid.write('ALLSEL\n' % ())
    fid.write('ESEL, S, SEC,,%i\n' % (ansysSecNumber))
    fid.write('*GET, enum, ELEM, 0, NUM, MIN, !  lowest element number in the selected set\n' % ())
    fid.write('*get, nelTemp, ELEM,0,count\n' % ())
    fid.write('*DIM, iel,ARRAY,nelTemp\n' % ())
    fname = 'section-'+str(ansysSecNumber)+'-faceAvgStresses'
    
    if os.path.isfile(fname+'.txt'):
        os.delete(fname+'.txt')
    
    fid.write('*CFOPEN, %s, txt,,APPEND\n' % (fname))
    #Create an array with the element numbers in the selected set
    fid.write('*DO, J, 1,nelTemp  !Loop through elements\n' % ())
    fid.write('iel(J)=enum\n' % ())
    fid.write('enum =ELNEXT(enum)  !Next higher element number above N in selected set\n' % ())
    fid.write('*ENDDO\n' % ())
    fid.write('\n' % ())
    fid.write('ALLSEL\n' % ())
    fid.write('*DO, J, 1,nelTemp  !Loop through elements\n' % ())
    fid.write('	S11a=0 !Initialize average stress variables for each element\n' % ())
    fid.write('	S22a=0\n' % ())
    fid.write('	S33a=0\n' % ())
    fid.write('	S23a=0\n' % ())
    fid.write('	S13a=0\n' % ())
    fid.write('	S12a=0\n' % ())
    fid.write('	*DO, I, 1,%i    !Loop through face layers\n' % (qty))
    fid.write('	    LAYER,I\n' % ())
    fid.write('		SHELL,MID   !Stress result at midlayer\n' % ())
    fid.write('		ESEL,S,ELEM,,iel(J)\n' % ())
    fid.write('		ETABLE,ERAS !Each element gets a new element table\n' % ())
    fid.write('	    ETABLE,S11,S,X,AVG !AVG - Store averaged element centroid value\n' % ())
    fid.write('	    ETABLE,S22,S,Y,AVG\n' % ())
    fid.write('		ETABLE,S33,S,Z,AVG\n' % ())
    fid.write('		ETABLE,S23,S,YZ,AVG\n' % ())
    fid.write('		ETABLE,S13,S,XZ,AVG\n' % ())
    fid.write('        ETABLE,S12,S,XY,AVG\n' % ())
    fid.write('	    *GET,tempS11, ELEM, iel(J), ETAB, S11\n' % ())
    fid.write('		*GET,tempS22, ELEM, iel(J), ETAB, S22\n' % ())
    fid.write('		*GET,tempS33, ELEM, iel(J), ETAB, S33\n' % ())
    fid.write('		*GET,tempS23, ELEM, iel(J), ETAB, S23\n' % ())
    fid.write('		*GET,tempS13, ELEM, iel(J), ETAB, S13\n' % ())
    fid.write('		*GET,tempS12, ELEM, iel(J), ETAB, S12\n' % ())
    fid.write('		S11a=S11a+tempS11\n' % ())
    fid.write('		S22a=S22a+tempS22\n' % ())
    fid.write('		S33a=S33a+tempS33\n' % ())
    fid.write('		S23a=S23a+tempS23\n' % ())
    fid.write('		S13a=S13a+tempS13\n' % ())
    fid.write('		S12a=S12a+tempS12\n' % ())
    fid.write('	*ENDDO\n' % ())
    fid.write('	S11a=S11a/%i\n' % (qty))
    fid.write('	S22a=S22a/%i\n' % (qty))
    fid.write('	S33a=S33a/%i\n' % (qty))
    fid.write('	S23a=S23a/%i\n' % (qty))
    fid.write('	S13a=S13a/%i\n' % (qty))
    fid.write('	S12a=S12a/%i\n' % (qty))
    fid.write('	ELNO=iel(J)  !It is needed to refer to ELNO in the command below\n' % ())
    fid.write('*VWRITE,ELNO,S11a,S22a,S33a,S23a,S13a,S12a\n' % ())
    fid.write('(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n' % ())
    fid.write('*ENDDO\n' % ())
    fid.write('*CFCLOS\n' % ())
    fid.write('\n' % ())
    fid.write('\n' % ())
    return
    

def writeAnsysLinearBuckling(blade, config, iLoad, fid, bucklingFilename): 
    fid.write('! BEGIN BUCKLE MACRO\n' % ())
    fid.write('allsel\n' % ())
    fid.write('/solu\n' % ())
    fid.write('irlf,-1\n' % ())
    fid.write('pstres,on\n' % ())
    fid.write('antype,buckle\n' % ())
    fid.write('bucopt,lanb,'+str(config.analysisFlags.globalBuckling)+',,,RANGE\n') % ()
    #fprintf(fid,strcat('MXPAND,',int2str(nmodes),',0,0,1\n'), nmodes); # Required for element stress/strain, etc..
    fid.write('solve\n' % ())
    fid.write('finish\n' % ())
    fid.write('/post1\n' % ())
    fid.write(np.array(['/output,',bucklingFilename,',out\n']) % ())
    fid.write('set,list\n' % ())
    fid.write('/output\n' % ())
    fid.write('finish\n' % ())
    fid.write('! END BUCKLE MACRO\n' % ())
    return


def writeAnsysLocalFields(blade, config, iLoad, fid
                          ): 
    ###################Outputs for fatigue analysis in MATLAB#################
    fid.write('! BEGIN LOCAL FIELD SCRIPT\n' % ())
    fid.write('allsel\n' % ())
    fid.write('/post1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('esel,all\n' % ())
    ### Element Stress ###
    fid.write('ALLSEL\n' % ())
    fid.write('ETABLE, zcent,CENT,Z\n' % ())
    fid.write('ETABLE, eps11,SMISC,9 \n' % ())
    fid.write('ETABLE, eps22,SMISC,10 \n' % ())
    fid.write('ETABLE, eps12,SMISC,11 \n' % ())
    fid.write('ETABLE, kapa11,SMISC,12 \n' % ())
    fid.write('ETABLE, kapa22,SMISC,13 \n' % ())
    fid.write('ETABLE, kapa12,SMISC,14 \n' % ())
    fid.write('ETABLE, gamma13,SMISC,15 \n' % ())
    fid.write('ETABLE, gamma23,SMISC,16 \n' % ())
    fid.write('/output,plateStrains-all-%s,txt\n' % (str(iLoad)))
    fid.write('PRETAB,zcent,eps11,eps22,eps12,kapa11,kapa22,kapa12,gamma12,gamma13,gamma23\n' % ())
    fid.write('ETABLE,ERAS\n\n' % ())
    #fprintf(fid,'ETABLE, zcent,CENT,Z\n');
    #fprintf(fid,'ETABLE, N11,SMISC,1 \n');
    #fprintf(fid,'ETABLE, N22,SMISC,2 \n');
    #fprintf(fid,'ETABLE, N12,SMISC,3 \n');
    #fprintf(fid,'ETABLE, M11,SMISC,4 \n');
    #fprintf(fid,'ETABLE, M22,SMISC,5 \n');
    #fprintf(fid,'ETABLE, M12,SMISC,6 \n');
    #fprintf(fid,'ETABLE, Q13,SMISC,7 \n');
    #fprintf(fid,'ETABLE, Q23,SMISC,8 \n');
    #fprintf(fid,'/output,plateExamplePlateForces-all-#s,txt\n',int2str(iLoad));
    #fprintf(fid,'/output,plateForces-all-#s,txt\n',int2str(iLoad));
    #fprintf(fid,'PRETAB,zcent,N11,N22,N12,M11,M22,M12,Q12,Q13,Q23\n');
    #fprintf(fid, 'ETABLE,ERAS\n\n');
    fid.write('finish\n' % ())
    return


def writeAnsysNonLinearBuckling(ansysFilename, ansys_path, ansys_product, config, ii, jj, ncpus, iLoad): 
    warnings.warn('output designvar. Currently does not work for nonlinear cases')
    script_name = 'commands3-'+str(ii)+'.mac'
    script_out = 'output3-'+str(ii)+'-'+str(jj)+'.txt'
    fid = open(script_name,'w+')
    fid.write('!************   MODE-%i   ************\n' % (ii))
    fid.write('/FILNAME,%s,1\n' % ansysFilename + '-Load' + str(iLoad))
    
    fid.write('resume, %s,db\n' % (ansysFilename+'-Load'+str(iLoad)))
    #Get Max displacement, UY
    #     fprintf(fid,'/POST1\n');
    #     fprintf(fid,'SET,1,#i\n',ii); #Read in results
    #     fprintf(fid,'nsel, all, node\n'); #Select all nodes
    #     fprintf(fid,'*get, Zncount, node,0,count\n'); #Find the number (quantity) of nodes selected
    #     fprintf(fid,'*dim,zNodeDisp,array,Zncount,1 \n'); #Allocate memory for an arry to hold nodal disp.
    
    #     #For each node populate array with Y-displacements
    #     fprintf(fid,'*DO, i,1,Zncount\n');
    #     fprintf(fid,'*VGET, zNodeDisp(i,1), NODE, i, U, Y\n');
    #     fprintf(fid,'*ENDDO\n');
    
    #     #Find the min/max disp. value
    #     fprintf(fid,'*VSCFUN,zMaxUY,max,zNodeDisp\n');
    #     fprintf(fid,'*VSCFUN,zMinUY,min,zNodeDisp\n');
    #     U0 = 1/400; #Dimple factor
    #     lg = 1;     #Largest horizontal dimension of buckle
    #     fprintf(fid,'zImperfectionSF=#f*#f*#f/max(abs(zMaxUY),abs(zMinUY))\n',config.analysisFlags.imperfection(jj), U0, lg);
    
    fid.write('zImperfectionSF=%f\n' % (config.analysisFlags.imperfection(jj)))
    #     U0 = 1/400; #Dimple factor
    #     lg = 1;     #Largest horizontal dimension of buckle
    #     zImperfectionSF = lg*U0;
    
    fid.write('/prep7\n' % ())
    fid.write('UPGEOM,zImperfectionSF,1,%i,%s,"rst"\n' % (ii, ansysFilename+'-Load'+str(iLoad)))
    
    fid.write('FINISH\n' % ())
    
    filename = ansysFilename+'-'+str(ii)+'-'+str(jj)
    
    fid.write('/FILNAME,%s,1\n' % (filename))
    
    #fprintf(fid,strcat('SAVE,''',filename,''',''db'',''',strrep(pwd,'\','\\'),'''\n'));
    
    fid.write('\n' % ())
    #Nonlinear Static Analysis
    fid.write('/solu\n' % ())
    fid.write('antype,0\n' % ())
    #fprintf(fid,'irlf,-1\n');
    fid.write('pstres,0\n' % ())
    fid.write('NLGEOM,1\n' % ())
    fid.write('TIME,%f\n' % (1))
    #     fprintf(fid,'AUTOTS,ON,\n');
    #     nsubstep = 20;
    #     fprintf(fid,'nsubstep=#f\n', nsubstep);
    #     fprintf(fid,'NSUBST,508,20,500\n');
    #     fprintf(fid,'NEQIT,200,\n');
    loadScaleFactor = 5
    
    #     fprintf(fid,'allsel\n');
    #     fprintf(fid,'esel,s,type,,33\n'); #Select all follower elements
    fid.write('NSEL, ALL\n' % ())
    fid.write('FSCALE,%f\n' % (loadScaleFactor))
    #     fprintf(fid,'RESCONTROL,DEFINE,ALL,1,\n');
    fid.write('OUTRES,NSOL,ALL\n' % ())
    #     fprintf(fid,'NROPT,UNSYM\n');
    
    #     fprintf(fid,'CUTCONTROL,PIVSTOP,2\n'); #Ends simulation once pivot becomes negative
    #     fprintf(fid,'PRED,OFF\n');
    
    fid.write('allsel\n' % ())
    fid.write('solve\n' % ())
    fid.write('FINISH\n' % ())
    # fid.write('SAVE,'',filename,'','db','',pwd,'\n') % ())
    fid.write('/EXIT,NOSAVE\n' % ())
    fid.close()
    ###### RUN ANSYS#######
    ansys_call = print('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',ansys_path,ansys_product,script_name,script_out,int2str(ncpus))
    # KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
    # simulations. EC
    
    subprocess.run(ansys_call)
    print('%s: Nonlinear Mode-%s Analysis Finished\n' % (datetime.now(),str(ii)))
    data = readANSYSoutputs(filename+'.mntr',11)
    a = data.shape
    nonlinearLoadFactors = data(a(1),7) * loadScaleFactor
    
    return nonlinearLoadFactors


def writeAnsysResultantVSSpan(blade, config, iLoad, fid): 
    fid.write('/POST1\n' % ())
    fid.write('set,LAST\n' % ())
    fid.write('RSYS,0\n' % ())
    
    #fprintf(fid,'seltol,0.05\n');
    #fprintf(fid,'*CFOPEN, resultantVSspan,txt\n');
    #for i=1:numel(blade.ispan)
    #fprintf(fid,'nsel,s,loc,z,0,#f  \n',blade.ispan(i));
    
    #if i==numel(blade.ispan)
    #fprintf(fid,'nsel,u,node,,z_master_node_number\n');
    #end
    
    #fprintf(fid,'spoint,0,#f,#f,#f\n',blade.sweep(i),blade.prebend(i),blade.ispan(i));
    #fprintf(fid,'nplot\n');
    #fprintf(fid,'FSUM\n');
    #fprintf(fid,'*GET, F1, FSUM, 0, ITEM,FX\n');
    #fprintf(fid,'*GET, F2, FSUM, 0, ITEM,FY\n');
    #fprintf(fid,'*GET, F3, FSUM, 0, ITEM,FZ\n');
    #fprintf(fid,'*GET, M1, FSUM, 0, ITEM,MX\n');
    #fprintf(fid,'*GET, M2, FSUM, 0, ITEM,MY\n');
    #fprintf(fid,'*GET, M3, FSUM, 0, ITEM,MZ\n');
    #fprintf(fid,'*VWRITE,#f,F1,F2,F3,M1,M2,M3\n',blade.ispan(i));
    #fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
    #fprintf(fid,'\n \n \n');
    #end
    #fprintf(fid,'*CFCLOS\n');
    #fprintf(fid,'finish\n');
    
    fid.write('/post1\n' % ())
    fid.write('elsize=%f\n' % (blade.mesh))
    fid.write('nz=nint(%f/elsize) !Integer number of points to output resultant loads\n' % (blade.ispan(end())))
    fid.write('zloc=0\n' % ())
    fid.write('delta=0.1\n' % ())
    fid.write('*CFOPEN, resultantVSspan,txt\n' % ())
    fid.write('*do,I,1,nz+1\n' % ())
    fid.write('allsel\n' % ())
    fid.write('nsel,s,loc,z,0,zloc+delta\n' % ())
    fid.write('spoint,0,0,0,zloc\n' % ())
    fid.write('!nplot\n' % ())
    fid.write('FSUM\n' % ())
    fid.write('*GET, F1, FSUM, 0, ITEM,FX\n' % ())
    fid.write('*GET, F2, FSUM, 0, ITEM,FY\n' % ())
    fid.write('*GET, F3, FSUM, 0, ITEM,FZ\n' % ())
    fid.write('*GET, M1, FSUM, 0, ITEM,MX\n' % ())
    fid.write('*GET, M2, FSUM, 0, ITEM,MY\n' % ())
    fid.write('*GET, M3, FSUM, 0, ITEM,MZ\n' % ())
    fid.write('*VWRITE,zloc,F1,F2,F3,M1,M2,M3\n' % ())
    fid.write('(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n' % ())
    fid.write('zloc=zloc+elsize\n' % ())
    fid.write('*ENDDO\n' % ())
    fid.write('*CFCLOS\n' % ())
    fid.write('finish\n' % ())
    return


def writeAnsysRupture(config, iLoad, fid, failureFilename): 
    fid.write('! BEGIN FAILURE SCRIPT\n' % ())
    fid.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' % ())
    fid.write('!Add for PLESOL and *get,findex,PLNSOL,0,MAX to work' % ())
    fid.write('/BATCH  \n' % ())
    fid.write('/COM,ANSYS RELEASE Release 18.1      BUILD 18.1      UP20170403       15:49:08\n' % ())
    fid.write('/GRA,POWER\n ' % ())
    fid.write('/GST,ON\n ' % ())
    fid.write('/PLO,INFO,3\n ' % ())
    fid.write('/GRO,CURL,ON\n ' % ())
    fid.write('/CPLANE,1   \n ' % ())
    fid.write('/REPLOT,RESIZE  \n ' % ())
    fid.write('WPSTYLE,,,,,,,,0\n ' % ())
    fid.write('/SHOW\n ' % ())
    fid.write('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n' % ())
    fid.write('/POST1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('allsel\n' % ())
    fid.write('RSYS,LSYS \n' % ())
    
    fid.write('layer,fcmax\n' % ())
    if not config.analysisFlags.failure.upper() in ['PUCK','LARC03','LARC04']:
        #Do this for the failure criteria that do not distinguish between fiber
        #and matrix failure
        fc = config.analysisFlags.failure.upper()
        fid.write('FCTYP,add,%s\n' % (fc))
        fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
        fid.write('*get,findex,PLNSOL,0,MAX\n' % ())
    else:
        #Do this for the failure criteria that do distinguish between fiber
        #and matrix failure
        if 'PUCK' == config.analysisFlags.failure.upper():
            #Fiber Failure
            fc = 'PFIB'
            fid.write('FCTYP,add,%s\n' % (fc))
            fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
            fid.write('*get,Ffindex,PLNSOL,0,MAX\n' % ())
            #Matrix Failure
            fc = 'PMAT'
            fid.write('FCTYP,add,%s\n' % (fc))
            fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
            fid.write('*get,Mfindex,PLNSOL,0,MAX\n' % ())
        else:
            if 'LARC03' == config.analysisFlags.failure.upper():
                #Fiber Failure
                fc = 'L3FB'
                fid.write('FCTYP,add,%s\n' % (fc))
                fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
                fid.write('*get,Ffindex,PLNSOL,0,MAX\n' % ())
                #Matrix Failure
                fc = 'L3MT'
                fid.write('FCTYP,add,%s\n' % (fc))
                fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
                fid.write('*get,Mfindex,PLNSOL,0,MAX\n' % ())
            else:
                if 'LARC04' == config.analysisFlags.failure.upper():
                    #Fiber Failure
                    fc = 'L4FB'
                    fid.write('FCTYP,add,%s\n' % (fc))
                    fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
                    fid.write('*get,Ffindex,PLNSOL,0,MAX\n' % ())
                    #Matrix Failure
                    fc = 'L4MT'
                    fid.write('FCTYP,add,%s\n' % (fc))
                    fid.write('PLESOL, FAIL,%s, 0,1.0\n' % (fc))
                    fid.write('*get,Mfindex,PLNSOL,0,MAX\n' % ())
        #Report the higher of the fiber failure index or the matrix
        fid.write('*IF, Ffindex, GT,Mfindex, THEN\n' % ())
        fid.write('findex=Ffindex\n' % ())
        fid.write('*ELSE\n' % ())
        fid.write('findex=Mfindex\n' % ())
        fid.write('*ENDIF\n' % ())
    
    fid.write(np.array(['/output,',failureFilename,',out\n']) % ())
    fid.write('*status,findex\n' % ())
    fid.write('/output\n' % ())
    ## EMA added:
    fid.write('/output,allElemFailureResults%s,out\n' % (str(iLoad)))
    fid.write('PRESOL,FAIL\n' % ())
    fid.write('/output\n' % ())
    ## END
    fid.write('finish\n' % ())
    fid.write('! END FAILURE SCRIPT\n' % ())
    return


def writeWrinklingForNonlinearBuckling(blade, coreMatName, settings, np, ansysFilename, i, j): 
    filename = ansysFilename+'-'+str(i)+'-'+str(j)
    
    #######        Generate Wrinkling Files        #############
    script_name = 'commands4-'+str(i)+'.mac'
    script_out = 'output4-'+str(i)+'-'+str(j)+'.txt'
    fid = open(script_name,'w+')
    fid.write('!************   MODE-%i   ************\n' % (i))
    fid.write('/FILNAME,%s,1\n' % (filename))
    
    fid.write('resume\n' % ())
    fid.write('/POST1\n' % ())
    fid.write('SET,LAST\n' % ())
    app,SkinAreas,compsInModel = writeAnsysGetFaceStresses(blade,fid,coreMatName)
    fid.write('/EXIT,NOSAVE\n' % ())
    fid.close()
    ansys_call = print('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',settings.ansys_path,settings.ansys_product,script_name,script_out,int2str(np))
    # KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
    # simulations. EC
    
    subprocess.run(ansys_call)
    print('%s: Nonlinear Mode-%s Analysis Finished\n' % (datetime.now(),str(i)))
    wrinklingLimitingElementData = writeAnsysFagerberWrinkling(app,SkinAreas,compsInModel,coreMatName)
    return wrinklingLimitingElementData
    

def writeAnsysNonLinearLocalBuckling(blade, config, iLoad, fid, ansysFilename, ii, jj): 
    #UNSUPPORTED AT THIS TIME
    #                 filename=strcat(ansysFilename,'-',int2str(ii),'-',int2str(jj)); #The name of the next job name
    #                 #######        Generate Wrinkling Files        #############
    #                 script_name=strcat('commands4-',int2str(ii),'.mac');
    #                 script_out=strcat('output4-',int2str(ii),'-',int2str(jj),'.txt');
    
    #                 fid=fopen(script_name,'w+');
    #                 fprintf(fid,strcat('!************   MODE-#i   ************\n'),ii);
    #                 fprintf(fid,'/FILNAME,''#s'',1\n',filename);   #From master, change the jobname
    #                 fprintf(fid,'resume\n');
    #                 fprintf(fid,'/POST1\n');
    #                 fprintf(fid,'SET,LAST\n');
    
    #                 [app,SkinAreas,compsInModel]=writeANSYSgetFaceStresses(blade,fid,config.analysisFlags.localBuckling);
    
    #                 fprintf(fid,'/EXIT,NOSAVE\n');
    # fid.close();
    
    #                 ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "#s" -b -p #s -I #s -o #s -np #s',settings.ansys_path,settings.ansys_product,script_name,script_out,int2str(np))    # KMP_STACKSIZE is 512k by default. This is not enough therefore SET
    #                 # KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
    #                 # simulations. EC
    #                 #
    
    
    #                 system(ansys_call)  # the windows system call to run the above ansys command
    #                 fprintf('#s: Nonlinear Mode-#s Analysis Finished\n',datestr(now),int2str(ii))
    return


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
    fid = open(fname+'.mac','w')
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


def writeforcefile(filename, forcemap, forcesums, maptype): 
    fid = open(filename,'wt')
    if (fid == - 1):
        raise Exception('Could not open file "%s"',filename)
    
    if ('forcesums' is not None):
        fid.write('!========== FORCE MAPPING SUMMARY ==========' % ())
        fid.write('\n!maptype = "%s"' % (maptype))
        s = '  %14.6e'*forcesums["Z"].size % tuple(forcesums["Z"])
        s = '\n!                Z =%s'% s +'      TOTAL     '
        fid.write(s)
        fid.write('\n!'+('-'*len(s)))
        fid.write('\n!Input          Fx =' % ())
        fid.write('  %14.6e'*forcesums["Fx"][:,0].size % tuple(forcesums["Fx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,0])))
        fid.write('\n!Output    sum(fx) =' % ())
        fid.write('  %14.6e'*forcesums["Fx"][:,1].size % tuple(forcesums["Fx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fx"][:,1])))
        fid.write('\n!'+('-'*len(s)))
        fid.write('\n!Input          Fy =' % ())
        fid.write('  %14.6e'*forcesums["Fy"][:,0].size % tuple(forcesums["Fy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,0])))
        fid.write('\n!Output    sum(fy) =' % ())
        fid.write('  %14.6e'*forcesums["Fx"][:,1].size % tuple(forcesums["Fx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["Fy"][:,1])))
        fid.write('\n!'+('-'*len(s)))
        fid.write('\n!Input           M =' % ())
        fid.write('  %14.6e'*forcesums["M"][:,0].size % tuple(forcesums["M"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,0])))
        fid.write('\n!sum(-y*fx + x*fy) =' % ())
        fid.write('  %14.6e'*forcesums["M"][:,1].size % tuple(forcesums["M"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["M"][:,1])))
        fid.write('\n!'+('-'*len(s)))
        fid.write('\n!Input        Z*Fy =' % ())
        fid.write('  %14.6e'*forcesums["RootMx"][:,0].size % tuple(forcesums["RootMx"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,0])))
        fid.write('\n!Output  sum(z*fy) =' % ())
        fid.write('  %14.6e'*forcesums["RootMx"][:,1].size % tuple(forcesums["RootMx"][:,1]))
        fid.write('  %14.6e' % (sum(forcesums["RootMx"][:,1])))
        fid.write('\n!'+('-'*len(s)))
        fid.write('\n!Input        Z*Fx =' % ())
        fid.write('  %14.6e'*forcesums["RootMy"][:,0].size % tuple(forcesums["RootMy"][:,0]))
        fid.write('  %14.6e' % (sum(forcesums["RootMy"][:,0])))
        fid.write('\n!Output  sum(z*fx) =' % ())
        fid.write('  %14.6e'*forcesums["RootMy"][:,1].size % tuple(forcesums["RootMy"][:,1]))
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


def writeAnsysShellModel(blade, filename, meshData, config):
    """ WRITE_SHELL7 Generate the ANSYS input file that creates the blade

    Parameters
    ----------

    Results
    -------

    """
    
    fcopts = [
        'EMAX','SMAX','TWSI','TWSR','HFIB',
        'HMAT','PFIB','PMAT','L3FB','L3MT',
        'L4FB','L4MT','USR1','USR2','USR3',
        'USR4','USR5','USR6','USR7','USR8',
        'USR9'
        ]
    config["FailureCriteria"] = dict.fromkeys(fcopts, 0)
    TotalStations = blade.ispan.size

    # write the macro file that builds the mesh
    
    fid = open(filename,'wt')
    fid.write('hide_warndlg_keyopt=1\n')
    fid.write('hide_warndlg_areas=1\n')
    fid.write('\n/nerr,500,50000\n')
    fid.write('\n/filename,master\n')
    fid.write('\n/prep7\n' )
    # DEFINE ELEMENT TYPES
    fid.write('\n! DEFINE ELEMENT TYPES =================================\n')
    # structural mass
    fid.write('\n   et,21,mass21,,,0')
    fid.write('\n   r,999,0.0,0.0,0.00001,0.0,0.0,0.0\n')
    # shell281, 8-node structural shell, store data for TOP, BOTTOM,
    # and MID for all layers
    fid.write('\n   et,11,shell281')
    fid.write('\n   keyopt,11,8,2')
    fid.write('\n*if,hide_warndlg_keyopt,eq,1,then')
    fid.write('\n   /UIS, MSGPOP, 3')
    fid.write('\n*endif')
    fid.write('\n   !Set keyopt(2)=1 for improved formulation in R12 & R13')
    fid.write('\n   keyopt,11,2,1')
    fid.write('\n*if,hide_warndlg_keyopt,eq,1,then')
    fid.write('\n   /UIS, MSGPOP, 2')
    fid.write('\n*endif\n')
    #jcb: I thought about checking the ansys version and making conditional
    #statements, but someone could share the shell7 with someone using a
    #different version of ansys
    #     if strncmp(ansys_version,'12',2) || strncmp(ansys_version,'13',2)
    #         # Releases 12 & 13 of ANSYS require keyopt(2)=1 with shell281
    #         #  to make use of the improved formulation
    #         fprintf(fid,'\n   keyopt,11,2,1\n');
    #     end
            # shell181, 4-node structural shell, store data for TOP, BOTTOM, and
    #    MID for all layers
    fid.write('\n   et,12,shell181')
    fid.write('\n   keyopt,12,8,2\n')
    fid.write('\n   keyopt,12,3,2\n')
    #tcl: Write material properties
    #tcl:    This changed dramatically on 2001 November 08
    #tcl:    Now only materials used in the model are written to shell7.src
    #tcl:    Also, material numbers are no longer recorded until write_shell7
    #tcl:    Two new local arrays ansysMPnumber and ansysRnumber are used within write_shell7
    fcfields = [
        'xten','xcmp','yten','ycmp','zten','zcmp',
        'xy','yz','xz','xycp','yzcp','xzcp','xzit',
        'xzic','yzit','yzic','g1g2','etal','etat','alp0'
    ]
    fcvalues = dict.fromkeys(fcfields,0)
    fid.write('\n! WRITE MATERIAL PROPERTIES ============================\n')
    fid.write('\n  ! FAILURE CRITERIA LIMIT TABLE LEGEND:')
    fid.write('\n  ! tb,fcli,<mat>,ntemp,npts,tbopt')
    fid.write('\n  !   (tbopt=1 for stress limits; default ntemp=1,npts=20)')
    fid.write('\n  ! tbdata,1,xten,xcmp,yten,ycmp,zten,zcmp')
    fid.write('\n  ! tbdata,7,xy,yz,xz,xycp,yzcp,xzcp')
    fid.write('\n  ! tbdata,13,xzit,xzic,yzit,yzic')
    fid.write('\n  ! tbdata,17,g1g2,etal,etat,alp0\n')

    for kmp in range(len(blade.materials)):
        mat = blade.materials[kmp]
        if mat.type == 'isotropic':
            fid.write('\n   ! %s' % (mat.name))
            fid.write('\n   mp,ex,%d,%g' % (kmp+1,mat.ex))
            fid.write('\n   mp,dens,%d,%g' % (kmp+1,mat.density))
            fid.write('\n   mp,nuxy,%d,%g' % (kmp+1,mat.prxy))
            xten = mat.uts
        elif 'orthotropic' == mat.type:
            fid.write('\n   ! %s' % (mat.name))
            fid.write('\n   mp,ex,%d,%g' % (kmp+1,mat.ex))
            fid.write('\n   mp,ey,%d,%g' % (kmp+1,mat.ey))
            fid.write('\n   mp,ez,%d,%g' % (kmp+1,mat.ez))
            fid.write('\n   mp,prxy,%d,%g' % (kmp+1,mat.prxy))
            fid.write('\n   mp,pryz,%d,%g' % (kmp+1,mat.pryz))
            fid.write('\n   mp,prxz,%d,%g' % (kmp+1,mat.prxz))
            fid.write('\n   mp,gxy,%d,%g' % (kmp+1,mat.gxy))
            fid.write('\n   mp,gyz,%d,%g' % (kmp+1,mat.gyz))
            fid.write('\n   mp,gxz,%d,%g' % (kmp+1,mat.gxz))
            fid.write('\n   mp,dens,%d,%g' % (kmp+1,mat.density))
        else:
            raise Exception('Unknown material type in database')
        if mat.type in ['isotropic','orthotropic']:
            # Note that entering a blank or a zero for XYCP,YZCP, or XZCP
            # triggers the default value of -1.0. To specify an effective zero,
            # use a small, nonzero value (such as 1E-14).
            #                 if isequal(0,mat.xycp), mat.xycp=1e-14; end
            #                 if isequal(0,mat.yzcp), mat.yzcp=1e-14; end
            #                 if isequal(0,mat.xzcp), mat.xzcp=1e-14; end
            # convert degrees to radians
            #mat.alp0 = mat.alp0 * pi/180;
            # read all of the failure criteria values
            #                 for kfc = 1:numel(fcfields)
            #                     fcname = fcfields{kfc}
            #                     fcvalues{kfc} = mat.(fcname)
            #                 end
            #                 if all(cellfun('isempty',fcvalues))
            #                     # do not print anything if all failure criteria
            #                     # properties are empty
            #                 else
            #                     tempArray=zeros(9,1);
            #Tesile Properties
            try:
                mat.uts.shape
                uts = mat.uts
            except AttributeError:
                uts = np.array([mat.uts])
            nStrenghts = uts.shape[0]
            if nStrenghts < 3:
                uts = fullyPopluateStrengthsArray(uts)
            try:
                mat.ucs.shape
                ucs = mat.ucs
            except AttributeError:
                ucs = np.array([mat.ucs])
            nStrenghts = ucs.shape[0]
            if nStrenghts < 3:
                ucs = fullyPopluateStrengthsArray(ucs)
            try:
                mat.uss.shape
                uss = mat.uss
            except AttributeError:
                uss = np.array([mat.uss])
            nStrenghts = uss.shape[0]
            if nStrenghts < 3:
                uss = fullyPopluateStrengthsArray(uss)
            fid.write('\n   tb,fcli,%d,1,20,1' % (kmp+1))
            fid.write('\n   tbdata,1,%g,%g,%g,%g,%g,%g' % (uts[0],ucs[0],uts[1],ucs[1],uts[2],ucs[2]))
            fid.write('\n   tbdata,7,%g,%g,%g,,,' % (uss[0],uss[1],uss[2]))
            fid.write('\n   tbdata,13,%g,%g,%g,%g' % (mat.xzit,mat.xzic,mat.yzit,mat.yzic))
            fid.write('\n   tbdata,17,%g,%g,%g,%g' % (mat.g1g2,mat.etal,mat.etat,mat.alp0))
            #                     for kf = 1:numel(fcvalues)
            #                         if ~isempty(fcvalues{kf})
            #                             fprintf(fid,'\n   tbdata,#d,#g',kf,fcvalues{kf});
            #                         end
            #                     end
            #end
        else:
            raise Exception('Unknown material type in database')
        fid.write('\n')
    fid.write('\n! WRITE THE COMPOSITE LAYUPS =================================\n')
    rCounter = 1
    if config["elementType"] in ['281','181']:
        ### Outer AeroShell
        nElements = len(meshData['sections'])
        for nelem in range(nElements):
            section = meshData['sections'][nelem]
            # Break from loop once shearweb elements are reached
            stackName = section['elementSet']
            strList = stackName.split('_')
            if strList[2] == "SW":
                swstart = nelem
                break
            
            layup = section['layup']
            secID = nelem+1
            fid.write('\n   ! %s' % (section['elementSet']))
            fid.write('\n   sectype,%d,shell' % (secID))
            for layer in layup:
                matid = layer[0]
                thickness = layer[1] / 1000
                angle = layer[2]
                fid.write('\n      secdata,%g,%d,%g,,' % (thickness,matid+1,angle))
            fid.write('\n   secoffset,bot\n')

        ### Web(s)
        for nelem in range(swstart, nElements):
            section = meshData['sections'][nelem]
            secID = nelem+1
            layup = section['layup']
            if layup:
                # secID = (webSectionIDstart + nstat+1) + (nweb * 10 ** orderOfMagnitude)
                fid.write('\n   ! %s' % (section['elementSet']))
                fid.write('\n   sectype,%d,shell' % (secID))
                for layer in layup:
                    matid = layer[0]
                    thickness = layer[1] / 1000
                    angle = layer[2]
                    fid.write('\n      secdata,%g,%d,%g,,' % (thickness,matid+1,angle))
                fid.write('\n   secoffset,mid\n')
    else:
        # logging.warning('Element System %s not yet available' % config["elementType"],'write_shell7 error')
        raise Exception('Element System %s not yet available',config["elementType"])
    # [~,jobtitle,~] = fileparts(blade.job_name);
    fid.write('\n/title,%s' % (config["dbname"]))
    fid.write('\nZrCount=%d\n' % (rCounter))
    #tcl: DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES
    #tcl:    THE LINES ARE PRODUCED WITH THREE DIFFERENT SPLINING MACROS
    fid.write('\n! DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES\n')
    # Create a coordinate system roughly in the fiber direction (+X down blade, +Z up toward LP side)
    # --> beginning with global csys, rotate -90 about local Z, then -90 about local Y
    fid.write('\nlocal,1000,CART,0,0,0, -90,0,-90\n')
    twistFlag = 1
    if blade.rotorspin == 1:
        twistFlag = - 1
    for kStation in range(TotalStations):
        # use the generating line to translate and rotate the coordinates
        # presweep ===========================================================
        if np.all(blade.isweep == 0):
            table = np.array([0,0,np.nan]).reshape((1,-1))
        else:
            N = len(blade.ispan)
            table = np.array([blade.ispan,blade.isweep,np.full([N,1],np.nan)])
        blade_struct = {}
        blade_struct["PresweepRef"] = {}
        blade_struct["PresweepRef"]["method"] = 'shear'
        blade_struct["PresweepRef"]["table"] = table
        blade_struct["PresweepRef"]["pptype"] = 'spline'
        # precurve ===========================================================
        if np.all(blade.iprebend == 0):
            table = np.array([0,0,np.nan])
        else:
            N = len(blade.ispan)
            table = np.vstack([blade.ispan,blade.iprebend,np.full((N,),np.nan)])
        blade_struct["PrecurveRef"] = {}
        blade_struct["PrecurveRef"]["method"] = 'shear'
        blade_struct["PrecurveRef"]["table"] = table
        blade_struct["PrecurveRef"]["pptype"] = 'spline'
        blade_struct = calcGenLinePP(blade_struct)
        if isinstance(blade_struct["PresweepRef"]["pp"],int):
            presweep_slope = 0
        else:
            presweep_slope = blade_struct["PresweepRef"]["pp"].__call__(blade.ispan[kStation], nu=1)
        if isinstance(blade_struct["PrecurveRef"]["pp"],int):
            precurve_slope = 0
        else:
            precurve_slope = blade_struct["PrecurveRef"]["pp"].__call__(blade.ispan[kStation], nu=1)
        presweepDeg = 180 / np.pi * np.arctan(presweep_slope * twistFlag)
        precurveDeg = 180 / np.pi * np.arctan(- precurve_slope)
        presweep_rot,precurve_rot = (0,0)
        if blade_struct["PresweepRef"]["method"]=='normal':
            presweep_rot = np.arctan(presweep_slope * twistFlag)
        if blade_struct["PrecurveRef"]["method"]=='normal':
            precurve_rot = np.arctan(- precurve_slope)
        transX = twistFlag * blade_struct["PresweepRef"]["pp"].__call__(blade.ispan[kStation])
        transY = blade_struct["PrecurveRef"]["pp"].__call__(blade.ispan[kStation])
        # ensure we are in csys0 and no keypoints are selected
        fid.write('\ncsys,0')
        # Create a coordinate system to be used later for aligning the fiber direction.
        # First, load the csys defined earlier (+X down blade, +Z up toward LP side)
        fid.write('\n   csys,1000')
        # Next, translate & rotate relative to this active csys (use CLOCAL, not LOCAL)
        # translation: global X,Y,Z => local y,z,x
        # rotation: presweep is local z rotation & precurve is local y rotation
        fid.write('\n   clocal,%d,CART,%g,%g,%g, %g,%g,%g\n' % \
            ((1000 + kStation+1), blade.ispan[kStation],transX,transY,presweepDeg,0,precurveDeg))
        # Create coordinate system at the tip
        fid.write('\nlocal,12,CART,%g,%g,%g, %g,%g,%g\n' % \
            (transX,transY,blade.ispan[kStation],0,precurve_rot * 180 / np.pi,presweep_rot * 180 / np.pi))
    fid.write('\n   csys,0')
    fid.write('\nksel,all\n')
    fid.write('\nallsel\n')
    nodes = meshData["nodes"]
    elements = meshData["elements"]
    nnodes = nodes.shape[0]
    nelements = elements.shape[0]
    fid.write('\n! DEFINE NODES =======================================\n')
    for iNode in range(nnodes):
        fid.write('n, %i, %f, %f, %f\n' % \
            (iNode+1,nodes[iNode,0],nodes[iNode,1],nodes[iNode,2]))
    #Set the element Type
    if '281' == config["elementType"]:
        fid.write('type, 11\n')
    else:
        if '181' == config["elementType"]:
            fid.write('type, 12\n')
        else:
            # errordlg(sprintf('Element System %s not yet available',config["elementType"]),'write_shell7 error')
            raise Exception('Element System %s not yet available',config["elementType"])
    dup = []
    fid.write('\n! DEFINE ELEMENTS =======================================\n')
    for iElement in range(nelements):
        if np.unique(elements[iElement,:]).size == 4:
            elem_data = (
                elements[iElement,0]+1,
                elements[iElement,1]+1,
                elements[iElement,2]+1,
                elements[iElement,3]+1,
                iElement+1
            )
            fid.write('e, %i, %i, %i, %i  !Element %i \n' % elem_data)
        else:
            dup.append(iElement)
        
    fid.write('\n! ASSIGN SECTIONS TO OUTER SHELL ELEMENTS =======================================\n')
    nElements = len(meshData['sections'])
    for nelem in range(nElements):
        stackName = meshData['sections'][nelem]['elementSet']
        strList = stackName.split('_')
        if strList[2] == "SW":
            swstart = nelem
            break
        nstat = int(strList[0])
        secID = nelem+1
        csID = 1000 + nstat+1
        elementList = meshData["sets"]["element"][nelem]["labels"]
        for iEl in range(len(elementList)):
            fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl]+1,secID))
            fid.write('   emodif,%i,esys,%i\n' % (elementList[iEl]+1,csID))

    fid.write('\n! ASSIGN SECTIONS TO SHEARWEB(S) SHELL ELEMENTS =======================================\n')
    for nelem in range(swstart, nElements):
        secID = nelem+1
        section = meshData['sections'][nelem]
        layup = section['layup']
        if layup:
            csID = 1000 + nstat
            elementList = meshData["sets"]["element"][nelem]["labels"]
            for iEl in range(len(elementList)):
                fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
                # fprintf(fid,'   emodif,#i,esys,#i\n',elementList(iEl),csID);
    
    fid.write('\n   ENSYM,,,,1,%i' % (nelements))
    #jcb: are these 2 lines necessary now that we have local coordinate
    #  systems to deal with presweep and precurve?
    fid.write('\n   local,11,CART,0,0,0,90,0,-90')
    fid.write('\n   esys,11')
    #LocationZ_lastStation = data.station(TotalStations).LocationZ;
    #fprintf(fid,'\n   local,12,cart,0,0,#f,0,0,0',LocationZ_lastStation);
    fid.write('\n   csys,12')
    fid.write('\n   nsel,none')
    fid.write('\n   n,,0.0,0.0,0.0')
    fid.write('\n   *get,z_master_node_number,node,,num,max')
    fid.write('\n   type,21')
    fid.write('\n   real,999')
    fid.write('\n   e,z_master_node_number')
    fid.write('\n   nsel,all')
    fid.write('\n   csys,0')
    fid.write('\n   allsel\n')
    fid.write('\n   nsll,s,1')
    fid.write('\n   nsel,a,node,,z_master_node_number')
    if config["elementType"] in ['91','99','281','181']:
        fid.write('\n   cerig,z_master_node_number,all,RXYZ\n')
    else:
        if config["elementType"] == '191':
            fid.write('\n   cerig,z_master_node_number,all,uxyz\n')
    if config["BoundaryCondition"]=='cantilevered':
        #jcb: FIXME - nsel could break with swept/bent blades
        fid.write('\n   nsel,s,loc,z,0')
        fid.write('\n   d,all,all')
        fid.write('\n   nsel,all\n')
    fid.write('\nallsel')
    fid.write('\n!   nummrg,all')
    fid.write('\n!   numcmp,node')
    fid.write('\ncsys,0\n')
    ### Material Properties ###
    fid.write('mpwrite,Materials,txt,,\n')
    #if ~all(cellfun('isempty',fcvalues))
    fid.write('/output,Strengths,txt,,\n')
    fid.write('TBLIST, ,ALL\n')
    fid.write('/output\n')
    #end
    # enter POST1 for postprocessing configuration commands
    fid.write('\nfinish')
    fid.write('\n/post1\n')
    fid.write('\nfctyp,dele,all   ! remove all material failure-criteria postprocessing\n')
    for kfc in config["FailureCriteria"].keys():
        if config["FailureCriteria"][kfc]:
            fid.write('fctyp,add,%s\n' % (config["FailureCriteria"][kfc]))
    fid.write('\nfinish\n')
    ### Material Properties ###
    fid.write('mpwrite,Materials,txt,,\n')
    # commenting below since fcvalues are not touched after initialization -kb
    # if not np.all(cellfun('isempty',fcvalues)) :
    #     fid.write('/output,Strengths,txt,,\n')
    #     fid.write('TBLIST, ,ALL\n')
    #     fid.write('/output\n')
    ### Section Properties ###
    fid.write('/output, Sections,txt\n')
    fid.write('SLIST,,,,FULL\n')
    #     fprintf(fid,'SLIST,\n');
    fid.write('/output\n')
    ### Element Properties ###
    fid.write('/output, Elements,txt\n')
    fid.write('elist,all,,,0,0 \n')
    fid.write('/output\n')
    # save database file
    fid.write('\nfinish')
    fid.write('\nsave')
    
    fid.close()
    print('The following file has been written %s\n',filename)
    
    return
    
    
def fullyPopluateStrengthsArray(strengthArray):
    nStrenghts = strengthArray.shape[0]

    if nStrenghts < 3:
        for i in range(3 - nStrenghts):
            strengthArray = np.concatenate([strengthArray,[strengthArray[i]]])
    
    return strengthArray


def write_ansys_loads(nodeData, loads, forcefilename, analysisConfig):
    """
    """
    maptype = 'map3D_fxM0'
    if ('FollowerForces' in analysisConfig["analysisFlags"]) and \
        not len(analysisConfig["analysisFlags"].FollowerForces)==0 \
        and analysisConfig["analysisFlags"].FollowerForces != 0 and \
        ('StaticNonlinear' in analysisConfig["analysisFlags"]) and not \
        len(analysisConfig["analysisFlags"].StaticNonlinear)==0  and \
        analysisConfig["analysisFlags"]["StaticNonlinear"] != 0:
        forcemap,forcesums = beamForceToAnsysShellFollower(nodeData,loads, maptype=maptype)
    else:
        forcemap,forcesums = beamForceToAnsysShell(nodeData,loads, maptype=maptype)

    writeforcefile(forcefilename+'.src',forcemap,forcesums,maptype)
    print('Forces mapped to ANSYS model')
    return