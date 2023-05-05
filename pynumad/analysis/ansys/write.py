import numpy as np
import os
import warnings
import subprocess
from datetime import datetime

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
    fid.write('SAVE,'',filename,'','db','',pwd.replace('\','\\'),''\n') % ())
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

