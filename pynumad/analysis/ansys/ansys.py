import numpy as np
import logging
import warnings
import os
from pynumad.utils.interpolation import calcGenLinePP

def writeANSYSshellModel(blade, filename, meshData, config, includeAdhesive):
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
    try:
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
            #Outer AeroShell
            nStationLayups,nStations = blade.stacks.shape
            maxSectionNumber = int(str(nStations+1)+str(nStationLayups+1))
            for iStation in range(nStations):
                for iPerimeter in range(nStationLayups):
                    secID = int(str(iStation+1) + str(iPerimeter+1))
                    currentStack = blade.stacks[iPerimeter,iStation]
                    fid.write('\n   ! %s' % (currentStack.name))
                    fid.write('\n   sectype,%d,shell' % (secID))
                    for iLayer in range(len(currentStack.plygroups)):
                        currentLayer = currentStack.plygroups[iLayer]
                        layerThickness = currentLayer.nPlies * currentLayer.thickness / 1000
                        layerLayupAngle = currentLayer.angle
                        matID = currentLayer.materialid
                        fid.write('\n      secdata,%g,%d,%g,,' % (layerThickness,matID+1,layerLayupAngle))
                    fid.write('\n   secoffset,bot\n')
            #Web(s)
            nWebs = len(blade.swstacks)
            #The following two lines help make unique IDs for web sections
            #based on the highes section already defined for aeroshell
            orderOfMagnitude = int(np.floor(np.log10(maxSectionNumber)))
            webSectionIDstart = np.ceil(maxSectionNumber / 10 ** orderOfMagnitude) * 10 ** orderOfMagnitude
            for iWeb in range(nWebs):
                nStations = len(blade.swstacks[iWeb])
                for iStation in range(nStations):
                    currentStack = blade.swstacks[iWeb][iStation]
                    if not len(currentStack.plygroups)==0 :
                        secID = (webSectionIDstart + iStation+1) + (iWeb * 10 ** orderOfMagnitude)
                        fid.write('\n   ! %s' % (currentStack.name))
                        fid.write('\n   sectype,%d,shell' % (secID))
                        for iLayer in range(len(currentStack.plygroups)):
                            currentLayer = currentStack.plygroups[iLayer]
                            layerThickness = currentLayer.nPlies * currentLayer.thickness / 1000
                            layerLayupAngle = currentLayer.angle
                            matID = currentLayer.materialid
                            fid.write('\n      secdata,%g,%d,%g,,' % (layerThickness,matID+1,layerLayupAngle))
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
        #jcb: as of 2011-05-26, the keypoints are transformed directly
        #     rather than with ansys commands
        #     fprintf(fid,'\n! ROTATE SECTIONS ======================================\n');
        #     fprintf(fid,'\ncsys,1');
        #     for kStation = 1:TotalStations
        #         fprintf(fid,'\n   lsel,s,loc,z,#g',data.station(kStation).LocationZ);
        #         fprintf(fid,'\n   lgen,2,all,,,,#g,,,,1\n',twistFlag*data.station(kStation).DegreesTwist);
        #     end
        #     fprintf(fid,'\ncsys,0');
        fid.write('\nallsel\n')
        nodes = meshData["nodes"]
        elements = meshData["elements"]
        nnodes = nodes.shape[0]
        nelements = elements.shape[0]
        fid.write('\n! DEFINE NODES =======================================\n')
        for iNode in range(nnodes):
            fid.write('n, %i, %f, %f, %f\n' % \
                (iNode+1,nodes[iNode,0],nodes[iNode,1],nodes[iNode,2]))
        #Set the elemetn Type
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
                    elements[iElement,0],
                    elements[iElement,1],
                    elements[iElement,2],
                    elements[iElement,3],
                    iElement+1
                )
                fid.write('e, %i, %i, %i, %i  !Element %i \n' % elem_data)
            else:
                dup.append(iElement)

        # Separate elements into outer shell and shear web
        outershell_els = []
        shearweb_els = []
        elements = meshData['sets']['element']
        sections = meshData['sections']
        for i in range(len(elements)):
            if "SW" in elements[i]["name"]:
                shearweb_els.append((elements[i], sections[i]))
            else:
                outershell_els.append((elements[i], sections[i]))
            
        fid.write('\n! ASSIGN SECTIONS TO OUTER SHELL ELEMENTS =======================================\n')
        for iStation in range(nStations):
            for iPerimeter in range(nStationLayups):
                secID = int(str(iStation+1) + str(iPerimeter+1))
                csID = 1000 + iStation
                elementList = meshData["outerShellElSets"][iPerimeter,iStation].elementList
                for iEl in range(len(elementList)):
                    fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
                    fid.write('   emodif,%i,esys,%i\n' % (elementList[iEl],csID))
        
        fid.write('\n! ASSIGN SECTIONS TO SHEARWEB(S) SHELL ELEMENTS =======================================\n')
        nWebs = len(blade.swstacks)
        for iWeb in range(nWebs):
            __,nStations = blade.swstacks[iWeb].shape
            for iStation in range(nStations):
                currentStack = blade.swstacks[iWeb][iStation]
                if not len(currentStack.plygroups)==0 :
                    secID = (webSectionIDstart + iStation + 1) + (iWeb * 10 ** orderOfMagnitude)
                    csID = 1000 + iStation
                    elementList = meshData["shearWebElSets"][iWeb][iStation].elementList
                    for iEl in range(len(elementList)):
                        fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
                        # fprintf(fid,'   emodif,#i,esys,#i\n',elementList(iEl),csID);
       

        #%tcl: reverse area normals if clockwise blade
        #%tcl:    shear web areas are reversed as well - not necessary, just easier
        #if blade.rotorspin == 1 % clockwise rotation
        #    fprintf(fid,'\n   areverse,all');
        #end 
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
        #jcb: TO BE DELETED
        #     fprintf(fid,'\n   nsel,s,loc,z,#.7f,#.7f\n',...
        #         (LocationZ_lastStation - 0.0000001),...
        #         (LocationZ_lastStation + 0.0000001));
        # select tip station lines and then nodes attached to those lines
        #jcb: I think this can be cleaned up by moving these after the
        #     'e,z_master_node_number' command above and changing 's' to 'a'
        #     below ('nsel' command is then unnecessary because z_master_node
        #     is already selected)
        #fprintf(fid,'\n   cmsel,s,tip_station_lines');
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
    finally:
        pass
    
    fid.close()
    msg = print('The following file has been written %s\n',filename)
    if 1:
        print(msg)
    else:
        # msgbox(msg,'Notification')
        pass
    
    return
    
    
def fullyPopluateStrengthsArray(strengthArray):
    nStrenghts = strengthArray.shape[0]

    if nStrenghts < 3:
        for i in range(3 - nStrenghts):
            strengthArray = np.concatenate([strengthArray,[strengthArray[i]]])
    
    return strengthArray


"""
# Attempting to write all elements together
        fid.write('\n! ASSIGN SECTIONS TO ELEMENTS =======================================\n')
        secID = 1
        for el in elements:
            # secID = int(str(iStation) + str(iPerimeter)) #not sure how to get this
            # secID = 0
            # csID = 1000 + iStation #not sure how to get this
            # csID = 1
            # elementList = meshData["outerShellElSets"][iPerimeter,iStation].elementList
            elementList = el['labels']
            for iEl in range(len(elementList)):
                fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
                # fid.write('   emodif,%i,esys,%i\n' % (elementList[iEl],csID))
            secID+=1
        # secID = 1
        # fid.write('\n! ASSIGN SECTIONS TO SHEARWEB(S) SHELL ELEMENTS =======================================\n')
        # for el, section in shearweb_els:
        #     if not len(section['layup'])==0:
        #         # secID = webSectionIDstart + iStation + (iWeb - 1) * 10 ** orderOfMagnitude #not how to get this
        #         secID = 0
        #         # csID = 1000 + iStation # not sure how to get this
        #         elementList = el['labels']
        #         for iEl in range(len(elementList)):
        #             fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
        #         secID+=1

        
        # #old version
        # nWebs = len(blade.swstacks)
        # for iWeb in range(nWebs):
        #     __,nStations = blade.swstacks[iWeb].shape
        #     for iStation in range(nStations):
        #         currentStack = blade.swstacks[iWeb][iStation]
        #         if not len(currentStack.plygroups)==0 :
        #             secID = webSectionIDstart + iStation + (iWeb - 1) * 10 ** orderOfMagnitude
        #             csID = 1000 + iStation
        #             elementList = meshData["shearWebElSets"][iWeb][iStation].elementList
        #             for iEl in range(len(elementList)):
        #                 fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
        
        #     #tcl: reverse area normals if clockwise blade
        #     #tcl:    shear web areas are reversed as well - not necessary, just easier
        #     if blade.rotorspin == 1 # clockwise rotation
        #         fprintf(fid,'\n   areverse,all');
        #     end
"""