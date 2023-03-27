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
                fid.write('\n   mp,ex,%d,%g' % (kmp,mat.ex))
                fid.write('\n   mp,dens,%d,%g' % (kmp,mat.density))
                fid.write('\n   mp,nuxy,%d,%g' % (kmp,mat.prxy))
                xten = mat.uts
            elif 'orthotropic' == mat.type:
                fid.write('\n   ! %s' % (mat.name))
                fid.write('\n   mp,ex,%d,%g' % (kmp,mat.ex))
                fid.write('\n   mp,ey,%d,%g' % (kmp,mat.ey))
                fid.write('\n   mp,ez,%d,%g' % (kmp,mat.ez))
                fid.write('\n   mp,prxy,%d,%g' % (kmp,mat.prxy))
                fid.write('\n   mp,pryz,%d,%g' % (kmp,mat.pryz))
                fid.write('\n   mp,prxz,%d,%g' % (kmp,mat.prxz))
                fid.write('\n   mp,gxy,%d,%g' % (kmp,mat.gxy))
                fid.write('\n   mp,gyz,%d,%g' % (kmp,mat.gyz))
                fid.write('\n   mp,gxz,%d,%g' % (kmp,mat.gxz))
                fid.write('\n   mp,dens,%d,%g' % (kmp,mat.density))
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
                fid.write('\n   tb,fcli,%d,1,20,1' % (kmp))
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
            maxSectionNumber = int(str(nStations)+str(nStationLayups))
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
                        fid.write('\n      secdata,%g,%d,%g,,' % (layerThickness,matID,layerLayupAngle))
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
                        secID = webSectionIDstart + iStation + (iWeb - 1) * 10 ** orderOfMagnitude
                        fid.write('\n   ! %s' % (currentStack.name))
                        fid.write('\n   sectype,%d,shell' % (secID))
                        for iLayer in range(len(currentStack.plygroups)):
                            currentLayer = currentStack.plygroups[iLayer]
                            layerThickness = currentLayer.nPlies * currentLayer.thickness / 1000
                            layerLayupAngle = currentLayer.angle
                            matID = currentLayer.materialid
                            fid.write('\n      secdata,%g,%d,%g,,' % (layerThickness,matID,layerLayupAngle))
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
                ((1000 + kStation), blade.ispan[kStation],transX,transY,presweepDeg,0,precurveDeg))
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
                (iNode,nodes[iNode,0],nodes[iNode,1],nodes[iNode,2]))
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
                    iElement
                )
                fid.write('e, %i, %i, %i, %i  !Element %i \n' % elem_data)
            else:
                dup = np.array([[dup],[iElement]])

        # Separate elements into outer shell and shear web
        # outershell_els = []
        # shearweb_els = []
        elements = meshData['sets']['element']
        sections = meshData['sections']
        # for i in range(len(elements)):
        #     if "SW" in elements[i]["name"]:
        #         shearweb_els.append((elements[i], sections[i]))
        #     else:
        #         outershell_els.append((elements[i], sections[i]))
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
        
#        fid.write('\n! ASSIGN SECTIONS TO OUTER SHELL ELEMENTS =======================================\n')
#        for iStation in range(nStations):
#            for iPerimeter in range(nStationLayups):
#                secID = int(str(iStation) + str(iPerimeter))
#                csID = 1000 + iStation
#                elementList = meshData["outerShellElSets"][iPerimeter,iStation].elementList
#                for iEl in range(len(elementList)):
#                    fid.write('   emodif,%i,secnum,%i\n' % (elementList(iEl),secID))
#                    fid.write('   emodif,%i,esys,%i\n' % (elementList(iEl),csID))
#        fid.write('\n! ASSIGN SECTIONS TO SHEARWEB(S) SHELL ELEMENTS =======================================\n')
#        nWebs = len(blade.swstacks)
#        for iWeb in range(nWebs):
#            __,nStations = blade.swstacks[iWeb].shape
#            for iStation in range(nStations):
#                currentStack = blade.swstacks[iWeb][iStation]
#                if not len(currentStack.plygroups)==0 :
#                    secID = webSectionIDstart + iStation + (iWeb - 1) * 10 ** orderOfMagnitude
#                    csID = 1000 + iStation
#                    elementList = meshData["shearWebElSets"][iWeb][iStation].elementList
#                    for iEl in range(len(elementList)):
#                        fid.write('   emodif,%i,secnum,%i\n' % (elementList[iEl],secID))
                        #fprintf(fid,'   emodif,#i,esys,#i\n',elementList(iEl),csID);
        #     #tcl: reverse area normals if clockwise blade
        #     #tcl:    shear web areas are reversed as well - not necessary, just easier
        #     if blade.rotorspin == 1 # clockwise rotation
        #         fprintf(fid,'\n   areverse,all');
        #     end
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
    

def writeANSYSinputFile(fid, mat, ansysSecNumber, coreMatName): 
    #####Find the face sheet####
    cellMat = np.array([])
    for i in range(len(mat.layer)):
        cellMat = np.array([[cellMat],[np.array([mat.layer[i].layerName])]])
    
    kbalsa = np.where((str(coreMatName) == 1) and (str(cellMat) == 1))
    iLayer = range(kbalsa - 1)
    
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
    fname = 'section-' + str(ansysSecNumber) + '-faceAvgStresses'
    
    if os.path.isfile(fname + '.txt'):
        os.delete(fname + '.txt')
    
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


import numpy as np
    
def writeAnsysDeflections(blade, config, iLoad, fid, deflectionFilename): 
    #Outer AeroShell
    nStationLayups,nStations = blade.stacks.shape
    maxSectionNumber = int(str(nStations) + str(nStationLayups))
    
    #The following two lines help make unique IDs for web sections
#based on the highes section already defined for aeroshell
    orderOfMagnitude = int(np.floor(np.log10(maxSectionNumber)))
    webSectionIDstart = np.ceil(maxSectionNumber / 10 ** orderOfMagnitude) * 10 ** orderOfMagnitude
    fid.write('/POST1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('RSYS,0\n' % ())
    
    fid.write('seltol,0.05\n' % ())
    for i in range(blade.ispan.size):
        fid.write('*CFOPEN, %s,out\n' % (deflectionFilename + '-' + str(i)))
        fid.write('ESEL,S,SEC,,1,%i   \n' % (webSectionIDstart))
        #fprintf(fid,'ESEL,S,SEC,,1,999   \n');    #Selects aero shell only
        fid.write('nsle,S,   \n' % ())
        fid.write('nsel,r,loc,z,%f  \n' % (blade.ispan(i)))
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
    
## Panel Stresses Analysis Script
    
def writeAnsysGetFaceStresses(blade, fid, coreMatName): 
    isoorthoInModel,compsInModel,SkinAreas,app = getMatrialLayerInfoWithOutGUI(blade)
    #fid=fopen('getFaceStresses.mac','w+');
    TotalStations = blade.ispan.size
    for kStation in range(TotalStations - 1):
        #kPanel=find(~cellfun('isempty',strfind([SkinAreas(kStation).Material],'PANEL'))); #Array that stores the kArea index that contains 'PANEL' in the name
        #for i=1:numel(kPanel)
        for kArea in range(SkinAreas[kStation].startIB.size):
            #See if the section contatins Balsa/core material name (i.e find
            #the sandwhich panels)
            n = str(SkinAreas[kStation].Material[kArea]) == str(app.matlist)
            mat = app.matdb[n]
            if coreMatName in [mat.layer.layerName]:
                ansysSecNumber = np.where(str(SkinAreas[kStation].Material[kArea]) == str(compsInModel) == 1)
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
        if coreMatName in [mat.layer.layerName]:
            ansysSecNumber = np.where(str([app.shearweb[kShearweb].Material]) == str(compsInModel) == 1)
            ansysSecNumber = ansysSecNumber + 1000
            writeANSYSinputFile(fid,mat,ansysSecNumber,coreMatName)
    
    fid.write('FINISH\n' % ())
    fid.write('allsel\n' % ())
    return app,SkinAreas,compsInModel

    
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
    fid.write('nz=nint(%f/elsize) !Integer number of points to output resultant loads\n' % (blade.ispan[-1]))
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
    

def writeAnsysLocalFields(blade, config, iLoad, fid): 
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
    if config.analysisFlags.failure.upper() not in ['PUCK','LARC03','LARC04']:
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
    

def writeAnsysLinearBuckling(blade, config, iLoad, fid, bucklingFilename): 
    fid.write('! BEGIN BUCKLE MACRO\n' % ())
    fid.write('allsel\n' % ())
    fid.write('/solu\n' % ())
    fid.write('irlf,-1\n' % ())
    fid.write('pstres,on\n' % ())
    fid.write('antype,buckle\n' % ())
    fid.write('bucopt,lanb,' + str(config.analysisFlags.globalBuckling) + ',,,RANGE\n' % ())
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
    

def writeAnsysNonLinearBuckling(ansysFilename = None,ansys_path = None,ansys_product = None,config = None,ii = None,jj = None,ncpus = None,iLoad = None): 
    warnings.warn('output designvar. Currently does not work for nonlinear cases')
    script_name = 'commands3-'+str(ii)+'.mac'
    script_out = 'output3-'+str(ii)+'-'+str(jj)+'.txt'
    fid = open(script_name,'w+')
    fid.write('!************   MODE-%i   ************\n') % (ii)
    fid.write('/FILNAME,'%s',1\n' % (strcat(ansysFilename,'-Load',int2str(iLoad))))
    
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
    fid.write('UPGEOM,zImperfectionSF,1,%i,'%s','rst'\n' % (ii,strcat(ansysFilename,'-Load',int2str(iLoad))))
    
    fid.write('FINISH\n' % ())
    
    filename = ansysFilename+'-'+str(ii)+'-',str(jj)
    
    fid.write('/FILNAME,'%s',1\n' % (filename))
    
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
    fid.write('SAVE','',filename,'','db','',pwd.replace('\','\\'),''\n') % ())
    fid.write('/EXIT,NOSAVE\n' % ())
    fid.close()
    ###### RUN ANSYS#######
    ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %s',ansys_path,ansys_product,script_name,script_out,int2str(ncpus))
    # KMP_STACKSIZE=2048k has been specifed. 2048k may not be enough for other
# simulations. EC
    
    system(ansys_call)
    print('%s: Nonlinear Mode-%s Analysis Finished\n' % (datestr(now),int2str(ii)))
    data = readANSYSoutputs(strcat(filename,'.mntr'),11)
    a = data.shape
    nonlinearLoadFactors = data(a(1),7) * loadScaleFactor
    
    return nonlinearLoadFactors


def writeAnsysFagerberWrinkling(app = None,SkinAreas = None,compsInModel = None,coreMatName = None): 
    #limitingElementData - [ansysSecNumber elno lf phicr]
    TotalStations = np.asarray(app.station).size
    TotalShearwebs = np.asarray(app.shearweb).size
    #################   Main loop #1: loop around aero shell.   #################
    LF = []
    for kStation in np.arange(1,TotalStations - 1+1).reshape(-1):
        for kArea in np.arange(1,np.asarray(SkinAreas(kStation).startIB).size+1).reshape(-1):
            #See if the section contatins Balsa/core material name (i.e find
#the sandwhich panels)
            n = str(SkinAreas(kStation).Material[kArea]) == str(app.matlist)
            mat = app.matdb(n)
            if contains(np.array([mat.layer.layerName]),coreMatName):
                ansysSecNumber = find(str(SkinAreas(kStation).Material(kArea)) == str(compsInModel) == 1)
                file = strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses.txt')
                avgFaceStress = txt2mat(file)
                os.delete(file)
                LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName)
    
    #################   Main loop #2: loop along web.   #################
    for kShearweb in np.arange(1,TotalShearwebs+1).reshape(-1):
        n = str(app.shearweb(kShearweb).Material) == str(app.matlist)
        mat = app.matdb(n)
        if contains(np.array([mat.layer.layerName]),coreMatName):
            ansysSecNumber = find(str(np.array([app.shearweb(kShearweb).Material])) == str(compsInModel) == 1)
            file = strcat('section-',int2str(ansysSecNumber + 1000),'-faceAvgStresses.txt')
            avgFaceStress = txt2mat(file)
            os.delete(file)
            LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber + 1000,avgFaceStress,app,mat,coreMatName)
    
    minLF,index = np.amin(LF(:,3))
    limitingElementData = LF(index,:)
    print('\n\n The minimum wrinkling LF is: %f, wrinkle angle: %.2f°' % (minLF,LF(index,4)))
    print('\n and occurs in section number %i, element number %i\n, ' % (LF(index,1),LF(index,2)))
    # [maxLF,index]=max(LF(:,3));
# fprintf('\n\n The maximum LF is: #f, wrinkle  angle: #.2f°' ,maxLF, LF(index,4))
# fprintf('\n and occurs in section number #i, element number #i, ',LF(index,1),LF(index,2))
    
    return limitingElementData


def postprocessANSYSfatigue(blade = None,meshData = None,wt = None,rccdata = None,IEC = None,loadsTable = None,config = None): 
    if np.any(contains(config.analysisFlags.fatigue.lower(),'all')):
        nSegments = 1
    else:
        nSegments = np.asarray(config.analysisFlags.fatigue).size
    
    # Order of the segment names in segmentNamesReference
# is very important. config.analysisFlags.fatigue can
# be any order
    segmentNamesReference = np.array(['HP_TE_FLAT','HP_TE_ReINF','HP_TE_PANEL','HP_SPAR','HP_LE_PANEL','HP_LE','LP_LE','LP_LE_PANEL','LP_SPAR','LP_TE_PANEL','LP_TE_REINF','LP_TE_FLAT'])
    nsegmentNamesReference = np.asarray(segmentNamesReference).size
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
        nGage = np.asarray(loadsTableTheta.input.rGage).size
        gageNumber = np.transpose((np.arange(1,nGage+1)))
        criticalElement,fatigueDamage,criticalLayerNo,criticalMatNo = deal(np.zeros((nGage,1)))
        criticalMat = cell(nGage,1)
        rGage = loadsTableTheta.input.rGage
        MrTheta = loadsTableTheta.input.Mrb
        MrThetaPlus90 = loadsTableThetaPlus90.input.Mrb
        plotFatigue = []
        fileNameTheta = np.array(['plateStrains-all-',int2str(kTheta),'.txt'])
        fileNameThetaPlus90 = np.array(['plateStrains-all-',int2str(kTheta + nSpace),'.txt'])
        print(fileNameTheta)
        print(fileNameThetaPlus90)
        #Used for reading element stresses
        pat = 'ELEM\s*ZCENT\s*EPS11\s*EPS22\s*EPS12\s*KAPA11\s*KAPA22\s*KAPA12\s*GAMMA13\s*GAMMA23'
        NCOLS = 10
        plateStrainsTheta = readANSYSElementTable(fileNameTheta,pat,NCOLS)
        plateStrainsThetaPlus90 = readANSYSElementTable(fileNameThetaPlus90,pat,NCOLS)
        for i in np.arange(1,nSegments+1).reshape(-1):
            iSegment = find(strcmpi(segmentNamesReference,config.analysisFlags.fatigue(i)) == 1)
            if np.any(contains(config.analysisFlags.fatigue.lower(),'all')):
                title = 'All segments'
            else:
                if not strcmpi(config.analysisFlags.fatigue(i),'webs') :
                    title = config.analysisFlags.fatigue(i)
                    __,nSpanRegions = meshData.outerShellElSets.shape
                    elementList = []
                    for iSpan in np.arange(1,nSpanRegions+1).reshape(-1):
                        elementList = np.array([elementList,meshData.outerShellElSets(iSegment,iSpan).elementList])
                else:
                    title = 'Webs'
                    __,nWebs = meshData.shearWebElSets.shape
                    elementList = []
                    for iWeb in np.arange(1,nWebs+1).reshape(-1):
                        __,nSpanRegions = meshData.shearWebElSets[iWeb].shape
                        for iSpan in np.arange(1,nSpanRegions+1).reshape(-1):
                            elementList = np.array([elementList,meshData.shearWebElSets[iWeb](iSpan).elementList])
                plateStrainsThetaSet = plateStrainsTheta(elementList,:)
                plateStrainsThetaPlus90Set = plateStrainsThetaPlus90(elementList,:)
            for chSpan in np.arange(1,nGage+1).reshape(-1):
                direction = int2str(theta)
                Ltheta = getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction)
                if theta + 90 < 180:
                    direction = int2str(theta + 90)
                else:
                    direction = int2str(theta - 90)
                LthetaPlus90 = getMomentMarkov(rccdata,wt,Yr,simtime,markovSize,chSpan,direction)
                Mtheta = interp1(rGage,MrTheta,rGage(chSpan))
                MthetaPlus90 = interp1(rGage,MrThetaPlus90,rGage(chSpan))
                zwidth = 0.75
                # at a blade gage location.
                z1 = rGage(chSpan) - zwidth / 2
                z2 = rGage(chSpan) + zwidth / 2
                binnedElements = intersect(find(plateStrainsThetaSet(:,2) < z2),find(plateStrainsThetaSet(:,2) > z1))
                fdData,plotFatigueChSpan = calcFatigue(blade,meshData,IEC,Ltheta,LthetaPlus90,Mtheta,MthetaPlus90,binnedElements,plateStrainsThetaSet,plateStrainsThetaPlus90Set,iSegment)
                plotFatigue = np.array([[plotFatigue],[plotFatigueChSpan]])
                criticalElement[chSpan] = fdData(1)
                fatigueDamage[chSpan] = fdData(2)
                criticalLayerNo[chSpan] = fdData(5)
                criticalMatNo[chSpan] = fdData(8)
                criticalMat[chSpan] = blade.materials(fdData(8)).name
            #         plotFatigueFileName=['plotFatigue-' int2str(kTheta)];
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
    n = len(plotFatigue(:,1))
    fid = open(np.array([fname,'.txt']),'w+')
    fid.write('Element fatigueDamage\n' % ())
    plotFatigue = sortrows(plotFatigue,1)
    for i in np.arange(np.arange(1,len(plotFatigue(,,1))+1)):
        fid.write('%8i  %6.5E\n' % (plotFatigue(i,1),plotFatigue(i,2)))
    
    fid.close()
    #Write plot commands
    fid = open(np.array([fname,'.mac']),'w+')
    fid.write('/post1\n' % ())
    fid.write('set,last\n' % ())
    fid.write('plnsol,u,sum\n' % ())
    fid.write('etab,test,u,X\n' % ())
    fid.write('*get,max_e,elem,0,count\n' % ())
    fid.write('*dim,d_res,array,max_e,3\n' % ())
    fid.write('!Column 1 = Element Number\n' % ())
    fid.write('!Column 2 = Where I'm putting result data\n' % ())
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
    
    return designVar

    
def writeAnsysFagerberWrinkling(app = None,SkinAreas = None,compsInModel = None,coreMatName = None): 
    #limitingElementData - [ansysSecNumber elno lf phicr]
    TotalStations = np.asarray(app.station).size
    TotalShearwebs = np.asarray(app.shearweb).size
    #################   Main loop #1: loop around aero shell.   #################
    LF = []
    for kStation in np.arange(1,TotalStations - 1+1).reshape(-1):
        for kArea in np.arange(1,np.asarray(SkinAreas(kStation).startIB).size+1).reshape(-1):
            #See if the section contatins Balsa/core material name (i.e find
            #the sandwhich panels)
            n = str(SkinAreas(kStation).Material[kArea]) == str(app.matlist)
            mat = app.matdb(n)
            if contains(np.array([mat.layer.layerName]),coreMatName):
                ansysSecNumber = find(str(SkinAreas(kStation).Material(kArea)) == str(compsInModel) == 1)
                file = strcat('section-',int2str(ansysSecNumber),'-faceAvgStresses.txt')
                avgFaceStress = txt2mat(file)
                os.delete(file)
                LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber,avgFaceStress,app,mat,coreMatName)
    
    #################   Main loop #2: loop along web.   #################
    for kShearweb in np.arange(1,TotalShearwebs+1).reshape(-1):
        n = str(app.shearweb(kShearweb).Material) == str(app.matlist)
        mat = app.matdb(n)
        if contains(np.array([mat.layer.layerName]),coreMatName):
            ansysSecNumber = find(str(np.array([app.shearweb(kShearweb).Material])) == str(compsInModel) == 1)
            file = strcat('section-',int2str(ansysSecNumber + 1000),'-faceAvgStresses.txt')
            avgFaceStress = txt2mat(file)
            os.delete(file)
            LF = getLoadFactorsForElementsWithSameSection(LF,ansysSecNumber + 1000,avgFaceStress,app,mat,coreMatName)
    
    minLF,index = np.amin(LF(:,3))
    limitingElementData = LF(index,:)
    print('\n\n The minimum wrinkling LF is: %f, wrinkle angle: %.2f°' % (minLF,LF(index,4)))
    print('\n and occurs in section number %i, element number %i\n, ' % (LF(index,1),LF(index,2)))
    # [maxLF,index]=max(LF(:,3));
    # fprintf('\n\n The maximum LF is: #f, wrinkle  angle: #.2f°' ,maxLF, LF(index,4))
    # fprintf('\n and occurs in section number #i, element number #i, ',LF(index,1),LF(index,2))
    
    return limitingElementData