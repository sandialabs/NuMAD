function [nodes,elements,outerShellElSets,shearWebElSets] = writeANSYSshellModel(blade,filename,fea)
%WRITE_SHELL7  Generate the ANSYS input file that creates the blade 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   write_shell7(app,blade,filename)
%     app = app data structure in NuMAD
%     blade = blade data structure in NuMAD
%     filename = ANSYS input file to be written (typically 'shell7.src')
%


global numadPath
if isequal(0,nargout)
    parent_pn = [numadPath '\numadObjects'];
    [success,message,~] = copyfile(fullfile(parent_pn,'FEA','macros','*.mac'),blade.paths.job);
    if ~success
        errordlg(message,'write_shell7: error copying macros');
        return;
    end
end


TotalStations = numel(blade.ispan);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% write the shell7.src file %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fid = fopen('shell7.src','wt');
fid = fopen(filename,'wt');
try
    fprintf(fid,'hide_warndlg_keyopt=1\n');
    fprintf(fid,'hide_warndlg_areas=1\n');
    
    fprintf(fid,'\n/nerr,500,50000\n');
    fprintf(fid,'\n/filename,master\n');
    
    fprintf(fid,'\n/prep7\n');
    
    %DEFINE ELEMENT TYPES
    fprintf(fid,'\n! DEFINE ELEMENT TYPES =================================\n');

    %structural mass
    fprintf(fid,'\n   et,21,mass21,,,0');
    fprintf(fid,'\n   r,999,0.0,0.0,0.00001,0.0,0.0,0.0\n');

    %shell281, 8-node structural shell, store data for TOP, BOTTOM,
    %     and MID for all layers
    fprintf(fid,'\n   et,11,shell281');
    fprintf(fid,'\n   keyopt,11,8,2');
    fprintf(fid,'\n*if,hide_warndlg_keyopt,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 3'); % turn warning dialog off
    fprintf(fid,'\n*endif');
    fprintf(fid,'\n   !Set keyopt(2)=1 for improved formulation in R12 & R13');
    fprintf(fid,'\n   keyopt,11,2,1');
    fprintf(fid,'\n*if,hide_warndlg_keyopt,eq,1,then');
    fprintf(fid,'\n   /UIS, MSGPOP, 2'); % turn warning dialog on
    fprintf(fid,'\n*endif\n');
    %jcb: I thought about checking the ansys version and making conditional
    %statements, but someone could share the shell7 with someone using a
    %different version of ansys
%     if strncmp(ansys_version,'12',2) || strncmp(ansys_version,'13',2)
%         % Releases 12 & 13 of ANSYS require keyopt(2)=1 with shell281
%         %  to make use of the improved formulation
%         fprintf(fid,'\n   keyopt,11,2,1\n');
%     end
  
    % shell181, 4-node structural shell, store data for TOP, BOTTOM, and 
    %    MID for all layers
    fprintf(fid,'\n   et,12,shell181');
    fprintf(fid,'\n   keyopt,12,8,2\n');
    fprintf(fid,'\n   keyopt,12,3,2\n');  % added by brr 10/25/2012

    %tcl: Write material properties
    %tcl:    This changed dramatically on 2001 November 08
    %tcl:    Now only materials used in the model are written to shell7.src
    %tcl:    Also, material numbers are no longer recorded until write_shell7
    %tcl:    Two new local arrays ansysMPnumber and ansysRnumber are used within write_shell7
    
    fcfields = {'xten','xcmp','yten','ycmp','zten','zcmp','xy','yz','xz',...
        'xycp','yzcp','xzcp','xzit','xzic','yzit','yzic','g1g2',...
        'etal','etat','alp0'};
    fcvalues = cell(numel(fcfields),1);
    
    fprintf(fid,'\n! WRITE MATERIAL PROPERTIES ============================\n');
    fprintf(fid,'\n  ! FAILURE CRITERIA LIMIT TABLE LEGEND:');
    fprintf(fid,'\n  ! tb,fcli,<mat>,ntemp,npts,tbopt');
    fprintf(fid,'\n  !   (tbopt=1 for stress limits; default ntemp=1,npts=20)');
    fprintf(fid,'\n  ! tbdata,1,xten,xcmp,yten,ycmp,zten,zcmp');
    fprintf(fid,'\n  ! tbdata,7,xy,yz,xz,xycp,yzcp,xzcp');
    fprintf(fid,'\n  ! tbdata,13,xzit,xzic,yzit,yzic');
    fprintf(fid,'\n  ! tbdata,17,g1g2,etal,etat,alp0\n');
    for kmp = 1:numel(blade.materials)

        mat = blade.materials(kmp);
        switch mat.type
            case 'isotropic'
                fprintf(fid,'\n   ! %s'           ,mat.name);
                fprintf(fid,'\n   mp,ex,%d,%g'    ,kmp,mat.ex);
                fprintf(fid,'\n   mp,dens,%d,%g'  ,kmp,mat.density);
                fprintf(fid,'\n   mp,nuxy,%d,%g',kmp,mat.prxy);
                xten=mat.uts;
                
            case 'orthotropic'
                fprintf(fid,'\n   ! %s'         ,mat.name);
                fprintf(fid,'\n   mp,ex,%d,%g'  ,kmp,mat.ex);
                fprintf(fid,'\n   mp,ey,%d,%g'  ,kmp,mat.ey);
                fprintf(fid,'\n   mp,ez,%d,%g'  ,kmp,mat.ez);
                fprintf(fid,'\n   mp,prxy,%d,%g',kmp,mat.prxy);
                fprintf(fid,'\n   mp,pryz,%d,%g',kmp,mat.pryz);
                fprintf(fid,'\n   mp,prxz,%d,%g',kmp,mat.prxz);
                fprintf(fid,'\n   mp,gxy,%d,%g' ,kmp,mat.gxy);
                fprintf(fid,'\n   mp,gyz,%d,%g' ,kmp,mat.gyz);
                fprintf(fid,'\n   mp,gxz,%d,%g' ,kmp,mat.gxz);
                fprintf(fid,'\n   mp,dens,%d,%g',kmp,mat.density);
            otherwise
                error('Unknown material type in database');
        end
        switch mat.type
            case {'isotropic','orthotropic'}
                % Note that entering a blank or a zero for XYCP,YZCP, or XZCP
                % triggers the default value of -1.0. To specify an effective zero,

                % use a small, nonzero value (such as 1E-14).
%                 if isequal(0,mat.xycp), mat.xycp=1e-14; end
%                 if isequal(0,mat.yzcp), mat.yzcp=1e-14; end
%                 if isequal(0,mat.xzcp), mat.xzcp=1e-14; end
                % convert degrees to radians
                %mat.alp0 = mat.alp0 * pi/180;
                % read all of the failure criteria values
                
                
%                 for kfc = 1:numel(fcfields)
%                     fcname = fcfields{kfc}
%                     fcvalues{kfc} = mat.(fcname)
%                 end
%                 if all(cellfun('isempty',fcvalues))
%                     % do not print anything if all failure criteria
%                     % properties are empty
%                 else
%                     
%                     tempArray=zeros(9,1);
                    
                    %Tesile Properties
    
                    uts=mat.uts;
                    nStrenghts=length(uts);
                    if nStrenghts <3
                        uts=fullyPopluateStrengthsArray(uts);
                    end
                    
                    ucs=mat.ucs;
                    nStrenghts=length(ucs);
                    if nStrenghts <3
                        ucs=fullyPopluateStrengthsArray(ucs);
                    end
                    
                    uss=mat.uss;
                    nStrenghts=length(uss);
                    if nStrenghts <3
                        uss=fullyPopluateStrengthsArray(uss);
                    end


                    fprintf(fid,'\n   tb,fcli,%d,1,20,1',kmp);
                    fprintf(fid,'\n   tbdata,1,%g,%g,%g,%g,%g,%g',uts(1),ucs(1),uts(2),ucs(2),uts(3),ucs(3));
                    fprintf(fid,'\n   tbdata,7,%g,%g,%g,,,',uss(1),uss(2),uss(3));
                    fprintf(fid,'\n   tbdata,13,%g,%g,%g,%g',mat.xzit,mat.xzic,mat.yzit,mat.yzic);
                    fprintf(fid,'\n   tbdata,17,%g,%g,%g,%g',mat.g1g2,mat.etal,mat.etat,mat.alp0);
%                     for kf = 1:numel(fcvalues)
%                         if ~isempty(fcvalues{kf})
%                             fprintf(fid,'\n   tbdata,%d,%g',kf,fcvalues{kf});
%                         end
%                     end
                %end
            otherwise
                error('Unknown material type in database');
        end
        fprintf(fid,'\n');
    end
    
    fprintf(fid,'\n! WRITE THE COMPOSITE LAYUPS =================================\n');
    rCounter=1;
    switch fea.ansys.ElementSystem
        %tcl:  first few lines are same for shell91 and shell99
        case {'281','181'}
            
            %Outer AeroShell
            [nStationLayups, nStations] = size(blade.stacks);
            maxSectionNumber=str2num([int2str(nStations) int2str(nStationLayups)]); %Used to start web section ids
            for iStation=1:nStations
                for iPerimeter=1:nStationLayups
                    secID=str2num([int2str(iStation) int2str(iPerimeter)]);
                    currentStack=blade.stacks(iPerimeter,iStation);
                    fprintf(fid,'\n   ! %s',currentStack.name);
                    fprintf(fid,'\n   sectype,%d,shell'  ,secID);
                    for iLayer=1:length(currentStack.plygroups)
                        currentLayer=currentStack.plygroups(iLayer);
                        layerThickness=currentLayer.nPlies*currentLayer.thickness/1000; %Divide by 1000 to convert to m.
                        layerLayupAngle=currentLayer.angle;
                        matID=currentLayer.materialid;
                        fprintf(fid,'\n      secdata,%g,%d,%g,,',layerThickness,matID,layerLayupAngle);
                    end
                    fprintf(fid,'\n   secoffset,bot\n');
                end
            end
            %Web(s)
            nWebs=length(blade.swstacks);
            
            %The following two lines help make unique IDs for web sections
            %based on the highes section already defined for aeroshell
            orderOfMagnitude=floor( log10(maxSectionNumber));
            webSectionIDstart=ceil(maxSectionNumber/10^orderOfMagnitude)*10^orderOfMagnitude;
            for iWeb=1:nWebs
                [~, nStations] = size(blade.swstacks{iWeb});
                for iStation=1:nStations

                    currentStack=blade.swstacks{iWeb}(iStation);
                    if ~isempty(currentStack.plygroups)
                        secID=webSectionIDstart+iStation+(iWeb-1)*10^orderOfMagnitude;
                        fprintf(fid,'\n   ! %s',currentStack.name);
                        fprintf(fid,'\n   sectype,%d,shell'  ,secID);
                        for iLayer=1:length(currentStack.plygroups)
                            currentLayer=currentStack.plygroups(iLayer);
                            layerThickness=currentLayer.nPlies*currentLayer.thickness/1000; %Divide by 1000 to convert to m.
                            layerLayupAngle=currentLayer.angle;
                            matID=currentLayer.materialid;
                            fprintf(fid,'\n      secdata,%g,%d,%g,,',layerThickness,matID,layerLayupAngle);
                        end
                        fprintf(fid,'\n   secoffset,mid\n');
                    end

                end
            end
        
        otherwise
            errordlg(sprintf('Element System %s not yet available',fea.ansys.ElementSystem),'write_shell7 error')
            error('Element System %s not yet available',fea.ansys.ElementSystem);
            
    end
    
    % [~,jobtitle,~] = fileparts(blade.job_name);
    fprintf(fid,'\n/title,%s',fea.ansys.dbname);
    fprintf(fid,'\nZrCount=%d\n',rCounter);

    %tcl: DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES
    %tcl:    THE LINES ARE PRODUCED WITH THREE DIFFERENT SPLINING MACROS    
    fprintf(fid,'\n! DEFINE KEYPOINTS FOR SECTIONS AND CONNECT KEYPOINTS WITH LINES\n');
    
    % Create a coordinate system roughly in the fiber direction (+X down blade, +Z up toward LP side)
    % --> beginning with global csys, rotate -90 about local Z, then -90 about local Y
    fprintf(fid,'\nlocal,1000,CART,0,0,0, -90,0,-90\n');
    twistFlag = 1;  % ccw rotor rotation
    if blade.rotorspin == 1
        twistFlag = -1;  % cw rotor rotation
    end
    for kStation = 1:TotalStations

        % use the generating line to translate and rotate the coordinates
        % presweep ===========================================================
        if all(blade.isweep==0)
            table = [0, 0, nan];
        else
            N = length(blade.ispan);
            table = [blade.ispan(:), blade.isweep(:), nan(N,1)];
        end
        blade_struct.PresweepRef.method = 'shear';
        blade_struct.PresweepRef.table = table;
        blade_struct.PresweepRef.pptype = 'spline';
        % precurve ===========================================================
        if all(blade.iprebend==0)
            table = [0, 0, nan];
        else
            N = length(blade.ispan);
            table = [blade.ispan(:), blade.iprebend(:), nan(N,1)];
        end
        blade_struct.PrecurveRef.method = 'shear';
        blade_struct.PrecurveRef.table = table;
        blade_struct.PrecurveRef.pptype = 'spline';
        
        blade_struct = calcGenLinePP(blade_struct);
        
        presweep_slope = ppval(blade_struct.PresweepRef.dpp,blade.ispan(kStation));
        precurve_slope = ppval(blade_struct.PrecurveRef.dpp,blade.ispan(kStation));
        presweepDeg = 180/pi*atan(presweep_slope*twistFlag);
        precurveDeg = 180/pi*atan(-precurve_slope);
        [presweep_rot, precurve_rot] = deal(0);
        if isequal(blade_struct.PresweepRef.method,'normal')
            presweep_rot = atan(presweep_slope*twistFlag);
        end
        if isequal(blade_struct.PrecurveRef.method,'normal')
            precurve_rot = atan(-precurve_slope);
        end
        transX = twistFlag*ppval(blade_struct.PresweepRef.pp,blade.ispan(kStation));
        transY = ppval(blade_struct.PrecurveRef.pp,blade.ispan(kStation));

        % ensure we are in csys0 and no keypoints are selected
        fprintf(fid,'\ncsys,0');

        % Create a coordinate system to be used later for aligning the fiber direction.
        % First, load the csys defined earlier (+X down blade, +Z up toward LP side)
        fprintf(fid,'\n   csys,1000');
        % Next, translate & rotate relative to this active csys (use CLOCAL, not LOCAL)
        % translation: global X,Y,Z => local y,z,x
        % rotation: presweep is local z rotation & precurve is local y rotation
        fprintf(fid,'\n   clocal,%d,CART,%g,%g,%g, %g,%g,%g\n',(1000+kStation),...
            blade.ispan(kStation),transX,transY, presweepDeg,0,precurveDeg);
        
        % Create coordinate system at the tip
        fprintf(fid,'\nlocal,12,CART,%g,%g,%g, %g,%g,%g\n',...
            transX,transY,blade.ispan(kStation), 0,precurve_rot*180/pi,presweep_rot*180/pi);
    end
    fprintf(fid,'\n   csys,0');
    fprintf(fid,'\nksel,all\n');
    
    %jcb: as of 2011-05-26, the keypoints are transformed directly
    %     rather than with ansys commands
%     fprintf(fid,'\n! ROTATE SECTIONS ======================================\n');
%     fprintf(fid,'\ncsys,1');
%     for kStation = 1:TotalStations
%         fprintf(fid,'\n   lsel,s,loc,z,%g',data.station(kStation).LocationZ);
%         fprintf(fid,'\n   lgen,2,all,,,,%g,,,,1\n',twistFlag*data.station(kStation).DegreesTwist);
%     end
%     fprintf(fid,'\ncsys,0');
    fprintf(fid,'\nallsel\n');
    

    
    blade.mesh=0.45;
    [nodes,elements,outerShellElSets,shearWebElSets]=blade.getShellMesh();

    [nnodes,~]=size(nodes);
    [nelements,~]=size(elements);
    fprintf(fid,'\n! DEFINE NODES =======================================\n');
    for iNode=1:nnodes
        fprintf(fid,'n, %i, %f, %f, %f\n', iNode,nodes(iNode,1),nodes(iNode,2),nodes(iNode,3));
    end

    
    %Set the elemetn Type
    switch fea.ansys.ElementSystem
        case '281'
            fprintf(fid,'type, 11\n');
        case '181'  
            fprintf(fid,'type, 12\n');
        otherwise
            errordlg(sprintf('Element System %s not yet available',fea.ansys.ElementSystem),'write_shell7 error')
            error('Element System %s not yet available',fea.ansys.ElementSystem);
    end
    dup=[];
    fprintf(fid,'\n! DEFINE ELEMENTS =======================================\n');
    for iElement=1:nelements
        if numel(unique(elements(iElement,:)))==4
        fprintf(fid,'e, %i, %i, %i, %i  !Element %i \n', elements(iElement,1),elements(iElement,2),elements(iElement,3),elements(iElement,4),iElement);
        else
            dup=[dup; iElement];
        end
    end
    
    
    fprintf(fid,'\n! ASSIGN SECTIONS TO OUTER SHELL ELEMENTS =======================================\n');
    for iStation=1:nStations
        for iPerimeter=1:nStationLayups
            secID=str2num([int2str(iStation) int2str(iPerimeter)]);
            csID=1000+iStation;
            elementList=outerShellElSets(iPerimeter,iStation).elementList;
            for iEl=1:numel(elementList)
                fprintf(fid,'   emodif,%i,secnum,%i\n',elementList(iEl),secID);
                fprintf(fid,'   emodif,%i,esys,%i\n',elementList(iEl),csID);
            end 
        end
    end
    
    
    fprintf(fid,'\n! ASSIGN SECTIONS TO SHEARWEB(S) SHELL ELEMENTS =======================================\n');
    nWebs=length(blade.swstacks);
    for iWeb=1:nWebs
        [~, nStations] = size(blade.swstacks{iWeb});
        for iStation=1:nStations
            currentStack=blade.swstacks{iWeb}(iStation);
            if ~isempty(currentStack.plygroups)
                secID=webSectionIDstart+iStation+(iWeb-1)*10^orderOfMagnitude;
                csID=1000+iStation;
                
                elementList=shearWebElSets{iWeb}(iStation).elementList;
                for iEl=1:numel(elementList)
                    fprintf(fid,'   emodif,%i,secnum,%i\n',elementList(iEl),secID);
                    %fprintf(fid,'   emodif,%i,esys,%i\n',elementList(iEl),csID);
                end 
            end

        end
    end
%     %tcl: reverse area normals if clockwise blade
%     %tcl:    shear web areas are reversed as well - not necessary, just easier
%     if blade.rotorspin == 1 % clockwise rotation
%         fprintf(fid,'\n   areverse,all');
%     end
    fprintf(fid,'\n   ENSYM,,,,1,%i',nelements);
    
    %jcb: are these 2 lines necessary now that we have local coordinate
    %  systems to deal with presweep and precurve?
    fprintf(fid,'\n   local,11,CART,0,0,0,90,0,-90');
    fprintf(fid,'\n   esys,11');
      
    
    %LocationZ_lastStation = data.station(TotalStations).LocationZ;
    %fprintf(fid,'\n   local,12,cart,0,0,%f,0,0,0',LocationZ_lastStation);
    fprintf(fid,'\n   csys,12');
    fprintf(fid,'\n   nsel,none');
    fprintf(fid,'\n   n,,0.0,0.0,0.0');
    fprintf(fid,'\n   *get,z_master_node_number,node,,num,max');
    fprintf(fid,'\n   type,21');
    fprintf(fid,'\n   real,999');
    fprintf(fid,'\n   e,z_master_node_number');
    fprintf(fid,'\n   nsel,all');
    fprintf(fid,'\n   csys,0');
    fprintf(fid,'\n   allsel\n');

%jcb: TO BE DELETED
%     fprintf(fid,'\n   nsel,s,loc,z,%.7f,%.7f\n',...
%         (LocationZ_lastStation - 0.0000001),...
%         (LocationZ_lastStation + 0.0000001));

    % select tip station lines and then nodes attached to those lines
    %jcb: I think this can be cleaned up by moving these after the
    %     'e,z_master_node_number' command above and changing 's' to 'a'
    %     below ('nsel' command is then unnecessary because z_master_node 
    %     is already selected)
    %fprintf(fid,'\n   cmsel,s,tip_station_lines');
    fprintf(fid,'\n   nsll,s,1');
    fprintf(fid,'\n   nsel,a,node,,z_master_node_number');
    
    
    switch fea.ansys.ElementSystem
        case {'91','99','281','181'}
            fprintf(fid,'\n   cerig,z_master_node_number,all,RXYZ\n');
        case '191'
            fprintf(fid,'\n   cerig,z_master_node_number,all,uxyz\n');
    end
               
    if isequal(fea.ansys.BoundaryCondition,'cantilevered')
        %jcb: FIXME - nsel could break with swept/bent blades
        fprintf(fid,'\n   nsel,s,loc,z,0');
        fprintf(fid,'\n   d,all,all');
        fprintf(fid,'\n   nsel,all\n');
    end
    
    fprintf(fid,'\nallsel');
    fprintf(fid,'\n!   nummrg,all');
    fprintf(fid,'\n!   numcmp,node');
    fprintf(fid,'\ncsys,0\n');
    
    %%% Material Properties %%%
    fprintf(fid,'mpwrite,Materials,txt,,\n');
    %if ~all(cellfun('isempty',fcvalues))
        fprintf(fid,'/output,Strengths,txt,,\n');
        fprintf(fid,'TBLIST, ,ALL\n');
        fprintf(fid,'/output\n');
    %end
    
    % enter POST1 for postprocessing configuration commands
    fprintf(fid,'\nfinish');
    fprintf(fid,'\n/post1\n');   
    
    fprintf(fid,'\nfctyp,dele,all   ! remove all material failure-criteria postprocessing\n');
    for kfc = 1:length(fea.ansys.FailureCriteria)
        if fea.ansys.FailureCriteria{kfc,2}
            fprintf(fid,'fctyp,add,%s\n',fea.ansys.FailureCriteria{kfc,1});
        end
    end
    fprintf(fid,'\nfinish\n')
    %%% Material Properties %%%
    fprintf(fid,'mpwrite,Materials,txt,,\n');
    if ~all(cellfun('isempty',fcvalues))
        fprintf(fid,'/output,Strengths,txt,,\n');
        fprintf(fid,'TBLIST, ,ALL\n');
        fprintf(fid,'/output\n');
    end
    %%% Section Properties %%%
    fprintf(fid,'/output, Sections,txt\n');
    fprintf(fid,'SLIST,,,,FULL\n');
%     fprintf(fid,'SLIST,\n');
    fprintf(fid,'/output\n');

    %%% Element Properties %%%
    fprintf(fid,'/output, Elements,txt\n');
    fprintf(fid,'elist,all,,,0,0 \n');
    fprintf(fid,'/output\n');
    
    %%% Make NLIST file %%%
    fprintf(fid,'!!! BEGIN MAKE_NLIST MACRO TEXT\n');
    fprintf(fid,'ESEL,S,SEC,,1,%i   \n',webSectionIDstart-1);
    % fprintf(fid,'ALLSEL   \n');
    fprintf(fid,'NSLE,S  \n');
    % fprintf(fid,'nsel,u,node,,z_master_node_number\n');
    fprintf(fid,'/output,NLIST,lis\n');
    fprintf(fid,'/page,1e6,,1e6,,\n');
    fprintf(fid,'NLIST,ALL, , ,XYZ,NODE,NODE,NODE\n');
    fprintf(fid,'/output,\n');
    fprintf(fid,'ALLSEL,ALL \n');

    fprintf(fid,'!!! END MAKE_NLIST TEXT\n');
    
    % save database file
    fprintf(fid,'\nfinish');
    fprintf(fid,'\nsave');

%    # unit tip load suite
%    if {$BPEstatus == 1} {
%       append shell7Contents "    
%    ...
%    ...


catch ME
    % The try..catch..end ensures the file gets closed
    % in case of a programming error.
    fclose(fid);
    rethrow(ME);
end
fclose(fid);
msg = sprintf('shell7.src written to %s',blade.paths.job);
if 1%ble: can make this smart: originally --> app.batchrun
    disp(msg)
else
    msgbox(msg,'Notification');
end

end
function strengthArray=fullyPopluateStrengthsArray(strengthArray)
    nStrenghts=length(strengthArray);
    if nStrenghts <3
        for i=1:(3-nStrenghts)
            strengthArray=[strengthArray strengthArray(i)];
        end
    end
end