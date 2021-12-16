function precomp = Blade2PreComp(blade,matdb)

% First, filter the full set of material properties and stacks
% contained in MatDBsi.txt down to only the materials and stacks that are
% used by the actual NuMAD model, as described by *.nmd.  fprintf(fid,'lays filter
% results to the screen and the saves final set of materials and stacks to
% Mat and Comp variables.

% Prep materials information
matids=[];
compids=[];
fid=fopen('PrepMat.txt','wt');

comp.list = {'**UNSPECIFIED**'};
for k=1:numel(blade.matdb)
    if isequal(blade.matdb(k).type,'composite')
        comp.list{end+1} = blade.matdb(k).name;
    end
end

for ii=2:length(comp.list)
    name=comp.list{ii};
    fprintf(fid,'****************\n');
    fprintf(fid,'%s\n',name);
    fprintf(fid,'        contains:\n');
    
    % find index to composite material
    for j=1:length(matdb)
        if strcmp(matdb(j).name,name)
            index1=j;
            break;
        end
    end
    compids=[compids index1];
    
    fprintf(fid,' LAYERS:\n');
    for j=1:length(matdb(index1).layer)
        name=matdb(index1).layer(j).layerName;
        fprintf(fid,'  %s\n',name);
        for k=1:length(matdb)
            if strcmp(matdb(k).name,name)
                index2=k;
                break;
            end
        end
        fprintf(fid,'      Material ID# %i\n',index2);
        matids=[matids index2];
    end
end

fclose(fid);

matkeepers=unique(matids);
compkeepers=unique(compids);

material=matdb(matkeepers);
composite=matdb(compkeepers);

% check that all matl properties are in place
for ii=1:length(material)
    % Applicable to the isotropic matls:
    if isempty(material(ii).ey)
        material(ii).ey=material(ii).ex;
    end
    if isempty(material(ii).gxy)
        material(ii).gxy=material(ii).ex/(2*(1+material(ii).nuxy));
    end
    if isempty(material(ii).prxy)
        material(ii).prxy=material(ii).nuxy;
    end        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

N_sections=length(blade.ispan);
Bl_length=blade.ispan(end);

for ii=1:N_sections
    precomp.station(ii).span=blade.ispan(ii)/Bl_length;
    precomp.station(ii).lepos=blade.xoffset(ii); %ble: these are different from NuMAD GUI values slightly
    precomp.station(ii).chlen=blade.ichord(ii);
    precomp.station(ii).twist=blade.idegreestwist(ii);
    precomp.station(ii).tetype=blade.TEtype{ii};
end

% Convert NuMAD Keypoints to PreComp airfoil coords and create shape files
% ble: save the airfoil coordinates ensuring that the region keypoints are
% saved as points in the data (non-dimensional)

for ii=1:N_sections
    x=blade.profiles(:,1,ii)';
    y=blade.profiles(:,2,ii)';
    % convert to chord-normalized coordinated of the airfoil nodes
    % the first node, a leading-edge node, must be (0,0)
    % find le node by finding the smallest x coordinate (should be zero)
    [~,le(ii)]=min(x);
    x=[x(le(ii):end) x(1:le(ii)-1)];
    y=[y(le(ii):end) y(1:le(ii)-1)];
    x=x-x(1);
    y(1)=0;
    
    % remove all TE points for flatback, all but one for sharp/round
    if strcmpi(blade.TEtype(ii),'flat')
        pointer=find(y==0);
        %remove trailing edge point(s) and replace later with shear web
        x=x([1:pointer(2)-1 pointer(end)+1:end]);
        y=y([1:pointer(2)-1 pointer(end)+1:end]);
    else % sharp or round        
        pointer=find(y==0);
        %remove trailing edge point(s) and replace later with shear web
        x=x([1:pointer(2) pointer(end)+1:end]);
        y=y([1:pointer(2) pointer(end)+1:end]);
    end
    % normalize chord to a value of 1
    r=1/max(x); 
    x=x*r;
    y=y*r;
    
    % check to see that points are not too close together
    % (could just resample the airfoil at this point using resampleAirfoil.m)
    newx=x(1);
    newy=y(1);
    
        
    for j=2:length(x)
        % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % change this so flatbacks can have multiple points specified at
        % the TE position of x=1
% %         %dist=sqrt( (x(j)-newx(end))^2 + (y(j)-newy(end))^2 );
        dist=sqrt( (x(j)-newx(end))^2 );

        if dist>0.015 % this worked for SDM2010 paper
            newx=[newx x(j)];
            newy=[newy y(j)];
        elseif x(j)==1 && newx(end)~=1 % ble
            newx=[newx(1:end-1) x(j)];
            newy=[newy(1:end-1) y(j)];
        elseif x(j)==1 && newx(end)==1 % ble
            newx=[newx x(j)];
            newy=[newy y(j)];
        end
        % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    end
    x=newx;
    y=newy;
    
    % The following are unnecessary if max of x is exactly 1?
    r=1/max(x);
    x=x*r;
    y=y*r;
    

    % make sure trailing edge is at 1,0 for sharp and round airfoils
    if ~strcmpi(blade.TEtype(ii),'flat')     
        pointer=find(x>(1-0.000001)&x<(1+0.000001));
        if numel(pointer)>1 && abs(y(pointer(1)) - y(pointer(2))) > 1e-10
            keyboard
            error('too many TE points')
        end
        y(pointer)=0;
    end
    
    % remove last point(s) to prevent errors in PreComp
    n_to_remove=1;
    x=x(1:end-n_to_remove);
    y=y(1:end-n_to_remove);
    
    % use this plot to check the airfoil shapes
    if 0
        figure(101)
        plot(x,y,'bo')
        axis equal
        title(['Section ' num2str(ii)])
        grid on
        pause(0.5)
%         pause
    end
    
    precomp.station(ii).x=x;
    precomp.station(ii).y=y;
end

% Set up Materials

for ii=1:length(material)
    precomp.material.id(ii)=ii;
    precomp.material.e1(ii)=material(ii).ex;
    precomp.material.e2(ii)=material(ii).ey;
    precomp.material.g12(ii)=material(ii).gxy;
    precomp.material.nu12(ii)=material(ii).prxy;
    precomp.material.dens(ii)=material(ii).dens;
    precomp.material.name(ii)={material(ii).name};
end

% Set up Layup information

for ii=1:length(composite)
    precomp.stack(ii).name={composite(ii).name};
    if ~strcmpi(composite(ii).thicknessType,'Constant')
        warning('PreComp is not set up for anything other than ''Constant'' thickness type')
    end
    precomp.stack(ii).n_layer=composite(ii).uniqueLayers;
    if ~strcmpi(composite(ii).symmetryType,'none')
        warning('PreComp is not set up for anything other than symmetry ''none''')
    end
    for j=1:precomp.stack(ii).n_layer
        precomp.stack(ii).layer(j).n=composite(ii).layer(j).quantity;
        precomp.stack(ii).layer(j).thickness=composite(ii).layer(j).thicknessA;
        precomp.stack(ii).layer(j).theta=composite(ii).layer(j).theta;
        precomp.stack(ii).layer(j).n=composite(ii).layer(j).quantity;
                
        % determine the material id for each layer
        for k=1:length(precomp.material.id)
            if strcmp(composite(ii).layer(j).layerName,precomp.material.name{k})
                precomp.stack(ii).layer(j).matid=k;
            end
        end
    end
end

% Set up Shear Web Layup information
% determine the composite id for each shear web
SW=blade.shearweb;
if ~isempty(SW)
    for j=1:length(SW)
        for k=1:length(precomp.stack)
            if strcmp(SW(j).Material,precomp.stack(k).name)
                SW(j).stackid=k;
            end
        end
    end
    
    % Find the ends of shear webs
    % If the end of one web is the beginning of another, then it isn't
    % counted
    non_terminal_web_id=[];
    for j=1:length(SW)
        sw_end_dps=SW(j).Corner([3 4]);
        sw_end_at=SW(j).EndStation;
        for k=1:length(SW)
            if k~=j
                if sw_end_at==SW(k).BeginStation
                    % getting here means that the station where the jth web ends at is a beginning station for another web
                    % need to find out whether the web is not the same continuous web
                    if sw_end_dps(1)==SW(k).Corner(2) && sw_end_dps(2)==SW(k).Corner(1)
                        % Note that this shear web should be counted as a
                        % terminal web
                        non_terminal_web_id=[non_terminal_web_id j];
                    end
                end
            end
        end
    end
    
    % set up indices to the SW's that are terminal (meaning non-nonterminal SW's)
    terminal_web_id=[];
    for ii=1:length(SW)
        if isempty(find(ii==non_terminal_web_id))
            terminal_web_id=[terminal_web_id ii];
        end
    end
else
    % what if SW doesn't exist (meaning there is no shear web)
    
end

% Set up station data
for ii=1:N_sections
    keycpos = blade.keycpos(:,ii);
    regionStack = {blade.stacks(:,ii).name}';          
    % check to make sure there is a leading edge dp in the section.  If
    % not, send an error message
    if isempty(find(keycpos==0))
        error('Each station must contain a material D.P. at the leading edge in order for NuMAD2PreComp to work.');
    end
    hp_pointer=find(keycpos<=0);
    hp_pointer=hp_pointer(end:-1:1);
    lp_pointer=find(keycpos>=0);
    perChord_hp=keycpos(hp_pointer)*-1;
    perChord_lp=keycpos(lp_pointer);
    N_divs_hp=length(perChord_hp)-1;
    N_divs_lp=length(perChord_lp)-1;
    mat_hp=regionStack(N_divs_hp:-1:1);
    mat_lp=regionStack(N_divs_hp+1:end);

    isflatback=strcmpi(blade.TEtype{ii},'flat');        
    
    if isflatback     % if flatback:
%         if ~strcmp(mat_lp(end),mat_hp(end))
%             warning(['Flatback isn''t defined of consistent material in section #' num2str(i) '.  Using material specification for bottom (HP) half of flatback.' ])
%         end
        mat_flat=mat_hp(end);
        mat_lp=mat_lp(1:end-1);
        mat_hp=mat_hp(1:end-1);
        perChord_lp=perChord_lp(1:end-1);
        perChord_hp=perChord_hp(1:end-1);
        N_divs_lp=length(perChord_lp)-1;
        N_divs_hp=length(perChord_hp)-1;
    end
    
    % Store in PreComp structure
    precomp.station(ii).hp.segmentloc=perChord_hp;
    precomp.station(ii).lp.segmentloc=perChord_lp;
    precomp.station(ii).hp.n_segments=N_divs_hp;
    precomp.station(ii).lp.n_segments=N_divs_lp;
    
    % determine stack id for the materials before saving them
    % LP Surface
    for j=1:length(mat_lp)
        for k=1:length(precomp.stack)
            if strcmp(mat_lp(j),precomp.stack(k).name)
                precomp.station(ii).lp.stackid(j)=k;
            end
        end
    end
    % HP Surface
    for j=1:length(mat_hp)
        for k=1:length(precomp.stack)
            if strcmp(mat_hp(j),precomp.stack(k).name)
                precomp.station(ii).hp.stackid(j)=k;
            end
        end
    end
    
    % if at tip of blade, assign same as next inboard station
    if ii==N_sections
        precomp.station(ii).hp.stackid=precomp.station(ii-1).hp.stackid;
        precomp.station(ii).lp.stackid=precomp.station(ii-1).lp.stackid;
    end
    
    % Do Shear web conversion stuff
    tmp=[];
    tmp2=[];
    if ~isempty(SW) % if there is a shearweb
        % Sort out Shear Web placement and material information
        for j=1:length(SW) % count number of shear webs at this station and find index to information
            if SW(j).BeginStation==ii
                tmp=[tmp j];
            end
        end
        
        % assign sw stack id's to appropriate stations
        for j=1:length(tmp)
            precomp.station(ii).sw.stackid(j)=SW(tmp(j)).stackid;
        end
        
        % find location of shearweb, fraction of chord, & perpendicular to chord
        for j=1:length(tmp)
            temp=SW(tmp(j)).Corner+1;
            chord_hp=keycpos(temp(2))*-1;
            chord_lp=keycpos(temp(1));
            precomp.station(ii).sw.chloc(j)=mean([chord_hp chord_lp]);
        end
        
        % identify any terminal shear webs
        for j=1:length(terminal_web_id)
            if SW(terminal_web_id(j)).BeginStation==ii
                disp(' ')
                disp(['NuMAD station #' num2str(ii) ' has the beginning of a terminal shear web with SW id #' num2str(terminal_web_id(j)) ])
                tmp2=[tmp2 terminal_web_id(j)];
            end
        end
        sw_end_id{ii+1}=tmp2;
        
    end  % end SW stuff, add in terminal shear web information later after this loop on stations
    
    % Deal with flatbacks as shear webs
    if isflatback
        precomp.station(ii).fb.flag=1;
        % find index to composite material        
        for k=1:length(precomp.stack)
            if strcmp(mat_flat,precomp.stack(k).name)
                precomp.station(ii).fb.stackid=k;                
                break;  % stop looping and searching
            end
        end
        % ble: <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        % if at tip of blade, assign same as next inboard station
        if ii==N_sections
            precomp.station(ii).fb.stackid=precomp.station(ii-1).fb.stackid;
        end
        % ble: >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    else
        precomp.station(ii).fb.flag=0;
        precomp.station(ii).fb.stackid=0;
    end    
end

% Now calculate more information for terminal webs
if ~isempty(SW)
    for ii=2:N_sections
        keycpos = blade.keycpos(:,ii);
        % figure out how many are present in the structure already
        try
            jj=length(precomp.station(ii).sw.stackid);
        catch
            jj=0;
        end
        
        % assign sw stack id's to appropriate stations
        for j=1:length(sw_end_id{ii})
            precomp.station(ii).sw.stackid(jj+j)=SW( sw_end_id{ii}(j) ).stackid;
        end
        
        % Calculate shear web locations for terminal webs
        tmp=sw_end_id{ii};
        for j=1:length(tmp)
            temp=SW(tmp(j)).Corner+1;
            chord_hp=keycpos(temp(3))*-1;
            chord_lp=keycpos(temp(4));
            precomp.station(ii).sw.chloc(j+jj)=mean([chord_hp chord_lp]);
        end
    end
end

% Make sure that webs are in order from LE to TE of blade
for ii=1:length(precomp.station)
    try
        tmp=[precomp.station(ii).sw.chloc ; precomp.station(ii).sw.stackid]';
        tmp=sortrows(tmp);
        precomp.station(ii).sw.chloc = tmp(:,1)';
        precomp.station(ii).sw.stackid = tmp(:,2)';
    catch
        % do nothing
    end
end


% pass precomp data structure into the PreComp file creation functions

writePreCompShape(precomp,'shape.inp')

writePreCompMaterial(precomp,'materials.inp')

writePreCompLayup(precomp,'layup.inp')

writePreCompInput(precomp,'input.pci')

fclose all;  % if not already, close all files


end