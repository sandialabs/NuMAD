function obj = YAML_to_BladeDef(obj,file)
% Convert YAML to prenumad blade object, Chris Kelley, in support of BAR Project, 4/18/19 %
% helper functions from Helena Canet at TU Munich <helena.canet@tum.de>

% Read yaml file from NREL
data = ReadYaml(file);
% folder name for blade object
out_folder = ['yaml2BladeDef_' strrep(file,'.yaml','')];
numad_af_folder = fullfile(out_folder,'airfoils/');
% Blade object output filename
blade_name = [ out_folder '/' strrep(file,'.yaml','') '_blade.mat'];
matdb_name = [ out_folder '/MatDBsi.txt'];
numad_name = [ out_folder '/' strrep(file,'.yaml','') '_numad.nmd'];
% Blade radius in meters (not in exactly in yaml file)
L = ceil(data.components.blade.outer_shape_bem.reference_axis.z.values{end});
R = L + data.components.hub.outer_shape_bem.diameter./2;


mkdir(out_folder)
mkdir([out_folder '/af_coords/'])
% mkdir([out_folder '/af_polars/'])
mkdir(numad_af_folder)

% create blade object
L = R - data.components.hub.outer_shape_bem.diameter./2;
obj.ispan = cell2mat(data.components.blade.outer_shape_bem.chord.grid)'.*L;

%% Aerodynamic properties
% using interp because yaml can have different r/R for twist and chord
obj.degreestwist = interp1(cell2mat(data.components.blade.outer_shape_bem.twist.grid)'.*L,cell2mat(data.components.blade.outer_shape_bem.twist.values)',obj.ispan).*180./pi;
obj.chord = interp1(cell2mat(data.components.blade.outer_shape_bem.chord.grid)'.*L,cell2mat(data.components.blade.outer_shape_bem.chord.values)',obj.ispan);

for i = 1:length(data.airfoils)
    af_dir_names{i} = data.airfoils{i}.name;
end

for i = 1:length(data.components.blade.outer_shape_bem.airfoil_position.labels)
    [C,IA,IAF(i)] = intersect(data.components.blade.outer_shape_bem.airfoil_position.labels{i},af_dir_names,'stable');
    tc(i) = data.airfoils{IAF(i)}.relative_thickness;
    tc_xL(i) = data.components.blade.outer_shape_bem.airfoil_position.grid{i};
    aero_cent(i) = data.airfoils{IAF(i)}.aerodynamic_center;
    xf_coords = [cell2mat(data.airfoils{IAF(i)}.coordinates.x)',cell2mat(data.airfoils{IAF(i)}.coordinates.y)'];
    
    % find coordinate direction (clockwise or counter-clockwise) Winding
    % Number. clockwise starting at (1,0) is correct
    if mean(gradient(atan(xf_coords(:,2)./xf_coords(:,1))), 'omitnan') < 0 % clockwise
        dir = 1;
    end
    if mean(gradient(atan(xf_coords(:,2)./xf_coords(:,1))), 'omitnan') > 0 % counter-clockwise
        dir = -1;
    end
        if dir == -1
            xf_coords = flipud(xf_coords);
        end
%     % add (1,0) at TE if necessary
%     if xf_coords(1,1) ~= 1
%         xf_coords = [1,0;xf_coords];
%     end    
%     if xf_coords(1,2) ~= 0
%         xf_coords = [1,0;xf_coords];
%     end
%     if xf_coords(end,1) ~= 1
%         xf_coords = [xf_coords;1,0];
%     end
%     if xf_coords(end,2) ~= 0
%         xf_coords = [xf_coords;1,0];
%     end
writeNuMADAirfoil(xf_coords,data.components.blade.outer_shape_bem.airfoil_position.labels{i},[out_folder '/af_coords/' data.components.blade.outer_shape_bem.airfoil_position.labels{i} '.txt'])
end
obj.span = obj.ispan;
obj.percentthick = interp1(cell2mat(data.components.blade.outer_shape_bem.airfoil_position.grid).*L,tc,obj.ispan).*100;
obj.aerocenter = interp1(cell2mat(data.components.blade.outer_shape_bem.airfoil_position.grid).*L,aero_cent,obj.span);
obj.chordoffset = interp1(cell2mat(data.components.blade.outer_shape_bem.pitch_axis.grid)'.*L,cell2mat(data.components.blade.outer_shape_bem.pitch_axis.values)',obj.span);
obj.naturaloffset = 0;
obj.prebend = interp1(cell2mat(data.components.blade.outer_shape_bem.reference_axis.x.grid)'.*L,cell2mat(data.components.blade.outer_shape_bem.reference_axis.x.values)',obj.span);
obj.sweep = interp1(cell2mat(data.components.blade.outer_shape_bem.reference_axis.y.grid)'.*L,cell2mat(data.components.blade.outer_shape_bem.reference_axis.y.values)',obj.span);

%
for  i = 1:length(tc)
    afc(i) = AirfoilDef(fullfile([out_folder '/af_coords/' data.components.blade.outer_shape_bem.airfoil_position.labels{i} '.txt']));
    obj.addStation(afc(i),tc_xL(i).*L);
end
afc.resample; % max 175 points for Ansys FEA
% afc.plot
% 
% 
% 
% %Write airfoil polars for AeroDyn
% Xq = fliplr(min(tc):0.001:max(tc)); % thickness of interpolated airfoils
% %afpi = afp.interp1(tc,Xq);  % perform interpolation
% afci = afc.interp1(tc,Xq,0,50);  % perform interpolation
% 
% afci.addModOpts('ModType','DynStall','NegStallCn',-0.8);
% afci.updateModData;
% figure
% afci.plotMod('alpha','CL')
% 
% % Write the WTPerf / AeroDyn polars
% opts = struct();
% opts.title = {'';''};
% opts.UseCM = true;
% for k=1:numel(Xq) % Write AeroDyn files if they don't exist
%     if exist([results_dir sprintf('af_%04d.txt',round(Xq(k)*1e3))],'file') == 2
%     else
%         opts.title = {sprintf('Interpolated thickness: %g',Xq(k));''};
%         polar = afci(k).writePolar(opts);
%         filename = sprintf('af_%04d.txt',round(Xq(k)*1e3));
%         filename = fullfile(results_dir,filename);
%         writeAeroDynPolar(polar,filename,'v13');
%     end
% end
%%


%% Material Properties
for i = 1:length(data.materials)
    mat(i).name = data.materials{i}.name;
    if data.materials{i}.orth == 1
        mat(i).type = 'orthotropic';
    else
        mat(i).type = 'isotropic';
    end
    % Add ply thickness option if ply thickness exists in yaml
    try
        mat(i).layerthickness = data.materials{i}.ply_t*1000;
    catch
        disp(['Warning! material ply thickness ' data.materials{i}.name  ' not defined, assuming 1 mm thickness'])
        mat(i).layerthickness = 1;
    end

    if size(data.materials{i}.E,2) == 1  % isotropic
        mat(i).ex = data.materials{i}.E(1);
        mat(i).ey = data.materials{i}.E(1);
        mat(i).ez = data.materials{i}.E(1);
        mat(i).gxy = data.materials{i}.G(1);
        mat(i).gxz = data.materials{i}.G(1);
        mat(i).gyz = data.materials{i}.G(1);
        mat(i).prxy = data.materials{i}.nu(1);
        mat(i).prxz = data.materials{i}.nu(1);
        mat(i).pryz = data.materials{i}.nu(1);
        mat(i).uts = data.materials{i}.Xt;
        mat(i).ucs = -data.materials{i}.Xc;
        mat(i).uss = data.materials{i}.S;
        mat(i).xzit = 0.3;
        mat(i).xzic = 0.25;
        mat(i).yzit = 0.3;
        mat(i).yzic = 0.25;
        mat(i).g1g2 = data.materials{i}.GIc/data.materials{i}.GIIc;
        mat(i).alp0 = data.materials{i}.alp0;
        mat(i).etat = LARCetaT(mat(i).alp0);
        mat(i).etal = LARCetaL(mat(i).uss,mat(i).ucs,mat(i).alp0);
    else % orthotropic
        mat(i).ex = data.materials{i}.E{1};
        mat(i).ey = data.materials{i}.E{2};
        mat(i).ez = data.materials{i}.E{3};
        mat(i).gxy = data.materials{i}.G{1};
        mat(i).gxz = data.materials{i}.G{2};
        mat(i).gyz = data.materials{i}.G{3};
        mat(i).prxy = data.materials{i}.nu{1};
        mat(i).prxz = data.materials{i}.nu{2};
        mat(i).pryz = data.materials{i}.nu{3};
        mat(i).uts = cell2mat(data.materials{i}.Xt);
        mat(i).ucs = -cell2mat(data.materials{i}.Xc);
        mat(i).uss = cell2mat(data.materials{i}.S);
        mat(i).xzit = 0.3;
        mat(i).xzic = 0.25;
        mat(i).yzit = 0.3;
        mat(i).yzic = 0.25;
        mat(i).g1g2 = data.materials{i}.GIc/data.materials{i}.GIIc;
        mat(i).alp0 = data.materials{i}.alp0;
        mat(i).etat = LARCetaT(mat(i).alp0);
        mat(i).etal = LARCetaL(mat(i).uss(1),mat(i).ucs(2),mat(i).alp0);
    end
    try 
        mat(i).m=data.materials{i}.m;
    catch
        warning('No fatigue exponent found for material: %s',data.materials{i}.name)
    end
    mat(i).density = data.materials{i}.rho;
    mat(i).dens = data.materials{i}.rho;
    mat(i).drydensity = data.materials{i}.rho;
    
    if isfield(data.materials{i},'description') && isfield(data.materials{i},'source')
        desc_sourc = [data.materials{i}.description ', ' data.materials{i}.source];
        mat(i).reference = desc_sourc;
    else
        mat(i).reference = [];
    end
    obj.addMaterial(mat(i));
end
obj.updateGeometry;

%%


%% Blade Components
N_layer_comp = length(data.components.blade.internal_structure_2d_fem.layers); % number of layer based components (gelcoat, sparcaps, skins, etc)
N_sweb_comp = length(data.components.blade.internal_structure_2d_fem.webs); % number of shear web components (gelcoat, sparcaps, skins, etc)
N_comp = length(data.components.blade.internal_structure_2d_fem.layers)+ length(data.components.blade.internal_structure_2d_fem.webs);


% Spar Cap Width and Offset
for i = 1:N_layer_comp
   if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'spar')'
        is_spar(i) = 1;
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.side]),'suc')'
        is_lpspar(i) = 1;
        end
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.side]),'pres')'
        is_hpspar(i) = 1;
        end
   else
       is_spar(i) = 0;
   end
end
layer_vec = 1:N_layer_comp;
I_spar = layer_vec(logical(is_spar));
I_spar_hp = layer_vec(logical(is_hpspar));
I_spar_lp = layer_vec(logical(is_lpspar));
% Because spar cap width must be constant, average the yaml file on
% pressure and suction surfaces across span
obj.sparcapwidth(1) =  mode(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_spar_hp}.width.values)).*1000;
obj.sparcapwidth(2) =  mode(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_spar_lp}.width.values)).*1000;
obj.sparcapoffset(1) = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_spar_hp}.offset_y_pa.values)).*1000;
obj.sparcapoffset(2) = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_spar_lp}.offset_y_pa.values)).*1000;

% TE and LE Bands
for i = 1:N_layer_comp
   if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'reinf')'
        is_reinf(i) = 1;
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'le')'
        LE(i) = 1;
        TE(i) = 0;
        else 
            LE(i) = 0;
            TE(i) = 1;
        end
   else
       is_reinf(i) = 0;
   end
end
layer_vec = 1:N_layer_comp;
web_vec = 1:N_sweb_comp;
I_reinf = layer_vec(logical(is_reinf));
I_LE = layer_vec(logical(LE));
I_TE = layer_vec(logical(TE));

% Leading and Trailing Edge bands are constants in millimeters
% LE Band
obj.leband = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_LE}.width.values)).*1000/2; %EC Divide by 2: Half of the le width to go on LP/HP surfaces respectively
% TE Band
obj.teband = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I_TE}.width.values)).*1000/2; %EC Divide by 2: Half of the te width to go on LP/HP surfaces respectively



% Create blade component structure (comp(i)) for addition to blade object
for i = 1:N_layer_comp
    comp(i).group = 0;
    comp(i).name = data.components.blade.internal_structure_2d_fem.layers{i}.name;
%   comp.material = data.components.blade.internal_structure_2d_fem.layers{i}.material;
    [C,IA,IB] = intersect({obj.materials.name},data.components.blade.internal_structure_2d_fem.layers{i}.material,'stable');
    comp(i).materialid = IA;
    try
        comp(i).fabricangle = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.fiber_orientation.values));
    catch
        comp(i).fabricangle = 0;
    end
    
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'spar')
    comp(i).imethod = 'pchip';
    else
     comp(i).imethod = 'linear';   
    end
    comp(i).cp(:,1) = cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.thickness.grid)';
    temp_n_layer = cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.thickness.values)'.*1000./mat(comp(i).materialid).layerthickness;
    I_round_up = find(temp_n_layer>0.05 & temp_n_layer< 0.5);
    comp(i).cp(:,2) = round(cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.thickness.values)'.*1000./mat(comp(i).materialid).layerthickness); % use for rounding n layers where each ply is not necessarily 1 mm
    if ~isempty(I_round_up)
        comp(i).cp(I_round_up,2) = 1; % increase n_layers from 0 to 1 for 0.05<n_layers<0.5
    else
    end
%     comp(i).cp(:,2) = cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.thickness.values)'.*1000;  % use when each material ply is 1 mm
    comp(i).pinnedends = 0;
    clear I_round_up
end


% uv coating
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'uv')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
comp(I).hpextents = {'le'  'te'};
comp(I).lpextents = {'le'  'te'};
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I

% Shell skin1
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'shell_skin_outer')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
comp(I).hpextents = {'le'  'te'};
comp(I).lpextents = {'le'  'te'};
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I


% Spar Caps (pressure and suction)
comp(I_spar_hp).hpextents = {'b'  'c'};
comp(I_spar_lp).lpextents = {'b'  'c'};


% LE Band
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'le_reinf')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
comp(I(i)).hpextents = {'le'  'a'};
comp(I(i)).lpextents = {'le'  'a'};
end 
obj.leband = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I}.width.values)).*1000/2; %EC Why assign obj.leband again? Divide by 2: Half of the le width to go on LP/HP surfaces respectively
clear I

% TE Band
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'te_reinf')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
comp(I(i)).hpextents = {'d'  'te'};
comp(I(i)).lpextents = {'d'  'te'};
end
obj.teband = mean(cell2mat(data.components.blade.internal_structure_2d_fem.layers{I}.width.values)).*1000/2; %EC Why assign obj.teband again? Divide by 2: Half of the te width to go on LP/HP surfaces respectively
clear I 


% Trailing edge suction surface panel
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'te_ss')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
								  
comp(I(i)).lpextents = {'c'  'd'};
end
clear I

% Leading edge suction surface panel
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'le_ss')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
								  
comp(I(i)).lpextents = {'a'  'b'};
end
clear I

% Leading edge pressure surface panel
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'le_ps')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
comp(I(i)).hpextents = {'a'  'b'};
								  
end
clear I

% Trailing edge pressure surface panel
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'te_ps')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
for i = 1:length(I)
comp(I(i)).hpextents = {'c'  'd'};
								  
end
clear I

% Shell skin2
clear I
for i = 1:N_layer_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'shell_skin_inner')'
        I(i) = 1;
    else
        I(i) = 0;
    end    
end
I = layer_vec(logical(I));
comp(I).hpextents = {'le'  'te'};
comp(I).lpextents = {'le'  'te'};
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I


% Forward Shear Web Skin1
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'fore')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'skin_le')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'fore')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000;% shear web offset from reference axis
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z+' sw_offset]};
% comp(I).lpextents = {['z+' sw_offset]};
comp(I).hpextents = {'b'};
comp(I).lpextents = {'b'};
comp(I).group = 1;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I I_web

% Forward Shear Web Filler
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'fore')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'filler')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'fore')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000; % shear web offset from max thickness
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z+' sw_offset]};
% comp(I).lpextents = {['z+' sw_offset]};
comp(I).hpextents = {'b'};
comp(I).lpextents = {'b'};
comp(I).group = 1;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
clear I I_web

% Forward Shear Web Skin2
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'fore')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'skin_te')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'fore')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000; % shear web offset from max thickness
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z+' sw_offset]};
% comp(I).lpextents = {['z+' sw_offset]};
comp(I).hpextents = {'b'};
comp(I).lpextents = {'b'};
comp(I).group = 1;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I I_web


% Rear Shear Web Skin1
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'rear')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'skin_le')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'rear')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000; % shear web offset from max thickness
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z-' sw_offset]};
% comp(I).lpextents = {['z-' sw_offset]};
comp(I).hpextents = {'c'};
comp(I).lpextents = {'c'};
comp(I).group = 2;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I I_web

% Rear Shear Web Filler
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'rear')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'filler')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'rear')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000; % shear web offset from max thickness
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z-' sw_offset]};
% comp(I).lpextents = {['z-' sw_offset]};
comp(I).hpextents = {'c'};
comp(I).lpextents = {'c'};
comp(I).group = 2;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
clear I I_web

% Rear Shear Web Skin2
clear I I_web
for i = 1:N_layer_comp
    if isfield(data.components.blade.internal_structure_2d_fem.layers{i},'web')        
        if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.web]),'rear')'
            if strfind(lower([data.components.blade.internal_structure_2d_fem.layers{i}.name]),'skin_te')'
        I(i) = 1;
    else
        I(i) = 0;
            end   
        end
    end
end
for i = 1:N_sweb_comp
    if strfind(lower([data.components.blade.internal_structure_2d_fem.webs{i}.name]),'rear')'
        I_web(i) = 1;
    else
        I_web(i) = 0;
    end    
end
I = layer_vec(logical(I));
I_web = web_vec(logical(I_web));
sw_offset = mean(cell2mat(data.components.blade.internal_structure_2d_fem.webs{I_web}.offset_y_pa.values)).*1000; % shear web offset from max thickness
xs =  (obj.sparcapwidth(1)./2 + sw_offset)/obj.sparcapwidth(1); % location of shear web relative to sparcap arc, b--c
% comp(I).hpextents = {[num2str(xs,2) 'b-c']};
% comp(I).lpextents = {[num2str(xs,2) 'b-c']};
% comp(I).hpextents = {['z-' sw_offset]};
% comp(I).lpextents = {['z-' sw_offset]};
comp(I).hpextents = {'c'};
comp(I).lpextents = {'c'};
comp(I).group = 2;
comp(I).name = [data.components.blade.internal_structure_2d_fem.layers{I}.name];
% CK Change me when yaml is fixed!!!!
comp(I).cp(:,2) = comp(I).cp(:,2);
clear I I_web



%%

for i = 1:N_layer_comp
obj.addComponent(comp(i))
end

obj.updateGeometry;
obj.updateKeypoints;
obj.updateBOM;
save(blade_name);
BladeDef_to_NuMADfile(obj,numad_name,matdb_name,numad_af_folder);


end
