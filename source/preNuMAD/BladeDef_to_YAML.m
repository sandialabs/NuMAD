function obj = BladeDef_to_YAML(obj,file)
% Convert YAML to prenumad blade object, Chris Kelley, in support of BAR Project, 4/18/19 %
% helper functions from Helena Canet at TU Munich <helena.canet@tum.de>

% Read original yaml file from NREL
file_orig = strrep(file,'_mod','');
data = ReadYaml(file_orig);
data_mod = data;
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

% update spar cap width and offset
N = length(data.components.blade.internal_structure_2d_fem.layers{I_spar_hp}.width.values);
data.components.blade.internal_structure_2d_fem.layers{I_spar_hp}.width.values = repmat({obj.sparcapwidth(1)./1000},1,N);
data.components.blade.internal_structure_2d_fem.layers{I_spar_lp}.width.values = repmat({obj.sparcapwidth(2)./1000},1,N);
data.components.blade.internal_structure_2d_fem.layers{I_spar_hp}.offset_x_pa.values = repmat({obj.sparcapoffset(1)./1000},1,N);
data.components.blade.internal_structure_2d_fem.layers{I_spar_lp}.offset_x_pa.values = repmat({obj.sparcapoffset(2)./1000},1,N);

% Update control points to control thicknesses of all components
for i = 1:N_layer_comp
    data.components.blade.internal_structure_2d_fem.layers{i}.thickness.grid = obj.components(i).cp(:,1)';
    data.components.blade.internal_structure_2d_fem.layers{i}.thickness.values = obj.components(i).cp(:,2)'.*obj.materials(obj.components(i).materialid).layerthickness./1000;
end
% comp(i).cp(:,2) = round(cell2mat(data.components.blade.internal_structure_2d_fem.layers{i}.thickness.values)'.*1000./); % use for rounding n layers where each ply is not necessarily 1 mm

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
N = length(data.components.blade.internal_structure_2d_fem.layers{I_LE}.width.values);
data.components.blade.internal_structure_2d_fem.layers{I_LE}.width.values = repmat({obj.leband./1000.*2},1,N);
% TE Band
N = length(data.components.blade.internal_structure_2d_fem.layers{I_TE}.width.values);
data.components.blade.internal_structure_2d_fem.layers{I_TE}.width.values = repmat({obj.teband./1000.*2},1,N);



file_name = strrep(file,'.yaml','')
data.name = file_name;


%% Write geometry (data) to yaml
WriteYaml(file,data);

end
