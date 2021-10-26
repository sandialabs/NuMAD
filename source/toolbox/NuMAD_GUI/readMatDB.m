function matdb = readMatDB(filename)
%READMATDB  Read a NuMAD material database file 
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
%   MATDB = readMatDB(FILENAME)
%     Read the NuMAD material database FILENAME (typically 'MatDBsi.txt')
%     and return structure MATDB.
%

%===== CREDITS & CHANGELOG ================================================
%Developed by Wind & Water Power Technologies, Sandia National Laboratories
%2010.12.09  JCB: first draft

%filename = 'MatDBsi.txt'; clear matdb;

% open the file for Reading in Text mode
fid = fopen(filename,'rt');
if (fid == -1)
    error('Could not open file "%s"',filename);
end
% read the entire file contents
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

% find the start and end of each material block
matS = regexp(filecontents, '<material>');
matE = regexp(filecontents, '</material>');
% adjust start and end indices to encompass only child elements
matS = matS + length('<material>');
matE = matE - 1;

% %%%%%%%%%%%%  PROGRAMMER'S NOTE %%%%%%%%%%%%%
% The following regular expression pattern matches any number of characters
% found between the opening and closing "reference" tags
%pattern = '<reference>(.*)</reference>';
%t = regexp(filecontents, pattern, 'tokens');
% t is a cell containing a cell array
% try
%     % jcb: is there a better way to extract the contents of t?
%     af.reference = cell2mat(t{1});
% catch me
%     af.reference = '';
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% process each material block
Nmat = length(matS);
% jcb: this structure must match the one in 'NuMAD_materials.m : cb_new()'
matdb(1:Nmat) = deal(struct('type',[],'name',[],'reference',[],...
    'dens',[],'nuxy',[],'ex',[],'ey',[],'ez',[],'gxy',[],'gyz',[],'gxz',[],...
    'prxy',[],'pryz',[],'prxz',[],'xten',[],'xcmp',[],'yten',[],'ycmp',[],...
    'zten',[],'zcmp',[],'xy',[],'yz',[],'xz',[],'xycp',[],'yzcp',[],'xzcp',[],...
    'xzit',[],'xzic',[],'yzit',[],'yzic',[],'g1g2',[],'etal',[],'etat',[],...
    'alp0',[],'thicknessType',[],'uniqueLayers',[],'symmetryType',[],'layer',[]));
for kmat = 1:Nmat
    matC = filecontents(matS(kmat):matE(kmat));
    
    matdb(kmat).type = regexp_str(matC,'<type>(.*)</type>');
    mattype = matdb(kmat).type;  % the material type is needed later
    matdb(kmat).name = regexp_str(matC,'<name>(.*)</name>');
    matdb(kmat).reference = regexp_str(matC,'<reference>(.*)</reference>');
    
    switch mattype
        case {'isotropic','orthotropic'}
            matdb(kmat).dens = regexp_dbl(matC,'<dens>(.*)</dens>');
            matdb(kmat).ex = regexp_dbl(matC,'<ex>(.*)</ex>');
            
            switch mattype
                case 'isotropic'
                    matdb(kmat).nuxy = regexp_dbl(matC,'<nuxy>(.*)</nuxy>');
                    
                case 'orthotropic'
                    matdb(kmat).ey = regexp_dbl(matC,'<ey>(.*)</ey>');
                    matdb(kmat).ez = regexp_dbl(matC,'<ez>(.*)</ez>');
                    matdb(kmat).gxy = regexp_dbl(matC,'<gxy>(.*)</gxy>');
                    matdb(kmat).gyz = regexp_dbl(matC,'<gyz>(.*)</gyz>');
                    matdb(kmat).gxz = regexp_dbl(matC,'<gxz>(.*)</gxz>');
                    matdb(kmat).prxy = regexp_dbl(matC,'<prxy>(.*)</prxy>');
                    matdb(kmat).pryz = regexp_dbl(matC,'<pryz>(.*)</pryz>');
                    matdb(kmat).prxz = regexp_dbl(matC,'<prxz>(.*)</prxz>');
            end
            
            matdb(kmat).xten = regexp_dbl(matC,'<xten>(.*)</xten>');
            matdb(kmat).xcmp = regexp_dbl(matC,'<xcmp>(.*)</xcmp>');
            matdb(kmat).yten = regexp_dbl(matC,'<yten>(.*)</yten>');
            matdb(kmat).ycmp = regexp_dbl(matC,'<ycmp>(.*)</ycmp>');
            matdb(kmat).zten = regexp_dbl(matC,'<zten>(.*)</zten>');
            matdb(kmat).zcmp = regexp_dbl(matC,'<zcmp>(.*)</zcmp>');
            matdb(kmat).xy = regexp_dbl(matC,'<xy>(.*)</xy>');
            matdb(kmat).yz = regexp_dbl(matC,'<yz>(.*)</yz>');
            matdb(kmat).xz = regexp_dbl(matC,'<xz>(.*)</xz>');
            matdb(kmat).xycp = regexp_dbl(matC,'<xycp>(.*)</xycp>');
            matdb(kmat).yzcp = regexp_dbl(matC,'<yzcp>(.*)</yzcp>');
            matdb(kmat).xzcp = regexp_dbl(matC,'<xzcp>(.*)</xzcp>');
            matdb(kmat).xzit = regexp_dbl(matC,'<xzit>(.*)</xzit>');
            matdb(kmat).xzic = regexp_dbl(matC,'<xzic>(.*)</xzic>');
            matdb(kmat).yzit = regexp_dbl(matC,'<yzit>(.*)</yzit>');
            matdb(kmat).yzic = regexp_dbl(matC,'<yzic>(.*)</yzic>');
            matdb(kmat).g1g2 = regexp_dbl(matC,'<g1g2>(.*)</g1g2>');
            matdb(kmat).etal = regexp_dbl(matC,'<etal>(.*)</etal>');
            matdb(kmat).etat = regexp_dbl(matC,'<etat>(.*)</etat>');
            matdb(kmat).alp0 = regexp_dbl(matC,'<alp0>(.*)</alp0>');
            
        case 'composite'
            matdb(kmat).thicknessType = regexp_str(matC,'<thicknessType>(.*)</thicknessType>');
            matdb(kmat).uniqueLayers = regexp_dbl(matC,'<uniqueLayers>(.*)</uniqueLayers>');
            matdb(kmat).symmetryType = regexp_str(matC,'<symmetryType>(.*)</symmetryType>');
            
            % find the start and end of each layer block
            layerS = regexp(matC, '<layer>');
            layerE = regexp(matC, '</layer>');
            % adjust start and end indices to encompass only child elements
            layerS = layerS + length('<layer>');
            layerE = layerE - 1;
            
            for klayer = 1:length(layerS)
                layerC = matC(layerS(klayer):layerE(klayer));
                matdb(kmat).layer(klayer).layerName = regexp_str(layerC,'<layerName>(.*)</layerName>');
                matdb(kmat).layer(klayer).thicknessA = regexp_dbl(layerC,'<thicknessA>(.*)</thicknessA>');
                matdb(kmat).layer(klayer).thicknessB = regexp_dbl(layerC,'<thicknessB>(.*)</thicknessB>');
                matdb(kmat).layer(klayer).quantity = regexp_dbl(layerC,'<quantity>(.*)</quantity>',1);
                matdb(kmat).layer(klayer).theta = regexp_dbl(layerC,'<theta>(.*)</theta>');
            end
    end
end
end % end function readMatDB


function str = regexp_str(input_str,pattern,default_str)
    if ~exist('default_str','var')
        default_str = '';
    end
    t = regexp(input_str, pattern, 'tokens');
    if isempty(t)
        str = default_str;
    else
        str = cell2mat(t{1});
    end
end

function dbl = regexp_dbl(input_str,pattern,default_val)
    if ~exist('default_val','var')
        default_val = [];
    end
    t = regexp(input_str, pattern, 'tokens');
    if isempty(t) || isempty(t{1});
        dbl = default_val;
    else
        dbl = str2double(cell2mat(t{1}));
    end
end
