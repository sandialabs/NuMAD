function loadNuMAD_matl(arg)
% LOADNUMAD_MATL  Convert MatDBsi.txt into Matlab data structure.
%
%  Important Notes/Limitations:
%   -?
%

% ========================================================
%   Written by Brian Resor, Sandia National Laboratories
%   Last update: 01/21/2011

fid=fopen('MatDBsi.txt');

mat=0;
comp=0;
flag=0;
layer=0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    
    if ~isempty(strfind(tline,'<material>'))
        % do nothing
    elseif ~isempty(strfind(tline,'</material>'))
        flag=0;
        layer=0;
    elseif ~isempty(strfind(tline,'<type>isotropic</type>')) || ~isempty(strfind(tline,'<type>orthotropic</type>'))
        mat=mat+1;
        flag=1;
    elseif ~isempty(strfind(tline,'<type>composite</type>'))
        comp=comp+1;
        flag=2;
    end
    
    switch flag
        case 1  % isotropic or orthotropic
            pat='<(?<name>\w+)>(?<setting>.+)</';
            n=regexp(tline,pat,'names');
            if ~isempty(n)
                n.setting=strrep(n.setting,'''','''''');
                try
                    eval(['Mat(' num2str(mat) ').' n.name '=' n.setting ';']);
                catch
                    eval(['Mat(' num2str(mat) ').' n.name '=''' n.setting ''';']);
                end
            end
        case 2 % composite material
            if ~isempty(strfind(tline,'<layer>'))
                layer=layer+1;
            elseif ~isempty(strfind(tline,'</layer>'))
                %do nothing
            else
                pat='<(?<name>\w+)>(?<setting>.+)</';
                n=regexp(tline,pat,'names');
                if ~isempty(n)
                    n.setting=strrep(n.setting,'''','''''');
                    if layer==0 % enter material information
                        try
                            eval(['Comp(' num2str(comp) ').' n.name '=' n.setting ';']);
                        catch
                            eval(['Comp(' num2str(comp) ').' n.name '=''' n.setting ''';']);
                        end
                    else % enter another layer
                        try
                            eval(['Comp(' num2str(comp) ').layer(' num2str(layer) ').' n.name '=' n.setting ';']);
                        catch
                            eval(['Comp(' num2str(comp) ').layer(' num2str(layer) ').' n.name '=''' n.setting ''';']);
                        end
                        
                    end
                end
            end
        otherwise
            % do nothing
    end
end


fclose(fid);

save SDataMat Mat Comp

end