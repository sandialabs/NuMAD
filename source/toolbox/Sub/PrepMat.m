function [Mat,Comp] = PrepMat(data,matdb)
% PREPMAT Filters the full set of material properties and stacks
% contained in MatDBsi.txt down to only the materials and stacks that are
% used by the actual NuMAD model, as described by *.nmd.  fprintf(fid,'lays filter
% results to the screen and the saves final set of materials and stacks to
% Mat and Comp variables.
%
%  Important Notes/Limitations:
%   -?
%

% ========================================================
%   Written by Brian Resor, Sandia National Laboratories
%   Last update: 01/26/2012

% Prep materials information
matids=[];
compids=[];
fid=fopen('PrepMat.txt','wt');

for ii=2:length(data.active.list)
    name=data.active.list{ii};
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

Mat=matdb(matkeepers);
Comp=matdb(compkeepers);

%% check that all matl properties are in place
for i=1:length(Mat)
    % Applicable to the isotropic matls:
    if isempty(Mat(i).ey)
        Mat(i).ey=Mat(i).ex;
    end
    if isempty(Mat(i).gxy)
        Mat(i).gxy=Mat(i).ex/(2*(1+Mat(i).nuxy));
    end
    if isempty(Mat(i).prxy)
        Mat(i).prxy=Mat(i).nuxy;
    end
        
end


end

