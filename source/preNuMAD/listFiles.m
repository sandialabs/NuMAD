function filelist = listFiles(directory,filefilter,depth)
%listFiles   List files in directory and subdirectories.
% **********************************************************************
% *                   Part of the SNL NuMAD Toolbox                    *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% *             See license.txt for disclaimer information             *
% **********************************************************************
% filelist = listFiles(DIRECTORY,[FILEFILTER],[DEPTH]);
%   Where:
%     FILEFILTER is an optional wildcard expression which filters the file
%       selection. Default is all files.
%     DEPTH is an optional parameter for number of subdirectories to
%       explore. Default is unlimited subdirectories.
%
%   Usage examples: 
%     filelist = listFiles(DIRECTORY);
%     filelist = listFiles(DIRECTORY,'*.txt');
%     filelist = listFiles(DIRECTORY,[],1);
%
%   See also dir, ls

    if ~exist('filefilter','var') || isempty(filefilter) || isequal(filefilter,'')
        filefilter = '*';
    end
    
    if ~exist('depth','var') || isempty(depth)
        depth = inf;
    end

    dlist = dir(directory);
    filter_regexp = regexptranslate('wildcard',filefilter);
    
    filelist = struct('name','','date','','bytes',0);
    for k = 1:numel(dlist)
        if any(strcmp(dlist(k).name,{'.','..'}))
            continue;  % skip over '.' and '..' directories
        else
            if dlist(k).isdir  % if file is directory
                % recursively process remaining directories, up to "depth"
                % number of sub-directories
                if depth > 0
                    subdir = fullfile(directory,dlist(k).name);
                    sublist = listFiles(subdir,filefilter,depth-1);
                    filelist = [filelist; sublist]; %#ok<AGROW>
                end
            else
                % add file to list if it matches filter
                filtermatch = regexp(dlist(k).name,filter_regexp);
                if filtermatch
                    sublist = struct('name','','date','','bytes',0);
                    sublist.name  = fullfile(directory,dlist(k).name);
                    sublist.date  = dlist(k).date;
                    sublist.bytes = dlist(k).bytes;
                    filelist = [filelist; sublist]; %#ok<AGROW>
                end
            end
        end 
    end
    
    if numel(filelist) >= 1
        % delete the first array position which was initialized empty
        filelist(1) = []; 
    end

end
