function writeBPEInput(segedges,precurve)
%WRITEBPEINPUT  Write BPE input file
% **********************************************************************
% *           Part of the SNL Wind Turbine Analysis Toolbox            *
% * Developed by Sandia National Laboratories Wind Energy Technologies *
% **********************************************************************
%   writeBPEInput(segedges,precurve)
%   Writes bpe.txt, which is an input file for the BPE algorithm.
%
%      segedges = span distance of BPE segment edges, meters
%      precurve = precurve offsets for swept blades
%
%      Outputs = Write file 'bpe.txt'
%

%===== CREDITS & CHANGELOG ================================================
% 2011.05.26  brr: Initial creation
% yyyy.mm.dd  initials: description

if length(precurve)~=length(segedges)
    error('Defined precurve values must match the desired segment edge locations');
    return
end

fid=fopen('bpe.txt','wt');   %try to open output_file for Writing in Text mode
if (fid == -1)
    error('Could not open file "%s"\n',output_file);
    return
end

wip(fid,length(segedges),         '%number of segment edges desired')
wipcsv(fid,segedges,      '% spanise coords of segment edges')
wipcsv(fid,precurve,'% spanise pre-curve values at segment edges')
wip(fid,1,'% mass weighting flag (0-no, 1-yes)');
wip(fid,6,'% no of terms in z displ expansion ,  minimum = 6');
wip(fid,0,'% conversion, 0 = global coords, 1 = local coords');

fclose(fid);
end

%==========================================================================
%===== FUNCTION DEFINITIONS ===============================================
%==========================================================================
function wip(fid,param,descrip)
% write input file parameter
if ~any(size(param)) && ~ischar(param)
    % do nothing if param = []
    % note: used ~any(size(param)) rather than isempty(param)
    %       so that unset parameters (size = [1 0]) will still
    %       get through to the following elseif statements
elseif ischar(param)
    fprintf(fid,'%-16s ',param);  %output string
elseif isfloat(param)
    if numel(param)==1
        fprintf(fid,'%-16g ',param);  %output single number
    else
        str = sprintf(' %g',param);        %create list of numbers
        str = sprintf('"%s"',str(2:end));  %quote the list of numbers
        fprintf(fid,'%-16s ',str);          %output the quoted list
    end
end
fprintf(fid,'%s\n',descrip);
end

function wipcsv(fid,param,descrip)
% write input file parameter list as comma separated values
str = '';
for k = 1:length(param)
    str = strcat(str,sprintf('%7.4f,',param(k)));
end
if length(str) > 1
    str(end) = [];
end
fprintf(fid,'%-16s %s\n',str,descrip);
end