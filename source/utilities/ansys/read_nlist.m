function nlist = read_nlist(filename)
if ~exist('filename','var')
    filename = 'NLIST.lis';
end

% Open the file and read the entire contents
fid = fopen(filename);
if (fid == -1)
    error('Could not open file "%s"',filename);
end
filecontents = fread(fid,inf,'uint8=>char')';
fclose(fid);

nlist = cell(1,4);
tbl_hdrs = regexp(filecontents,'NODE\s*X\s*Y\s*Z\s*');
for kTbl = 1:numel(tbl_hdrs)-1
    tbl = filecontents(tbl_hdrs(kTbl):tbl_hdrs(kTbl+1)-1);  % grab table
    data = regexprep(tbl,'\s*NODE\s*X\s*Y\s*Z\s*','');  % remove header
    nlist = [nlist; textscan(data,'%f %f %f %f')];  %#ok (growing data)
end

% get the last table
kTbl=0;  %brr  updated for new make_nlist macro with huges "pages"
tbl = filecontents(tbl_hdrs(kTbl+1):end);  % grab table
data = regexprep(tbl,'\s*NODE\s*X\s*Y\s*Z\s*','');  % remove header
nlist = [nlist; textscan(data,'%f %f %f %f')];

nlist = cell2mat(nlist);
end