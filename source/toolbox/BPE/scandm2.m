% reads a single line as a string up to the first "%", then truncates blanks
function str=scandm2(fid)
line =fgetl(fid);
str =[];
str2 =[];
if isstr(line)==1
  match = findstr(line,'%');
  str = line(1:match(1)-1);
  str = deblank(str);

end
