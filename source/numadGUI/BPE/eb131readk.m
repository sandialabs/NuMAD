function eb131readk(filename)
%  eb26readk.m
%  reads in pairs of specified k matrices
%  latest update:  17 march 03
global givenk ngivenk npairk  zgivenk zpairk
%  'start of eb26readk'

%  givenk = zeros((npairk-1)*16-2,6);

givenk = load(filename);
% eval(['load ' filename]);                % loads all file into temp
% [nr,nc] = size(filename);
% newname = trimDM(filename);
% [nr,nc] = size(newname);		          % assumes that string end with single quote sign
% name = newname(:,1:nc-5);
% givenk = eval(name);
%
[nrgivenk,ncgivenk] = size(givenk);


for ipairk = 1:npairk
    zpairk(ipairk,1) = zgivenk((ipairk-1)*2+1);      % store start and stop z values as pairs
    zpairk(ipairk,2) = zgivenk((ipairk-1)*2+2);
end
