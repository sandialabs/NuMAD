function fst = runIEC_setIC(fst,w, params_parDir)

% ble: add this check to make parallel instance compatible with standard.
if ~exist('params_parDir','var')
    params_parDir = '';
end
% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

% fst = FAST data structure
% w = wind speed at which to find initial conditions

params_parDir = ''; % read from local out file

out=loadFASTOutData([params_parDir 'out\IECSweep_ramp.out']);

fds={'RotSpeed','OoPDefl','IPDefl','TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS'};
chs={'RotSpeed','OoPDefl1','IPDefl1','TeetDefl','Azimuth','NacYaw','TTDspFA','TTDspSS'};
pointerWind=strcmpi(out.list,'WindVxi');

for i = 1:length(fds)
    pointer1=strcmpi(out.list,chs{i});
    pointer2=find(out.data(:,pointerWind)>w, 1 );
    tmp=out.data(pointer2,pointer1);
    if ~isempty(tmp)
        eval(sprintf('fst.Init.%s = %g;',fds{i},tmp));
    end
end



