function PerformANSYS_LinearBucklingAnalysis(blade, config, iLoad, fid, bucklingFilename)
    fprintf(fid,'! BEGIN BUCKLE MACRO\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'/solu\n');
    fprintf(fid,'irlf,-1\n');
    fprintf(fid,'pstres,on\n');
    fprintf(fid,'antype,buckle\n');
    fprintf(fid,strcat('bucopt,lanb,',int2str(config.ansys.analysisFlags.globalBuckling),',,,RANGE\n'));
    %fprintf(fid,strcat('MXPAND,',int2str(nmodes),',0,0,1\n'), nmodes); % Required for element stress/strain, etc..
    fprintf(fid,'solve\n');
    fprintf(fid,'finish\n');
    fprintf(fid,'/post1\n');
    fprintf(fid,['/output,' bucklingFilename ',out\n']);
    fprintf(fid,'set,list\n');
    fprintf(fid,'/output\n');
    fprintf(fid,'finish\n');
    fprintf(fid,'! END BUCKLE MACRO\n'); 
end
