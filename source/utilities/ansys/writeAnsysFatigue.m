function writeAnsysFatigue(fid,iLoad)
    %%%%%%%%%%%%%%%%%%%Outputs for fatigue analysis in MATLAB%%%%%%%%%%%%%%%%%
    fprintf(fid,'! BEGIN FATIGUE SCRIPT\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'esel,all\n'); 
        
    fprintf(fid,'allsel\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'esel,all\n');
    fprintf(fid,'esel,u,type,,21  \n');
    fprintf(fid,'/POST1\n');
    fprintf(fid,'set,LAST\n');
    fprintf(fid,'RSYS,SOLU\n'); %Result in the element coordinate system

    %%% Element strains and curvatures %%%
    fprintf(fid, 'ALLSEL\n');
    fprintf(fid,'ETABLE, zcent,CENT,Z\n'); 
    fprintf(fid,'ETABLE, eps11,SMISC,9 \n');
    fprintf(fid,'ETABLE, eps22,SMISC,10 \n');
    fprintf(fid,'ETABLE, eps12,SMISC,11 \n');
    fprintf(fid,'ETABLE, kapa11,SMISC,12 \n');
    fprintf(fid,'ETABLE, kapa22,SMISC,13 \n');
    fprintf(fid,'ETABLE, kapa12,SMISC,14 \n');
    fprintf(fid,'ETABLE, gamma13,SMISC,15 \n');
    fprintf(fid,'ETABLE, gamma23,SMISC,16 \n');
    fprintf(fid,'/output,plateStrains-all-%s,txt\n',int2str(iLoad));
    fprintf(fid,'PRETAB,zcent,eps11,eps22,eps12,kapa11,kapa22,kapa12,gamma12,gamma13,gamma23\n');
    fprintf(fid, 'ETABLE,ERAS\n\n');   
    

    fprintf(fid,'finish\n');
    fprintf(fid,'! END FATIGUE OUTPUT SCRIPT\n');
end