function CreateANSYS_LocalFieldsResults(blade, config, iLoad, fid)
    %%%%%%%%%%%%%%%%%%%Outputs for fatigue analysis in MATLAB%%%%%%%%%%%%%%%%%
    fprintf(fid,'! BEGIN LOCAL FIELD SCRIPT\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'/prep7\n');
    fprintf(fid,'esel,all\n');
    %% RC Addition
    %High Pressure Spar
    fprintf(fid,'esel,s,sect,,3,39,12\n');
    fprintf(fid,'esel,a,sect,,53,389,14\n');
    
    %Low Pressure Spar
    fprintf(fid,'esel,a,sect,,9,45,12\n');
    fprintf(fid,'esel,a,sect,,59,395,14\n');
    %% RC Addition End
    fprintf(fid,'esel,u,type,,21  \n');
    fprintf(fid,'/POST1\n');
    fprintf(fid,'set,LAST\n');
    fprintf(fid,'RSYS,SOLU\n'); %Result in the element coordinate system

    %%% Element Stress %%%
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

    %fprintf(fid,'ETABLE, zcent,CENT,Z\n'); 
    %fprintf(fid,'ETABLE, N11,SMISC,1 \n');
    %fprintf(fid,'ETABLE, N22,SMISC,2 \n');
    %fprintf(fid,'ETABLE, N12,SMISC,3 \n');
    %fprintf(fid,'ETABLE, M11,SMISC,4 \n');
    %fprintf(fid,'ETABLE, M22,SMISC,5 \n');
    %fprintf(fid,'ETABLE, M12,SMISC,6 \n');
    %fprintf(fid,'ETABLE, Q13,SMISC,7 \n');
    %fprintf(fid,'ETABLE, Q23,SMISC,8 \n');
    %fprintf(fid,'/output,plateExamplePlateForces-all-%s,txt\n',int2str(iLoad));
    %fprintf(fid,'/output,plateForces-all-%s,txt\n',int2str(iLoad));
    %fprintf(fid,'PRETAB,zcent,N11,N22,N12,M11,M22,M12,Q12,Q13,Q23\n');
    %fprintf(fid, 'ETABLE,ERAS\n\n'); 
    fprintf(fid,'finish\n');
end
    