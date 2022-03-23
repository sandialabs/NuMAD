function PerformANSYS_ResultantVSSpanAnalysis(blade, config, iLoad, fid)
    fprintf(fid,'/POST1\n');
    fprintf(fid,'set,LAST\n');
    fprintf(fid,'RSYS,0\n'); %global coordinates

    %fprintf(fid,'seltol,0.05\n');
    %fprintf(fid,'*CFOPEN, resultantVSspan,txt\n');
    %for i=1:numel(blade.ispan)
        %fprintf(fid,'nsel,s,loc,z,0,%f  \n',blade.ispan(i));
 
        %if i==numel(blade.ispan)
            %fprintf(fid,'nsel,u,node,,z_master_node_number\n');
        %end
 
        %fprintf(fid,'spoint,0,%f,%f,%f\n',blade.sweep(i),blade.prebend(i),blade.ispan(i));
        %fprintf(fid,'nplot\n');
        %fprintf(fid,'FSUM\n');
        %fprintf(fid,'*GET, F1, FSUM, 0, ITEM,FX\n');
        %fprintf(fid,'*GET, F2, FSUM, 0, ITEM,FY\n');
        %fprintf(fid,'*GET, F3, FSUM, 0, ITEM,FZ\n');
        %fprintf(fid,'*GET, M1, FSUM, 0, ITEM,MX\n');
        %fprintf(fid,'*GET, M2, FSUM, 0, ITEM,MY\n');
        %fprintf(fid,'*GET, M3, FSUM, 0, ITEM,MZ\n');
        %fprintf(fid,'*VWRITE,%f,F1,F2,F3,M1,M2,M3\n',blade.ispan(i));
        %fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
        %fprintf(fid,'\n \n \n');
    %end
    %fprintf(fid,'*CFCLOS\n');
    %fprintf(fid,'finish\n');            
            
    fprintf(fid,'/post1\n');
    fprintf(fid,'elsize=%f\n',blade.mesh);
    fprintf(fid,'nz=nint(%f/elsize) !Integer number of points to output resultant loads\n',blade.ispan(end));
    fprintf(fid,'zloc=0\n');
    fprintf(fid,'delta=0.1\n');
    fprintf(fid,'*CFOPEN, resultantVSspan,txt\n');
    fprintf(fid,'*do,I,1,nz+1\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'nsel,s,loc,z,0,zloc+delta\n');
    fprintf(fid,'spoint,0,0,0,zloc\n');
    fprintf(fid,'!nplot\n');
    fprintf(fid,'FSUM\n');
    fprintf(fid,'*GET, F1, FSUM, 0, ITEM,FX\n');
    fprintf(fid,'*GET, F2, FSUM, 0, ITEM,FY\n');
    fprintf(fid,'*GET, F3, FSUM, 0, ITEM,FZ\n');
    fprintf(fid,'*GET, M1, FSUM, 0, ITEM,MX\n');
    fprintf(fid,'*GET, M2, FSUM, 0, ITEM,MY\n');
    fprintf(fid,'*GET, M3, FSUM, 0, ITEM,MZ\n');
    fprintf(fid,'*VWRITE,zloc,F1,F2,F3,M1,M2,M3\n');
    fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
    fprintf(fid,'zloc=zloc+elsize\n');
    fprintf(fid,'*ENDDO\n');
    fprintf(fid,'*CFCLOS\n');
    fprintf(fid,'finish\n');
end