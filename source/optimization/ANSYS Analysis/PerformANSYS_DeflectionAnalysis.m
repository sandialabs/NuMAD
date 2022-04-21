function performANSYS_DeflectionAnalysis(blade, config, iLoad, fid, deflectionFilename)
    fprintf(fid,'/POST1\n');
    fprintf(fid,'set,last\n');
    fprintf(fid,'RSYS,0\n');  %global coordinates
    
    fprintf(fid,'seltol,0.05\n');
    for i=1:numel(blade.ispan)
        fprintf(fid,'*CFOPEN, %s,out\n',[deflectionFilename '-' int2str(i)]);
        fprintf(fid,'ESEL,S,SEC,,1,999   \n');    %Selects aero shell only
        fprintf(fid,'nsle,S,   \n');    %Selects aero shell only
        fprintf(fid,'nsel,r,loc,z,%f  \n',blade.ispan(i));
        %fprintf(fid,'nsll,s,,\n');
        if i==numel(blade.ispan)
            fprintf(fid,'nsel,u,node,,z_master_node_number\n');
        end
        %fprintf(fid,'nplot\n');
        fprintf(fid,'*GET, NsectionNodes, NODE,0,COUNT   !Get the number of nodes in the set\n');
        fprintf(fid,'*GET, node_num, NODE,0,NUM,MIN        !Get the smallest number node in the set\n');
        fprintf(fid,'*DO, i, 1, NsectionNodes                 !loop through all nodes in cross section\n');
        fprintf(fid,'*GET, xpos, NODE,node_num,loc,X\n');
        fprintf(fid,'*GET, ypos, NODE,node_num,loc,Y\n');
        fprintf(fid,'*GET, zpos, NODE,node_num,loc,Z\n');
        fprintf(fid,'*GET, u1, NODE,node_num,U,X\n');
        fprintf(fid,'*GET, u2, NODE,node_num,U,Y\n');
        fprintf(fid,'*GET, u3, NODE,node_num,U,Z\n');
        fprintf(fid,' *VWRITE,node_num,xpos,ypos,zpos,u1,u2,u3\n');
        fprintf(fid,'(E20.12,E20.12,E20.12,E20.12,E20.12,E20.12,E20.12)\n');
        fprintf(fid,'node_num=NDNEXT(node_num)             !Get the next higher node number in the set\n');
        fprintf(fid,'*ENDDO\n');
        fprintf(fid,'*CFCLOS\n');
        fprintf(fid,'\n \n \n');        
    end
    fprintf(fid,'finish\n');
end