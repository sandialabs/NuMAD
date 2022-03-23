function PerformANSYS_FailureAnalysis(blade, config, iLoad, fid, failureFilename)
    fprintf(fid,'! BEGIN FAILURE SCRIPT\n');
    fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    fprintf(fid,'!Add for PLESOL and *get,findex,PLNSOL,0,MAX to work');
    fprintf(fid,'/BATCH  \n');
    fprintf(fid,'/COM,ANSYS RELEASE Release 18.1      BUILD 18.1      UP20170403       15:49:08\n');
    fprintf(fid,'/GRA,POWER\n ');
    fprintf(fid,'/GST,ON\n ');
    fprintf(fid,'/PLO,INFO,3\n ');
    fprintf(fid,'/GRO,CURL,ON\n ');
    fprintf(fid,'/CPLANE,1   \n ');
    fprintf(fid,'/REPLOT,RESIZE  \n ');
    fprintf(fid,'WPSTYLE,,,,,,,,0\n ');
    fprintf(fid,'/SHOW\n ');
    fprintf(fid,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n');
    fprintf(fid,'/POST1\n');
    fprintf(fid,'set,last\n');
    fprintf(fid,'allsel\n');
    fprintf(fid,'esel,all\n');
    fprintf(fid,'RSYS,LSYS \n'); %Result in the layer coordinate system
    fprintf(fid,'layer,fcmax\n');

    if ~contains(upper(config.ansys.analysisFlags.failure),{'PUCK' 'LARC03' 'LARC04'})
        %Do this for the failure criteria that do not distinguish between fiber and matrix failure
        fc=upper(config.ansys.analysisFlags.failure);
        fprintf(fid,'FCTYP,add,%s\n',fc);
        fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
        fprintf(fid,'*get,findex,PLNSOL,0,MAX\n');
    else
        %Do this for the failure criteria that do distinguish between fiber and matrix failure
        switch upper(config.ansys.analysisFlags.failure)
            case 'PUCK'
                %Fiber Failure
                fc='PFIB';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                %Matrix Failure
                fc='PMAT';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
            case 'LARC03'
                %Fiber Failure
                fc='L3FB';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                %Matrix Failure
                fc='L3MT';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
            case 'LARC04'
                %Fiber Failure
                fc='L4FB';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Ffindex,PLNSOL,0,MAX\n');

                %Matrix Failure
                fc='L4MT';
                fprintf(fid,'FCTYP,add,%s\n',fc);
                fprintf(fid,'PLESOL, FAIL,%s, 0,1.0\n',fc);
                fprintf(fid,'*get,Mfindex,PLNSOL,0,MAX\n');
       end
       %Report the higher of the fiber failure index or the matrix
       fprintf(fid,'*IF, Ffindex, GT,Mfindex, THEN\n');
       fprintf(fid,'findex=Ffindex\n');
       fprintf(fid,'*ELSE\n');
       fprintf(fid,'findex=Mfindex\n');
       fprintf(fid,'*ENDIF\n');
    end
            
    fprintf(fid,['/output,' failureFilename ',out\n']);
    fprintf(fid,'*status,findex\n');
    fprintf(fid,'/output\n');
    fprintf(fid,'finish\n');
    fprintf(fid,'! END FAILURE SCRIPT\n');
end
    