function [freq] = layupDesign_ANSYSfrequency(config)
    global ansysPath;
    ansys_product = 'ANSYS';
    script_name = 'freqAnalysis.mac';
    script_out = 'output.txt';
    
    if(isfield(config.ansys,'np'))
        np = config.ansys.np;
    else
        np = 1;
    end
    
    % File name base name for ansys analysis files
    if isfield(config.ansys,'analysisFileName') && ~isempty(config.ansys.analysisFileName)
        ansysFilename = [config.ansys.analysisFileName];
    else
        ansysFilename=['FEmodel'];
    end
    
    if(isfield(config.ansys,'nFrequencyModes'))
        nModes = config.ansys.nFrequencyModes;
    else
        nModes = 10;
    end
    
    fid=fopen(script_name,'wt');

    fprintf(fid,'/NERR,,99999999\n');
    fprintf(fid,'resume,master,db\n');
    fprintf(fid,'/FILNAME,''%s'',1\n',ansysFilename);
    
    if(isfield(config.ansys,'rpm'))
        fprintf(fid,'\n');
        fprintf(fid,'\n/SOLU');
        fprintf(fid,'\nz_rpm=%f',config.ansys.rpm);
        fprintf(fid,'\nANTYPE,0');
        fprintf(fid,'\nCGLOC,0,0,-3');
        fprintf(fid,'\nOMEGA,0,z_rpm/60*2*3.14,0');
        fprintf(fid,'\nSOLVE');
        fprintf(fid,'\nFINISH');
    end
    fprintf(fid,'\n/SOLU');
    fprintf(fid,'\nANTYPE,2');
    if(isfield(config.ansys,'rpm'))
        fprintf(fid,'\nPSTRES,ON');
    end
    fprintf(fid,'\nMODOPT,LANB,%i ',nModes); 
    fprintf(fid,'\nEQSLV,SPAR  ');
    fprintf(fid,'\nMXPAND,%i ',nModes); 
    fprintf(fid,'\nLUMPM,0 ');
    fprintf(fid,'\n/STATUS,SOLU');
    fprintf(fid,'\nSOLVE   ');
    fprintf(fid,'\nFINISH  ');
    fprintf(fid,'\n/post1');
    fprintf(fid,'\n/output,results_natFreq,out');
    fprintf(fid,'\nset,list');
    fprintf(fid,'\n/output');
    fprintf(fid,'\nfinish');
    fprintf(fid,'\n');
    
    fclose(fid);
    
    delete 'file.lock';
    ansys_call = sprintf('SET KMP_STACKSIZE=2048k & "%s" -b -p %s -I %s -o %s -np %i',...
            ansysPath,ansys_product,script_name,script_out,np);
    disp(ansys_call)
        
    status = 1;
    while(status ~= 0)
        [status,~] = system(ansys_call);  % the windows system call to run the above ansys command
        if(status==0)
            disp('ANSYS frequency analysis completed')
        else
            fprintf('%s: Waiting for ANSYS frequency analysis...\n',datestr(now))
            pause(3);
        end
    end
    
    freq = readAnsysFreq('results_natFreq.out');
    
end

