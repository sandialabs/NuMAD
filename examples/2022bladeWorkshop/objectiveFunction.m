function objVal = objectiveFunction(XVar,blade,defLoadsTable,config)
    %% Set the thickness of the spar caps based on the current value of XVar
    blade.components(3).cp(:,2) = XVar*blade.components(3).cp(:,2); % Suction side spar cap
    blade.components(4).cp(:,2) = XVar*blade.components(4).cp(:,2); % Pressure side spar cap
    
    blade.updateBlade();
    
    includeAdhesive = 0;
    meshData=blade.generateShellModel('ansys',includeAdhesive);
    
    %% Run ANSYS to determine tip deflection
    ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,config);
    
    flapDef = ansysResult.deflection{1}(end,2);
    
    %% Calculate the value of the objective function
    objVal = (flapDef - 20.0)^2;
    
    fid = fopen('objectiveHistory.txt','a');
    fprintf(fid,'X: %f \n',XVar);
    fprintf(fid,'flapDef: %f \n',flapDef);
    fprintf(fid,'objVal: %f \n',objVal);
    fprintf(fid,'\n');
    fclose(fid);
    
    blade.components(3).cp(:,2) = (1/XVar)*blade.components(3).cp(:,2); 
    blade.components(4).cp(:,2) = (1/XVar)*blade.components(4).cp(:,2);
    
    blade.updateBlade();
end