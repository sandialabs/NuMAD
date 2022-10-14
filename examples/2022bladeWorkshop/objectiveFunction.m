function objVal = objectiveFunction(XVar,blade,meshData,defLoadsTable,config)
    %% Set the thickness of the spar caps based on the current values of XVar
    blade.components(3).cp(:,2) = XVar(1)*blade.components(3).cp(:,2); % Suction side spar cap
    blade.components(4).cp(:,2) = XVar(2)*blade.components(4).cp(:,2); % Pressure side spar cap
    
    blade.updateBlade();
    
    %% Run ANSYS to determine tip deflection
    ansysResult = mainAnsysAnalysis(blade,meshData,defLoadsTable,config);
    
    flapDef = ansysResult.deflection{1}(end,2);
    
    %% Calculate the value of the objective function
    objVal = (flapDef - 20.0)^2;
end

