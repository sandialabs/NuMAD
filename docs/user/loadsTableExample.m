output=layupDesign_FASTanalysis(blade,DLCoptions,runFASTAnal,useParallel) 
loadsTable =FastLoads4ansys(output,fast_gage,r)
beamForceToAnsysShell(maptype,nlistfile, loadsTable{1},outfile)
