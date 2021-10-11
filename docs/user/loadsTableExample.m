output=layupDesign_FASTanalysis(blade,DLCoptions,runFASTAnal,useParallel) 
loadsTable =FastLoads4ansys(output,fast_gage,r)
ad2ansys(maptype,nlistfile, loadsTable{1},outfile)
