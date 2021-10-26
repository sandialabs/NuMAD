function updateFASTFromBLADEDef(params,blade)
    %% Rewright main FAST file and blade file so that properties are consistent with those
    %% in the given blade object.
    fst=readFastMain([params.fstfn '.fst']);
    bld=readFastBlade(fst.BldFile{1}(2:end-1));
    
    fst.TurbConf.TipRad = fst.TurbConf.HubRad + blade.ispan(end);
    bld.prop.PrecrvRef = interp1(blade.ispan,blade.prebend,bld.prop.BlFract.*blade.ispan(end));
    bld.prop.PreswpRef = interp1(blade.ispan,blade.sweep,bld.prop.BlFract.*blade.ispan(end));
    bld.prop.StrcTwst = interp1(blade.ispan,blade.degreestwist,bld.prop.BlFract.*blade.ispan(end));
    
    writeFastBlade(bld,fst.BldFile{1}(2:end-1));
    writeFastMain(fst,[params.fstfn '.fst']);
end

