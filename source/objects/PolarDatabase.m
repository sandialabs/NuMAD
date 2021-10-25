function polardb = PolarDatabase(directory,filter,depth)

    filelist = listFiles(directory,filter,depth);
    
    N = numel(filelist);
    polardb(N) = PolarDef;
    skipped_files = [];
    for k = 1:N
        try
            poldat = PolarDef(filelist(k));
            polardb(k) = poldat;
        catch ME
            id = ME.identifier;
            if isequal(id,'PolarDef:fileNotRecognized')
                fprintf('Skipping file: %s\n',filelist(k).name);
                skipped_files(end+1) = k; %#ok<AGROW>
            else
                rethrow(ME);
            end
        end
    end
    polardb(skipped_files) = [];
    
end