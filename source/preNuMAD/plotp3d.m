function plotp3d(file)

fid = fopen(file,'rt');
if (fid == -1)
    error('Could not open file "%s"',file);
end

try
    nblocks = fscanf(fid,'%d',1);
    nijk = zeros(nblocks,3);
    for m=1:nblocks
        nijk(m,:) = fscanf(fid,'%d',3);
    end
    
    
    [x, y, z] = deal(cell(1,nblocks));
    for m=1:nblocks
        ni=nijk(m,1);
        nj=nijk(m,2);
        nk=nijk(m,3);
        n_points = ni*nj*nk;
        chunk = fscanf(fid,'%g',n_points);
        x{m} = reshape(chunk,ni,nj,nk);
        chunk = fscanf(fid,'%g',n_points);
        y{m} = reshape(chunk,ni,nj,nk);
        chunk = fscanf(fid,'%g',n_points);
        z{m} = reshape(chunk,ni,nj,nk);
    end
    
    blockcolors = lines(nblocks);
    for m=1:nblocks
        hold on;
        mesh(z{m},x{m},y{m},'edgecolor',blockcolors(m,:));
        hold off;
    end
    fclose(fid);
catch ME
    fclose(fid);
    rethrow(ME);
end

end