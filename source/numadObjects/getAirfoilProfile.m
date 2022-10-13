function profile = getAirfoilProfile(thickness,percentthick,camber,c)
        %jcb: note that I'm using max thickness about camber
        %instead of overall thickness of airfoil. We may need to
        %change this definition.
    [maxthick,~] = max(thickness);
    tratio = percentthick / (maxthick * 100);
    thick = thickness * tratio;
    hp = camber - 0.5*thick;
    lp = camber + 0.5*thick;

    profile(:,1) = [ c(end); flipud(c) ;  c(2:end);  c(end)];
    profile(:,2) = [ 0     ; flipud(hp); lp(2:end);  0     ];
end