function damage = shiftedGoodman(markov,XTEN,XCMP,m,SFs,SFf)

    num_range=size(markov,1)-1;
    num_mean=size(markov,2)-1;
    

    damage=0;
    FSloads = 1.0; % from IEC design standard for DLC 1.2
    for i=1:num_range
        for j=1:num_mean
            sa=markov(i+1,1);
            sm=markov(1,j+1);
            n=markov(i+1,j+1);
            N=((XTEN+abs(XCMP)-abs(2*sm*SFs*FSloads-XTEN+abs(XCMP)))/(2*sa*(SFf)*FSloads))^m;
            damage=damage+n/N;
        end
    end
end