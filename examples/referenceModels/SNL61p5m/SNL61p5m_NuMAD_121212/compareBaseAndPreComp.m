geezerplot

bld1=readFASTBlade('NRELOffshrBsline5MW_Blade.dat');
bld2=readFASTBlade('FASTBlade_precomp.dat');

fields=fieldnames(bld1.prop);
logs=[4 5 6 7 8 10 11];
x1=getfield(bld1.prop,fields{1});
x2=getfield(bld2.prop,fields{1});
for j=2:length(fields)
    y1=getfield(bld1.prop,fields{j});
    y2=getfield(bld2.prop,fields{j});
    
    figure(1)
    if find(logs==j)
        semilogy(x1,y1,'k:o',...
            x2,y2,'r-v')
    else
        plot(x1,y1,'k:o',...
            x2,y2,'r-v')
    end
    grid on
    legend(strrep('NRELOffshrBsline5MW_Blade.dat','_','\_'),'SNL 61.5m')
    %     legend(strrep('bpe.dat','_','\_'),'precomp')
    xlabel(fields{1})
    ylabel(fields{j})
    
    if 0
        pause
    else
        fn=sprintf('../figs/fig_%s',fields{j});
        print('-dpng',fn);
        
    end
    
end