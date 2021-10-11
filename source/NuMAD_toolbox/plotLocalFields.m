function plotLocalFields(figureNumbers,myYlabel,result,scaleFactor,sym)



    figure(figureNumbers(1))
    m=2;
    n=3;
    %Plot stress in first figure
    subplot(m,n,1)
    plot(result.sig11./scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Stress, $\sigma_{11}$ [MPa]','interpreter','latex','FontWeight','bold')

    subplot(m,n,2)
    plot(result.sig22./scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Stress, $\sigma_{22}$ [MPa]','interpreter','latex','FontWeight','bold')


    subplot(m,n,3)
    plot(result.sig12./scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Stress, $\sigma_{12}$ [MPa]','interpreter','latex','FontWeight','bold')




    %figure(figureNumbers(2))
    
    %Plot strains in second figure
    set(gca,'TickLabelInterpreter','latex')
    subplot(m,n,4)
    hold on
    plot(result.eps11.*scaleFactor,result.x3*1000,sym)
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Strain, $\epsilon_{11}$ [$\mu\epsilon$]','interpreter','latex','FontWeight','bold')

    subplot(m,n,5)
    plot(result.eps22.*scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Strain, $\epsilon_{22}$ [$\mu\epsilon$]','interpreter','latex','FontWeight','bold')


    subplot(m,n,6)
    plot(result.eps12.*scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Strain, $\epsilon_{12}$ [$\mu\epsilon$]','interpreter','latex','FontWeight','bold')

    figure(figureNumbers(2))
    plot(result.eps33.*scaleFactor,result.x3*1000,sym)
    hold on
    ylabel(myYlabel,'interpreter','latex','FontWeight','bold')
    xlabel('Strain, $\epsilon_{33}$ [$\mu\epsilon$]','interpreter','latex','FontWeight','bold')
    

end