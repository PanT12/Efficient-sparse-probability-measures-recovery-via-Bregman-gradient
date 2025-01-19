function [] = portfolio_plot(path,k, titlename)
    standard = readtable([path,'/port',num2str(k),'mean-variance.txt'],VariableNamingRule = 'preserve');
    efficient = readtable([path,'/port',num2str(k),'mean-variance_efficient.txt'],VariableNamingRule = 'preserve');
    standard = table2array(standard);
    efficient = table2array(efficient);
    
    figure;
    set(gcf, 'PaperSize', [25 25])
    plot(standard(:,2),standard(:,3),'Linewidth', 1.5);hold on;
    plot(efficient(:,2),efficient(:,3),'o','Linewidth', 1.5,"MarkerFaceColor",'red');hold off;
    grid on
    set(gca,'FontSize',20);
    xlabel('Variance of return','Fontname','Times New Roman','FontSize',28)
    ylabel("Mean of return",'Fontname','Times New Roman','FontSize',28)
    legend({'SEF','GEF(sparsity 10)'}, ...
        'location','southeast' ,'Fontname','Times New Roman','FontSize',28)
    % xlim([0e-3,3e-3])
    title(titlename,'FontSize',28)
end