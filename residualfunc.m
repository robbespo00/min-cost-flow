function residualfunc(x, y)

    
    y = log10(y);
    scatter(x,y,'blue','filled');
    hold on
    p = polyfit(x,y,1);
    fplot(@(t) p(1)*t+p(2), 'red', 'LineWidth', 1.7);
    xlabel('Number of iterations');
    ylabel('Relative residual (log)');