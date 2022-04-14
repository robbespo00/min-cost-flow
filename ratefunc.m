function ratefunc(x, y)

    scatter(x,y,'green','filled');

    hold on
    p = polyfit(x,y,1);
    fplot(@(t) p(1)*t+p(2), 'red', 'LineWidth', 1.7);

    xlabel('Number of iterations');
    ylabel('Convergence rate');
    ylim([0.75, 1.02]);
end
