function PlotCosts(best_Fit,it)

% description
% function of plotting cost
    
    plot(it,best_Fit,'r.','MarkerSize',4);
    xlabel('iteration');
    ylabel('Best Fitness value');
    title('optimization result');
    grid on;

end