function PlotHeat(corrCoef,string)
    figure;
    yvalues = {string,'Wave_{m}','Wave_{Hz}','|cos(\gamma_{Wave})|',...
        'Wind_{u}', 'Current_{u}'};
    xvalues = {string,'Wave_{m}','Wave_{Hz}','|cos(\gamma_{Wave})|',...
        'Wind_{u}', 'Current_{u}'};
    h = heatmap(xvalues,yvalues,corrCoef);
    h.Title = 'Correlation Matrix';
end

