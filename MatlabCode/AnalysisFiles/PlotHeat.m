function PlotHeat(corrCoef,string)
    figure;
    yvalues = {string,'ForecastWaveSize','ForecastWaveFreq','abs(\gamma_{wave})',...
        'WindSurgeSpeed', 'CurrentSurgeSpeed'};
    xvalues = {string,'ForecastWaveSize','ForecastWaveFreq','abs(\gamma_{wave})',...
        'WindSurgeSpeed', 'CurrentSurgeSpeed'};
    h = heatmap(xvalues,yvalues,corrCoef);
    h.Title = 'Correlation Matrix';
end

