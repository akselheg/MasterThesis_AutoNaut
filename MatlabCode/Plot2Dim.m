function Plot2Dim(speed,variable, YAxis,xAxis)
    figure;
    scatter(variable, speed)
    hold on 
    p = polyfit(variable, speed, 1);
    x1 = linspace(min(variable), max(variable), length(variable));
    y1 = polyval(p,x1);
    plot(x1,y1)
    xlabel(YAxis),ylabel(xAxis);
    hold off
end

