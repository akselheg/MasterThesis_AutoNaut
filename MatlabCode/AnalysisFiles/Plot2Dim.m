function Plot2Dim(data1,data2, strXLabel, strYLabel)
    figure;
    scatter(data2, data1)
    hold on 
    p = polyfit(data2, data1, 1);
    x1 = linspace(min(data2), max(data2), length(data2));
    y1 = polyval(p,x1);
    plot(x1,y1)
    xlabel(strXLabel),ylabel(strYLabel);
    hold off
end

