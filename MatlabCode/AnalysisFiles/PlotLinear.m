function PlotLinear(data, w, X, string)
    diff1 = [];
    sog_MSE = 0;
    output = [];

    for i = 1:length(data)
        out = w'*X(i, :)';
        output = cat(1, output, out);
    end
    mean_output = mean(output);
    sog_summ1 = 0;
    sog_summ2 = 0;
    sog_summ3 = 0;
    for i = 1:length(data)
        out = w'*X(i, :)';
        sog_summ1 = sog_summ1 + (out-mean_output)*(data(i) - mean(data));
        sog_summ2 = sog_summ2 + (out-mean_output)^2;
        sog_summ3 = sog_summ3 + (data(i) - mean(data))^2;
        diff1 = cat(1,diff1, out - data(i));
        sog_MSE = sog_MSE + (data(i) - out)^2;
    end
    sog_r = sog_summ1/sqrt(sog_summ2*sog_summ3);
    sog_RMSE = sqrt(sog_MSE/length(data));
    figure;
    scatter(linspace(1,1,length(diff1)),diff1)
    hold on
    boxplot(diff1)
    tittel  = join(['Error between guessed ', string,  ' and actual ', string]);
    title(tittel)
    hold off
    %%
    figure;
    scatter(data, output) 
    hold on
    [p] = polyfit(data,output,1);
    x1 = linspace(min(data),max(data), length(data));
    y1 = polyval(p,x1);
    plot(x1,y1)
    plot(x1,x1, 'k--')
    legend('Data',' Fit', 'Y = T', 'Location', 'NorthWest')
    tittel = join(['RMSE = ', num2str(sog_RMSE,4),  ', R = ', num2str(sog_r,4)]);
    xlabel(join(['Actual ', string, ' [m/s]']))
    ylabel(join(['Predicted ', string, ' [m/s]']))
    title(tittel)
end

