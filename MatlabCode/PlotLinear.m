function PlotLinear(test_sog_data, w1, X_test, string)
diff1 = [];
sog_MSE = 0;
output = [];

for i = 1:length(test_sog_data)
    out = w1'*X_test(i, :)';
    output = cat(1, output, out);
end
mean_output = mean(output);
sog_summ1 = 0;
sog_summ2 = 0;
sog_summ3 = 0;
for i = 1:length(test_sog_data)
    out = w1'*X_test(i, :)';
    sog_summ1 = sog_summ1 + (out-mean_output)*(test_sog_data(i) - mean(test_sog_data));
    sog_summ2 = sog_summ2 + (out-mean_output)^2;
    sog_summ3 = sog_summ3 + (test_sog_data(i) - mean(test_sog_data))^2;
    diff1 = cat(1,diff1, out - test_sog_data(i));
    sog_MSE = sog_MSE + (test_sog_data(i) - out)^2;
end
sog_r = sog_summ1/sqrt(sog_summ2*sog_summ3);
sog_RMSE = sqrt(sog_MSE/length(test_sog_data));
figure;
scatter(linspace(1,1,length(diff1)),diff1)
hold on
boxplot(diff1)
tittel  = join(['Error between guessed', string,  ' and actual ', string]);
title(tittel)
hold off
%%
figure;
scatter(test_sog_data, output) 
hold on
[p] = polyfit(test_sog_data,output,1);
x1 = linspace(min(test_sog_data),max(test_sog_data), length(test_sog_data));
y1 = polyval(p,x1);
plot(x1,y1)
plot(x1,x1, 'k--')
legend('Data',' Fit', 'Y = T', 'Location', 'NorthWest')
string = join(['RMSE = ', num2str(sog_RMSE,4),  ', R = ', num2str(sog_r,4)]);
xlabel 'Actual Vg [m/s]'
ylabel 'Predicted Vg [m/s]'
title(string)
end

