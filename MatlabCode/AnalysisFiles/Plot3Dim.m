function Plot3Dim(data1, data2, data3, lowthresh, highthresh, Flag, string1, string2, string3)
    table1 = [];table2 = [];table3 = [];
    for i = 1: length(data3)
        if data3(i) < lowthresh
            table1 = cat(1,table1,[data2(i) data1(i)]);
        elseif data3(i) < highthresh
            table2 = cat(1,table2,[data2(i) data1(i)]);
        else
            table3 = cat(1,table3,[data2(i) data1(i)]);
        end
    end
    figure;
    scatter(table1(:,1), table1(:,2))
    hold on 
    scatter(table2(:,1), table2(:,2))
    scatter(table3(:,1), table3(:,2))
    legend(join([string1, ' < ', num2str(lowthresh)]),...
        join([num2str(lowthresh),' < ' ,string1, ' < ', num2str(highthresh)]), ...
        join([string1, ' > ', num2str(highthresh)]))
    xlabel(string2);ylabel(string3);
    if Flag
        p = polyfit(data2,data1,1);
        x1 = linspace(min(data2),max(data2), length(data2));
        y1 = polyval(p,x1);
        plot(x1,y1)
        legend(join([string1, ' < ', num2str(lowthresh)]),...
        join([num2str(lowthresh),' < ' ,string1, ' < ', num2str(highthresh)]), ...
        join([string1, ' > ', num2str(highthresh)]), 'Linear Regression')
    end
    hold off
    
end

