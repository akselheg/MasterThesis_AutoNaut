function Plot3Dim(speed, Var1, Var2, lowthresh, highthresh, Flag, YAxis, XAxis, legString)
    table1 = [];table2 = [];table3 = [];
    for i = 1: length(Var2)
        if Var2(i) < lowthresh
            table1 = cat(1,table1,[Var1(i) speed(i)]);
        elseif Var2(i) < highthresh
            table2 = cat(1,table2,[Var1(i) speed(i)]);
        else
            table3 = cat(1,table3,[Var1(i) speed(i)]);
        end
    end
    figure;
    scatter(table1(:,1), table1(:,2))
    hold on 
    scatter(table2(:,1), table2(:,2))
    scatter(table3(:,1), table3(:,2))
    legend(join([legString, ' < ', num2str(lowthresh)]),...
        join([num2str(lowthresh),' < ' ,legString, ' < ', num2str(highthresh)]), ...
        join([legString, ' > ', num2str(highthresh)]))
    xlabel(XAxis);ylabel(YAxis);
    if Flag
        p = polyfit(Var1,speed,1);
        x1 = linspace(min(Var1),max(Var1), length(Var1));
        y1 = polyval(p,x1);
        plot(x1,y1)
        legend(join([legString, ' < ', num2str(lowthresh)]),...
        join([num2str(lowthresh),' < ' ,legString, ' < ', num2str(highthresh)]), ...
        join([legString, ' > ', num2str(highthresh)]), 'Linear Regression')
    end
    hold off
    
end

