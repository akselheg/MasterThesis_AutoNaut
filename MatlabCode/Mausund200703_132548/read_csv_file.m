function out = read_csv_file(file_name)
    out = struct;
    T = readtable([file_name]);
    column_names = T.Properties.VariableNames;
    for i=1:length(column_names)
        out.(column_names{i}) = T{:,column_names{i}};
    end
    
end