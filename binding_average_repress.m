function mean_repression  = binding_average_repress(repression_vals, b_nb, type)

    if (type == "nb")
        indices = (b_nb == '');
    elseif (type == "b")
        indices = (b_nb ~= '');
    else
        disp("You didnt supply enough parameters!")
        exit
    end
        
    temp_values = table2cell(repression_vals);
    temp_values = cell2mat(temp_values(:, 2:end));

    no_vals = isnan(temp_values);

    temp_values(no_vals) = 0;

    temp_values(indices == 0) = 0;
    
    if (type == "nb")
        num_divide = sum(indices == 1) - sum(no_vals);    
    elseif (type == "b")
        num_divide = sum(indices == 1);
    end
    
    mean_repression = sum(temp_values)/num_divide;

end