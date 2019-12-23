function x_and_y = create_usable_data(feature_mat, dependent_mat)
    
    if size(feature_mat) ~= size(dependent_mat)
        disp("Error! There is a diasgreement between dependent and independent variables.")
        return
    end
    
    [rows, cols] = size(feature_mat);
    
    nonempty_elements = sum(sum(~isnan(feature_mat))) + sum(sum(~isempty(feature_mat)));
    
    x_and_y = zeros(2, nonempty_elements);
    
    count = 0;
    
    for i = 1:rows
        for j = 1:cols
            if ~isempty(feature_mat(i, j)) && ~isnan(feature_mat(i, j))
                count = count + 1;
                x_and_y(1, count) = feature_mat(i, j);
                x_and_y(2, count) = dependent_mat(i, j);
            end
        end
    end
    
    temp = x_and_y;
    temp(isnan(temp)) = 0;
    mean_replace = mean(temp(2,:));
    
    x_and_y(isnan(x_and_y)) = mean_replace;
    
     

end