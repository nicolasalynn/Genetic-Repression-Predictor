function [X, y] = create_usable_data(feature_mat, dependent_mat)
    
    if size(feature_mat) ~= size(dependent_mat)
        disp("Error! There is a diasgreement between dependent and independent variables.")
        return
    end
    
    [rows, cols] = size(feature_mat);
    
    truth = feature_mat + dependent_mat;
    truth(~isnan(truth) & ~isempty(truth)) = 1;
    truth(isnan(truth) | isempty(truth)) = 0;

    nonempty_elements = sum(sum(truth));
       
    X = zeros(1, nonempty_elements);
    y = zeros(1, nonempty_elements);
    
    count = 0;

    for i = 1:rows
        for j = 1:cols
            if ~isempty(feature_mat(i, j)) && ~isnan(feature_mat(i, j)) && ~isnan(dependent_mat(i, j)) && ~isempty(dependent_mat(i, j))
                count = count + 1;
                X(count) = feature_mat(i, j);
                y(count) = dependent_mat(i, j);
            end
        end
    end
        
   
    if nonempty_elements ~= count
        disp("error!!! disagreement in non-nan elements")
    end

end