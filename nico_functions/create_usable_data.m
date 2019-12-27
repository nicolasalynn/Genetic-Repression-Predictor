function [X, y] = create_usable_data(feature_mat, dependent_mat)
    
    if size(feature_mat) ~= size(dependent_mat)
        disp("Error! There is a diasgreement between dependent and independent variables.")
        return
    end
    
    [rows, cols] = size(feature_mat);
    
    truth1 = feature_mat;
    truth1_1 = ~isnan(truth1);      %matrix of 1 where value is not nan
    truth1_2 = ~isempty(truth1);    %matrix of 1 where value is not empty
    truth1 = truth1_1 + truth1_2;   %matrix where 2 is where it is neither empty or nan, 1 is where it is one, 
    truth1(truth1 < 2) = 0;
    truth1(truth1 == 2) = 1;
    
    truth2 = dependent_mat;         
    truth2_1 = ~isnan(truth2);      
    truth2_2 = ~isempty(truth2);
    truth2 = truth2_1 + truth2_2;
    truth2(truth2 < 2) = 0;
    truth2(truth2 == 2) = 1;
  
    
    truth = truth1 + truth2;        %if there is a both a vlue in the depnent and independent matrix, this will be 2 
    truth(truth < 2) = 0;
    truth(truth == 2) = 1;
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
   

end