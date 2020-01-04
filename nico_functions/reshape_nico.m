function reshaped_vector = reshape_nico(some_matrix, types)


    size_matrix = size(some_matrix);

    if length(size_matrix) > 2

        reshaped_vector = cell(1, size_matrix(3));


        for i = 1:size_matrix(3)
            count = 0;
            if types == "num"                
                elements = sum(sum(~isnan(some_matrix(:,:,i))));
                new_matrix = zeros(1, elements);
                
                for j = 1:size_matrix(1)
                    for k = 1:size_matrix(2)
                        if ~isnan(some_matrix(j, k, i))
                            count = count + 1;
                            new_matrix(1, count) = some_matrix(j, k, i);

                        end
                    end
                end
                reshaped_vector{i} = new_matrix;
                
                
                
                
            elseif types == "str"
                elements = sum(sum(~ismissing(some_matrix(:,:,i))));
                new_matrix = strings(1, elements);
                
                for j = 1:size_matrix(1)
                    for k = 1:size_matrix(2)
                        if ~ismissing(some_matrix(j, k, i))
                            count = count + 1;
                            new_matrix(1, count) = some_matrix(j, k, i);
                        end
                    end
                end
                reshaped_vector{i} = new_matrix;
            end 
            

        end
    end 
end