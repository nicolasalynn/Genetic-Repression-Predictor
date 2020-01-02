function X = reshape_with_nan(data, type)

    [x, y, z] = size(data);
    if type == "str"
        X = strings(z, x*y);
        for i = 1:z
            count = 0;
            for j = 1:x
                for k = 1:y
                    count = count + 1;
                    if ~ismissing(data(j, k, i)) || ~isempty(data(j, k, i))
                        X(i, count) = data(j, k, i);
                    else
                        X(i, count) = '';
                    end
                    
                end
            end
        end
        
    elseif type == "num"
        X = zeros(z, x*y);
        for i = 1:z
            count = 0;
            for j = 1:x
                for k = 1:y
                    count = count + 1;
                    if ~isnan(data(j, k, i)) || ~isempty(data(j, k, i))
                        X(i, count) = data(j, k, i);
                    else
                        X(i, count) = 0;
                    end
                    
                end
            end
        end

    else
        disp("error")
    end

end