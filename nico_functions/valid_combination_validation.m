function index = valid_combination_validation(primary, second_one, second_two)


            if (length(primary) == 1 && isempty(second_one) && isempty(second_two))
                index = primary(1);
                
            else
                index = NaN;
            end




end