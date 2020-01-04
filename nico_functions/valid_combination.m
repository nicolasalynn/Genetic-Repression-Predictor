function [index, repress] = valid_combination(primary, second_one, second_two, additional)


            if (length(primary) == 1 && isempty(second_one) && isempty(second_two)) && (~isnan(additional))
                index = primary(1);
                repress = additional;
                
            else
                index = NaN;
                repress = NaN;
            end




end