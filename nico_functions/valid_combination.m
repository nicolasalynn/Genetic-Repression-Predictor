function [index, repress, ok] = valid_combination(primary, second_one, second_two, additional, trial)

   if trial == "noutr5"
    
        if length(primary) == 1 && isempty(second_one) && (~isnan(additional))
                index = primary(1);
                repress = additional;
                ok = 1;
                
        else
                index = NaN;
                repress = NaN;
                ok = 0;
        end
        
   elseif trial == "normal"

            if (length(primary) == 1 && isempty(second_one) && isempty(second_two)) && (~isnan(additional))
                index = primary(1);
                repress = additional;
                ok = 1;
                
            else
                index = NaN;
                repress = NaN;
                ok = 0;
            end
   end
    
    
end




end