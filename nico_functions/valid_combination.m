function [index, repress, mer7_count, ok] = valid_combination(primary, second_one, second_two, mer7, additional, trial)

   if trial == "noutr5"
    
        if length(primary) == 1 && isempty(second_one) && (~isnan(additional))
                index = primary(1);
                repress = additional;
                ok = 1;
                mer7_count = max(0.00001, mer7);
                
        else
                index = NaN;
                repress = NaN;
                ok = 0;
                mer7_count = NaN;
        end
        
   elseif trial == "normal"

            if (length(primary) == 1 && isempty(second_one) && isempty(second_two)) && (~isnan(additional))
                index = primary(1);
                repress = additional;
                ok = 1;
                mer7_count = max(0.00001, mer7);

                
            else
                index = NaN;
                repress = NaN;
                ok = 0;
                mer7_count = NaN;

            end
   end
    
    
end

