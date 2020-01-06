function folding_energies = find_folding_energies(nt_windows, method)
    
    f = waitbar(0, "Calculating Window Energies...");
    
    temp_size = size(nt_windows);  %1, 3
    
    if (temp_size(1) == 1 || temp_size(2) ==1)
      
        folding_energies = cell(1, temp_size(2));
        
        for dimension = 1:temp_size(2)
            
            wait_message = "Looping through series " + num2str(dimension); 
            temp_matrix = zeros(1, length(nt_windows{1, dimension}));
            
            for i = 1:length(nt_windows{1, dimension})
                
                waitbar(i/length(nt_windows{1, dimension}), f, wait_message);
                
                sequence = char(nt_windows{1, dimension}(1, i));
                [~, temp_matrix(1, i), ~] = rnafold(sequence);
                
                
            end
            
            folding_energies{1, dimension} = temp_matrix;
        end
    end
    
    if method == "training"
        save('data_sets/feature_data/folding_energies.mat', 'folding_energies')
    elseif method == "validation"
        save('data_sets/validation_data/folding_energies.mat', 'folding_energies')
    end
        close(f)


end