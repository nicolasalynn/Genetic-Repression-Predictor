function folding_energies = find_folding_energies(nt_windows, method)
    
    f = waitbar(0, "Calculating Window Energies...");
    temp_size = size(nt_windows);  %1, 3
    if (temp_size(1) == 1 || temp_size(2) ==1)
        disp("successfully entered")
        folding_energies = cell(1, temp_size(2));
        
        for dimension = 1:temp_size(2)
            wait_message = "Looping through series " + num2str(dimension); 
            temp_matrix = zeros(1, length(nt_windows{1, dimension}));
            for i = 1:length(nt_windows{1, dimension})
                waitbar(i/length(nt_windows{1, dimension}), f, wait_message);
                [~, temp_matrix(1, i), ~] = rnafold(nt_windows{1, dimension}(1, i));
            end
            folding_energies{1, dimension} = temp_matrix;
        end
        
    else
     
        [mi_length, ge_length, c] = size(nt_windows);

        folding_energies = zeros(mi_length, ge_length, c);
        
        nt_windows(isempty(nt_windows)) = '';
        nt_windows(ismissing(nt_windows)) = '';
        
        for k = 1:c
            for i = 1:mi_length
                waitbar(i/mi_length, f, "Looping through windows...")

                for j = 1:ge_length
                    waitbar(j*i*k/(mi_length*ge_length*c), f, strcat("Looping through dimesion",num2str(c)));

                    if (nt_windows(i, j) == '')
                        folding_energies(i, j) = NaN;
                    else
                        sequence = char(nt_windows(i, j));
                        [~, folding_energies(i, j), ~] = rnafold(sequence);
                    end

                end
            end
        end
        
    end
    
    if method == "training"
        save('data_sets/feature_data/folding_energies.mat', 'folding_energies')
    elseif method == "validation"
        save('data_sets/validation_data/folding_energies.mat', 'folding_energies')
    end
        close(f)


end