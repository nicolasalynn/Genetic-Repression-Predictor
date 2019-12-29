function folding_energies = find_folding_energies(nt_windows, dim)
    
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
        mi_length = size(nt_windows, 1);
        ge_length = size(nt_windows, 2);
        folding_energies = zeros(mi_length, ge_length);

        nt_windows(isempty(nt_windows)) = '';
        nt_windows(ismissing(nt_windows)) = '';

       nt_windows(1, 112)

        for i = 1:mi_length
            waitbar(i/mi_length, f, "Looping through windows...")

            for j = 1:ge_length

                if (nt_windows(i, j) == '')
                    folding_energies(i, j) = NaN;
                else
                    sequence = char(nt_windows(i, j));
                    [discard, folding_energies(i, j), discard2] = rnafold(sequence);
                end
                clearvars discard discard2

            end
        end
    end
    
    if dim == 1
            save('data_sets/feature_data/folding_energies_utr5.mat', 'folding_energies')        
    elseif dim == 2
            save('data_sets/feature_data/folding_energies_orf.mat', 'folding_energies')
    elseif dim == 3
            save('data_sets/feature_data/folding_energies_utr3.mat', 'folding_energies')
    else
            save('data_sets/feature_data/folding_energies.mat', 'folding_energies')
    end
        close(f)


end