%% Returns the 63 nucleotide window around the first index of the seed area

function get_gene_windows(gene_list, indices, window_width, method)
    
    f = waitbar(0, "Calculating Window Energies...");
    
    conservation_vals = gene_list(:, 5);
    indices(isnan(indices)) = 0;
    
    if method == "training"

        path = 'data_sets/feature_data/';
        
    elseif method == "validation"

        path = 'data_sets/validation_data/';

    end
        
        
    window = window_width;
    orfs = table2array(gene_list(:, 3));
    utr5s = table2array(gene_list(:, 2));    
    utr3s = table2array(gene_list(:, 4));    
    
 
    [num_mirnas, num_genes, dim] = size(indices);
    true_nt_windows = strings(num_mirnas, num_genes, dim);
    total_lengths = zeros(num_mirnas, num_genes, dim);
    whole_sequence = strings(num_mirnas, num_genes, dim);
    average_conservation = zeros(num_mirnas, num_genes, dim);
    whole_conservation = zeros(num_mirnas, num_genes, dim);
    corresponding_orf = strings(num_mirnas, num_genes, dim);
    corresponding_utr3 = strings(num_mirnas, num_genes, dim);
    corresponding_utr5 = strings(num_mirnas, num_genes, dim);
    
    corf_lengths = zeros(num_mirnas, num_genes, dim);
    cutr5_lengths = zeros(num_mirnas, num_genes, dim);
    cutr3_lengths = zeros(num_mirnas, num_genes, dim); 

        for mirna = 1:num_mirnas
            
            for gene = 1:num_genes
                
                waitbar(mirna/num_mirnas, f, "Looping through indices...")


                orf = cell2mat(orfs(gene));
                utr5 = cell2mat(utr5s(gene));   
                utr3 = cell2mat(utr3s(gene));   


            %for mirna = 1:num_mirnas
                
                
                index_val_utr5 = indices(mirna, gene, 1);
                index_val_orf = indices(mirna, gene, 2);
                index_val_utr3 = indices(mirna, gene, 3); 


                utr5_length = strlength(utr5);
                orf_length = strlength(orf);
                utr3_length = strlength(utr3);
                
                conservation_vector = cell2mat(conservation_vals{gene, 1})';
                
                if (index_val_utr5 ~= 0)
                    total_lengths(mirna, gene, 1) = strlength(utr5);
                    whole_sequence(mirna, gene, 1) = utr5;
                    whole_conservation(mirna, gene, 1) = mean(conservation_vector);
                    
                    corresponding_orf(mirna, gene, 1) = orf;
                    corresponding_utr3(mirna, gene, 1) = utr3;
                    corresponding_utr5(mirna, gene, 1) = utr5;
                    corf_lengths(mirna, gene, 1) = strlength(orf);
                    cutr5_lengths(mirna, gene, 1) = strlength(utr5);
                    cutr3_lengths(mirna, gene, 1) = strlength(utr3);
                    
                    if strlength(utr5) <= window
                        true_nt_windows(mirna, gene, 1) = utr5;
                        average_conservation(mirna, gene, 1) = mean(conservation_vector(1:utr5_length));
                    elseif (index_val_utr5 <= window/2)
                        true_nt_windows(mirna, gene, 1) = utr5(1:window);
                        average_conservation(mirna, gene, 1) = mean(conservation_vector(1:utr5_length));
                    elseif (index_val_utr5 > strlength(utr5) - window/2)
                        true_nt_windows(mirna, gene, 1) = utr5(length(utr5)-window:end);
                        average_conservation(mirna, gene, 1) = mean(conservation_vector(utr5_length - window:utr5_length));
                    else
                        true_nt_windows(mirna, gene, 1) = utr5(index_val_utr5 - window/2:index_val_utr5 + window/2);
                        average_conservation(mirna, gene, 1) = mean(conservation_vector(index_val_utr5 - window/2:index_val_utr5 + window/2));
                    end
                else
                    true_nt_windows(mirna, gene, 1) = NaN;
                    whole_sequence(mirna, gene, 1) = NaN;
                    average_conservation(mirna, gene, 1) = NaN;
                    corresponding_orf(mirna, gene, 1) = NaN;
                    corresponding_utr3(mirna, gene, 1) = NaN;
                    corresponding_utr5(mirna, gene, 1) = NaN;
                    corf_lengths(mirna, gene, 1) = NaN;
                    cutr5_lengths(mirna, gene, 1) = NaN;
                    cutr3_lengths(mirna, gene, 1) = NaN;
                end
                
                
                if (index_val_orf ~= 0)
                    total_lengths(mirna, gene, 2) = strlength(orf);
                    whole_sequence(mirna, gene, 2) = orf;
                    whole_conservation(mirna, gene, 2) = mean(conservation_vector);

                    corresponding_orf(mirna, gene, 2) = orf;
                    corresponding_utr3(mirna, gene, 2) = utr3;
                    corresponding_utr5(mirna, gene, 2) = utr5;
                    
                    
                    corf_lengths(mirna, gene, 2) = strlength(orf);
                    cutr5_lengths(mirna, gene, 2) = strlength(utr5);
                    cutr3_lengths(mirna, gene, 2) = strlength(utr3);
                    
                    if strlength(orf) <= window
                        true_nt_windows(mirna, gene, 2) = orf;
                        average_conservation(mirna, gene, 2) = mean(conservation_vector(utr5_length:utr5_length + orf_length));
                    elseif (index_val_orf <= window/2)
                        true_nt_windows(mirna, gene, 2) = orf(1:window);
                        average_conservation(mirna, gene, 2) = mean(conservation_vector(utr5_length:utr5_length + orf_length));
                    elseif (index_val_orf > strlength(orf) - window/2)
                        true_nt_windows(mirna, gene, 2) = orf(length(orf)-window:end);
                        average_conservation(mirna, gene, 2) = mean(conservation_vector(utr5_length + orf_length - window:utr5_length + orf_length));
                    else
                        true_nt_windows(mirna, gene, 2) = orf(index_val_orf - window/2:index_val_orf + window/2);
                        average_conservation(mirna, gene, 2) = mean(conservation_vector(utr5_length + index_val_orf - window/2:utr5_length + index_val_orf + window/2));
                    end
                else
                    true_nt_windows(mirna, gene, 2) = NaN;                    
                    whole_sequence(mirna, gene, 2) = NaN;
                    average_conservation(mirna, gene, 2) = NaN;
                    corresponding_orf(mirna, gene, 2) = NaN;
                    corresponding_utr3(mirna, gene, 2) = NaN;
                    corresponding_utr5(mirna, gene, 2) = NaN;
                    corf_lengths(mirna, gene, 2) = NaN;
                    cutr5_lengths(mirna, gene, 2) = NaN;
                    cutr3_lengths(mirna, gene, 2) = NaN;
                end
                    
                    
                if (index_val_utr3~= 0)
                    total_lengths(mirna, gene, 3) = strlength(utr3);
                    whole_sequence(mirna, gene, 3) = utr3;
                    whole_conservation(mirna, gene, 3) = mean(conservation_vector);

                    
                    corresponding_orf(mirna, gene, 3) = orf;
                    corresponding_utr3(mirna, gene, 3) = utr3;
                    corresponding_utr5(mirna, gene, 3) = utr5;
                    
                    corf_lengths(mirna, gene, 3) = strlength(orf);
                    cutr5_lengths(mirna, gene, 3) = strlength(utr5);
                    cutr3_lengths(mirna, gene, 3) = strlength(utr3); 
                    
                    if strlength(utr3) <= window
                        true_nt_windows(mirna, gene, 3) = utr3;
                        average_conservation(mirna, gene, 3) = mean(conservation_vector(utr5_length + orf_length:utr5_length + orf_length + utr3_length));
                    elseif (index_val_utr3 <= window/2)
                        true_nt_windows(mirna, gene, 3) = utr3(1:window);
                        average_conservation(mirna, gene, 3) = mean(conservation_vector(utr5_length + orf_length:utr5_length + orf_length + utr3_length));
                    elseif (index_val_utr3 > strlength(utr3) - window/2)
                        true_nt_windows(mirna, gene, 3) = utr3(length(utr3)-window:end);
                        average_conservation(mirna, gene, 3) = mean(conservation_vector(utr5_length + orf_length + utr3_length - window:utr5_length + orf_length + utr3_length));                        
                    else
                        true_nt_windows(mirna, gene, 3) = utr3(index_val_utr3 - window/2:index_val_utr3 + window/2);
                        average_conservation(mirna, gene, 3) = mean(conservation_vector(utr5_length + orf_length + index_val_utr3 - window/2:utr5_length + orf_length + index_val_utr3 + window/2));                   
                    end
                else
                    true_nt_windows(mirna, gene, 3) = NaN;
                    whole_sequence(mirna, gene, 3) = NaN;
                    average_conservation(mirna, gene, 3) = NaN;
                    
                    corresponding_orf(mirna, gene, 3) = NaN;
                    corresponding_utr3(mirna, gene, 3) = NaN;
                    corresponding_utr5(mirna, gene, 3) = NaN;
                    
                    corf_lengths(mirna, gene, 3) = NaN;
                    cutr5_lengths(mirna, gene, 3) = NaN;
                    cutr3_lengths(mirna, gene, 3) = NaN; 
                end
                    
                    
                          
            end
        end
    
    whole_conservation(whole_conservation == 0) = NaN;
    whole_conservations_reshaped = reshape_nico(whole_conservation, "num");
    windows_reshaped = reshape_nico(true_nt_windows, "str");
    total_lengths(total_lengths == 0) = NaN;
    lengths_reshaped = reshape_nico(total_lengths, "num");
    whole_reshaped = reshape_nico(whole_sequence, "str");
    %average_conservation(average_conservation == 0) = NaN;
    conservation = reshape_nico(average_conservation, "num");
    corresponding_orf = reshape_nico(corresponding_orf, "str");
    corresponding_utr5 = reshape_nico(corresponding_utr5, "str");
    corresponding_utr3 = reshape_nico(corresponding_utr3 , "str");
    cutr5_lengths = reshape_nico(cutr5_lengths, "num");
    cutr3_lengths = reshape_nico(cutr3_lengths, "num");
    corf_lengths = reshape_nico(corf_lengths, "num");

    clear binding_indices
    
    save(strcat(path, 'windows_reshaped', '.mat'), 'windows_reshaped')
    save(strcat(path, 'lengths_reshaped', '.mat'), 'lengths_reshaped')
    save(strcat(path, 'whole_reshaped', '.mat'), 'whole_reshaped')
    save(strcat(path, 'conservation', '.mat'), 'conservation');
    save(strcat(path, 'whole_conservations_reshaped', '.mat'), 'whole_conservations_reshaped')
    save(strcat(path, 'corresponding_orf','.mat'), 'corresponding_orf')
    save(strcat(path, 'corresponding_utr3.mat'), 'corresponding_utr3')
    save(strcat(path, 'corresponding_utr5.mat'), 'corresponding_utr5')
    save(strcat(path, 'corf_lengths.mat'), 'corf_lengths')
    save(strcat(path, 'cutr5_lengths.mat'), 'cutr5_lengths')
    save(strcat(path, 'cutr3_lengths.mat'), 'cutr3_lengths')

    close(f)

end
