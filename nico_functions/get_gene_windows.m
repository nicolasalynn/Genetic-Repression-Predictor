%% Returns the 63 nucleotide window around the first index of the seed area

function get_gene_windows(gene_list, indices, file_save_name, window_width)
    
    f = waitbar(0, "Calculating Window Energies...");

    indices(isnan(indices)) = 0;
    
    file_name_1 = strcat('data_sets/feature_data/', char(file_save_name), '.mat');
    file_name_2 = strcat('data_sets/feature_data/reshaped_', char(file_save_name), '.mat');

    
    window = window_width;
    orfs = table2array(gene_list(:, 3));
    utr5s = table2array(gene_list(:, 2));    %
    utr3s = table2array(gene_list(:, 4));    %
    
 
    [num_mirnas, num_genes, dim] = size(indices);
    true_nt_windows = strings(num_mirnas, num_genes, dim);

        for gene = 1:num_genes
            waitbar(gene/num_genes, f, "Looping through indices...")


            orf = cell2mat(orfs(gene));
            utr5 = cell2mat(utr5s(gene));   
            utr3 = cell2mat(utr3s(gene));   


            for mirna = 1:num_mirnas

                index_val_utr5 = indices(mirna, gene, 1);
                index_val_orf = indices(mirna, gene, 2);
                index_val_utr3 = indices(mirna, gene, 3); 


                if (index_val_utr5 ~= 0)
                    if strlength(utr5) < window
                        true_nt_windows(mirna, gene, 1) = utr5;
                    elseif (index_val_utr5 <= window/2)
                        true_nt_windows(mirna, gene, 1) = utr5(1:window);
                    elseif (index_val_utr5 > strlength(utr5) - window/2)
                        true_nt_windows(mirna, gene, 1) = utr5(length(utr5)-window:end);
                    else
                        true_nt_windows(mirna, gene, 1) = utr5(index_val_utr5 - window/2:index_val_utr5 + window/2);
                    end
                else
                    true_nt_windows(mirna, gene, 1) = NaN;
                end
                
                
                if (index_val_orf ~= 0)
                    if strlength(orf) < window
                        true_nt_windows(mirna, gene, 2) = orf;
                    elseif (index_val_orf <= window/2)
                        true_nt_windows(mirna, gene, 2) = orf(1:window);
                    elseif (index_val_orf > strlength(orf) - window/2)
                        true_nt_windows(mirna, gene, 2) = orf(length(orf)-window:end);
                    else
                        true_nt_windows(mirna, gene, 2) = orf(index_val_orf - window/2:index_val_orf + window/2);
                    end
                else
                    true_nt_windows(mirna, gene, 2) = NaN;
                end
                    
                    
                if (index_val_utr3~= 0)
                    if strlength(utr3) < window
                        true_nt_windows(mirna, gene, 3) = utr3;
                    elseif (index_val_utr3 <= window/2)
                        true_nt_windows(mirna, gene, 3) = utr3(1:window);
                    elseif (index_val_utr3 > strlength(utr3) - window/2)
                        true_nt_windows(mirna, gene, 3) = utr3(length(utr3)-window:end);
                    else
                        true_nt_windows(mirna, gene, 3) = utr3(index_val_utr3 - window/2:index_val_utr3 + window/2);
                    end
                else
                    true_nt_windows(mirna, gene, 3) = NaN;
                end
                    
                    
                          
            end
        end
    
    
    windows_reshaped = reshape_nico(true_nt_windows, "str");
    
    clear binding_indices
    save(file_name_1, 'true_nt_windows')
    save(file_name_2, 'windows_reshaped')
    
    close(f)

end
