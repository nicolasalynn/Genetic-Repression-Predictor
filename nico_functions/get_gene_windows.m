%% Returns the 63 nucleotide window around the first index of the seed area

function get_gene_windows(gene_list, indices)
    
    window = 74;
    orfs = table2array(gene_list(:, 3));
    utr5s = table2array(gene_list(:, 2));    %
    utr3s = table2array(gene_list(:, 4));    %
    
 
    [num_mirnas, num_genes, dim] = size(indices);
    nt_windows = strings(num_mirnas, num_genes);
    true_nt_windows = strings(num_mirnas, num_genes, dim);

        for gene = 1:num_genes


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
                    if strlength(utr5) < window
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
    
    
    
    
    clear binding_indices
    save('data_sets/feature_data/true_nt_windows.mat', 'true_nt_windows')
end
