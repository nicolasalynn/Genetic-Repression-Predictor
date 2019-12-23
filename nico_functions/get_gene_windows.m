%% Returns the 63 nucleotide window around the first index of the seed area

function get_gene_windows(gene_list)

    load('data_sets/feature_data/binding_indices.mat')
    
    orfs = table2array(gene_list(:, 3));
    utr5s = table2array(gene_list(:, 2));    %
    utr3s = table2array(gene_list(:, 4));    %
    
    num_mirnas = size(first_indices, 1);
    num_genes = size(first_indices, 2);
    nt_windows = strings(num_mirnas, num_genes);
    
    for gene = 1:num_genes
        
        
        orf = cell2mat(orfs(gene));
        utr5 = cell2mat(utr5s(gene));   %
        utr3 = cell2mat(utr3s(gene));   %
        
        
        for mirna = 1:num_mirnas
            
            index_val_orf = first_indices(mirna, gene, 2);
            index_val_utr5 = first_indices(mirna, gene, 1); %
            index_val_utr3 = first_indices(mirna, gene, 3); %
            
            if (index_val_utr5 ~= 0)
                if strlength(utr5) < 63
                    nt_windows(mirna, gene) = utr5;
                elseif (index_val_utr5 <= 26)
                    nt_windows(mirna, gene) = utr5(1:63);
                elseif (index_val_utr5 > strlength(utr5) - 35)
                    nt_windows(mirna, gene) = utr5(length(utr5)-63:length(utr5));
                else
                    nt_windows(mirna, gene) = utr5(index_val_utr5 - 26:index_val_utr5 + 35);
                end
                
            elseif (index_val_utr5 == 0)
                
                if (index_val_orf ~= 0)
                    
                    if (index_val_orf <= 26)
                        nt_windows(mirna, gene) = orf(1:63);
                    elseif (index_val_orf > strlength(orf) - 35)
                        nt_windows(mirna, gene) = orf(length(orf)-63:length(orf));
                    else
                        nt_windows(mirna, gene) = orf(index_val_orf - 26:index_val_orf + 35);
                    end
                    
                elseif (index_val_orf == 0)
                    
                    if (index_val_utr3 ~= 0)
                        if strlength(utr3) < 63
                            nt_windows(mirna, gene) = utr3;
                        elseif (index_val_utr3 <= 26)
                            nt_windows(mirna, gene) = utr3(1:63);
                        elseif (index_val_utr3 > strlength(utr3) - 35)
                            nt_windows(mirna, gene) = utr3(length(utr3)-63:length(utr3));
                        else
                            nt_windows(mirna, gene) = utr3(index_val_utr3 - 26:index_val_utr3 + 35);
                    
                        end
                        
                        
                    elseif (index_val_utr3 == 0)
                        continue
                        
                    end
                end
            end
            
        end
    end

    clear binding_indices
    save('data_sets/feature_data/nt_windows.mat', 'nt_windows')
end
