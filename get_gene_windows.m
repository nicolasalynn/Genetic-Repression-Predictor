%% Returns the 63 nucleotide window around the first index of the seed area

function get_gene_windows(gene_list)

    load('binding_indices.mat')
    orfs = table2array(gene_list(:, 3));
    num_genes = size(site_indices, 1);
    num_mirnas = size(site_indices, 2);
    nt_windows = strings(num_mirnas, num_genes);
    
    for gene = 1:num_genes
        orf = cell2mat(orfs(gene));
        for mirna = 1:num_mirnas
          
            index_val = site_indices(gene, mirna);
            if (index_val == 0)
                continue
            else
                if (index_val <= 26)
                    nt_windows(mirna, gene) = orf(1:63);
                elseif (index_val > strlength(orf) - 35)
                    nt_windows(mirna, gene) = orf(length(orf)-63:length(orf));
                else
                    nt_windows(mirna, gene) = orf(index_val - 26:index_val + 35);
                end
            end
            
        end
    end

    clear binding_indices
    save('nt_windows.mat', 'nt_windows')
end
