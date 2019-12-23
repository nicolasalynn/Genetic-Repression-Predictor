function folding_energies = find_folding_energies(nt_windows)
    
    gene_length = size(nt_windows, 1);
    mirna_length = size(nt_windows, 2);
    folding_energies = zeros(gene_length, mirna_length);
    disp(size(nt_windows));
    disp(size(folding_energies));
    for i = 1:gene_length
          
        for j = 1:mirna_length
            
            if (nt_windows(i, j) == '')
                continue
            else
                sequence = char(nt_windows(i, j));
                [discard, folding_energies(i, j), discard2] = rnafold(sequence);
            end
            clearvars discard discard2
        end
    end

    save('folding_energies.mat', 'folding_energies')


end