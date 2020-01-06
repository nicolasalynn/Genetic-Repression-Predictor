function codon_counts = codon_count(sequences)

    w = waitbar(0, 'Calculating Codon Counts');
    
    [~, dimensions] = size(sequences);
    
    codon_counts = cell(1, 3);
    
    for dim = 1:dimensions

        num_sequences = length(sequences{dim});
        current_sequences = sequences{1, dim};
        waitbar(dim/dimensions, w)
        
        for i = 1:num_sequences
            sequence = current_sequences(1, i);
            codons = struct2cell(codoncount(sequence))';
            total_codons = sum(cell2mat(codons));
            temp = zeros(1, size(codons, 2));
            for j = 1:size(codons, 2)
                temp(j) = (codons{j}/total_codons);
            end

         
            codon_counts{1, dim} = [codon_counts{1, dim}; temp];
        end
        codon_counts{1, dim} = codon_counts{1, dim}';
        
        

    end

    close(w)

end