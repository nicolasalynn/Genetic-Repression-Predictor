%% Will return the first occurance of particular miRNA seed complement in each gene

function binding_indices(mirs_training, gene_training)
    site_indices = zeros(length(mirs_training), size(gene_training, 1));
    orfs = table2array(gene_training(:, 3));

    for i = 1:length(mirs_training)
        mirna_seq = char(mirs_training(2, i));
        seed = mirna_seq(2:8);
        mer_site = seqrcomplement(seed);
        mer_site = strcat(mer_site,'A');

        for j = 1:size(gene_training, 1)
            str_of_int = dna2rna(string(orfs{j}));
            temp = regexp(str_of_int, mer_site, 'once');
            if isempty(temp)
                site_indices(i, j) = 0;
            else
                site_indices(i, j) = temp(1);
            end
        end

    end

    save('binding_indices.mat', 'site_indices')
end

