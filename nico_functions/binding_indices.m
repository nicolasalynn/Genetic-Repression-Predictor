%% Will return the first occurance of particular miRNA seed complement in each gene

function binding_indices(mirs_training, gene_training, repress, path)

    f = waitbar(0, "Calculating Window Energies...");

    first_indices = zeros(length(mirs_training), size(gene_training, 1), 3); %74 rows, 3947 columns
    valid_repress = zeros(length(mirs_training), size(gene_training, 1), 3);
    num_mer7 = zeros(length(mirs_training), size(gene_training, 1), 3);
    
    utr5 = table2array(gene_training(:,2));
    orfs = table2array(gene_training(:, 3));
    utr3 = table2array(gene_training(:,4));
    
    for i = 1:length(mirs_training)                 
        
        waitbar(i/length(mirs_training), f, "Looping through miRNAs...")
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mirna_seq = char(mirs_training(1, i));          
        seed = mirna_seq(2:8);                          
        mer_site_7 = seqrcomplement(seed);              
        mer_site_8 = rna2dna(strcat(mer_site_7,'A'));            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        for j = 1:size(gene_training, 1)            % j = 1:3947
            
            str_of_utr5 = (string(utr5{j}));
            str_of_orf = (string(orfs{j}));          
            str_of_utr3 = (string(utr3{j}));
            
            temp_utr5 = regexp(str_of_utr5, mer_site_8);
            temp_orf = regexp(str_of_orf, mer_site_8);            
            temp_utr3 = regexp(str_of_utr3, mer_site_8);

            [first_indices(i, j, 1), valid_repress(i, j, 1), num_mer7(i, j, 1), ok1] = valid_combination(temp_utr5, temp_orf, temp_utr3, length(regexp(str_of_utr5, mer_site_7)), repress(i, j), "noutr5");
            [first_indices(i, j, 2), valid_repress(i, j, 2), num_mer7(i, j, 2), ok2] = valid_combination(temp_orf, temp_utr3, temp_utr5, length(regexp(str_of_orf, mer_site_7)), repress(i, j), "noutr5");
            [first_indices(i, j, 3), valid_repress(i, j, 3), num_mer7(i, j, 3), ok3] = valid_combination(temp_utr3, temp_orf, temp_utr5, length(regexp(str_of_utr3, mer_site_7)), repress(i, j), "noutr5");
      
            
            if (ok2 + ok3) > 1
                disp('error')
            end
            
        end

    end

    true_indices = first_indices;
    true_indices(isnan(true_indices)) = 0;
    
    reshaped_repress = reshape_nico(valid_repress, "num");
    reshaped_indices = reshape_nico(first_indices, "num");
    reshaped_mer7 = reshape_nico(num_mer7, "num");
    
    save(strcat(path, 'reshaped_repress.mat'), 'reshaped_repress');
    save(strcat(path, 'reshaped_indices.mat'), 'reshaped_indices');
    save(strcat(path, 'true_indices.mat'), 'true_indices');
    save(strcat(path, 'reshaped_mer7.mat'), 'reshaped_mer7');
    
    close(f)

end
