%% Will return the first occurance of particular miRNA seed complement in each gene

function reshaped_indices = binding_indices_validation(mirs_training, gene_training, path)


    f = waitbar(0, "Calculating Window Energies...");

    utr5 = table2array(gene_training(:,2));
    orfs = table2array(gene_training(:, 3))';
    utr3 = table2array(gene_training(:,4));

    first_indices = zeros(length(mirs_training), size(gene_training, 1), 3); %74 rows, 3947 columns
    all_indices = zeros(length(mirs_training), size(gene_training, 1), 3); 
        

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mirna_seq = char(mirs_training(1, 1));          
    seed = mirna_seq(2:8);                        
    mer_site_7 = seqrcomplement(seed);             
    mer_site_8 = strcat(mer_site_7, 'A')             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    i = 1;
    temp_orf = zeros(1, size(gene_training, 2));
    for j = 1:5805 

        str_of_utr5 = (string(utr5{j})); 
        str_of_orf = (string(orfs{j})); 
        str_of_utr3 = (string(utr3{j}));

        temp_utr5 = regexp(str_of_utr5, mer_site_8); 
        temp_orf = regexp(str_of_orf, mer_site_8);            %finding indices ineach segment 
        temp_utr3 = regexp(str_of_utr3, mer_site_8);

        all_indices(i, j, 1) = length(temp_utr5); 
        all_indices(i, j, 2) = length(temp_orf); 
        all_indices(i, j, 3) = length(temp_utr3);

        % first_index has a value 0 if there is no binding index, num %if there is a singe index and NaN if there are multiple % indices

        
        first_indices(i, j, 1) = valid_combination_validation(temp_utr5, temp_orf, temp_utr3);
        first_indices(i, j, 2) = valid_combination_validation(temp_orf, temp_utr5, temp_utr3);
        first_indices(i, j, 3) = valid_combination_validation(temp_utr3, temp_utr5, temp_orf);
            
            
    end

    size(all_indices)
    
    index_truths = all_indices; 
    index_truths(index_truths ~= 1) = 0;
    
   
    % usability is an array that tells you which elements in miRNA x gene
    % tables you can actually use
    
    true_indices = first_indices;
    
        
    
    save(strcat(path, 'true_indices.mat'), 'true_indices') %usable
    save(strcat(path, 'all_indices.mat'), 'all_indices')   %num of binding sites 

    reshaped_indices = reshape_nico(true_indices, "num");
    reshaped_indices_one = reshaped_indices{1, 1};
    reshaped_indices_one = reshaped_indices_one(reshaped_indices_one ~= 0);
    reshaped_indices{1, 1} = reshaped_indices_one;
    reshaped_indices_one = reshaped_indices{1, 2};
    reshaped_indices_one = reshaped_indices_one(reshaped_indices_one ~= 0);
    reshaped_indices{1, 2} = reshaped_indices_one;
    reshaped_indices_one = reshaped_indices{1, 3};
    reshaped_indices_one = reshaped_indices_one(reshaped_indices_one ~= 0);
    reshaped_indices{1, 3} = reshaped_indices_one;

    save(strcat(path, 'reshaped_indices.mat'), 'reshaped_indices');
    
    
   close(f)


end

