%% Will return the first occurance of particular miRNA seed complement in each gene
%{
    This function will return a matrix with the first occurance of a
    bidning side in each of the 3 segments of code for the gene, all in a
    3D array. Additionally, this function will return the number of
    occurances of binding sites in each region. This infomation may prove
    useful as the number of bidning sites may increase ease of binding and
    increase repression.
%}


function binding_indices(mirs_training, gene_training, repress, path)

    f = waitbar(0, "Calculating Window Energies...");

    first_indices = zeros(length(mirs_training), size(gene_training, 1), 3); %74 rows, 3947 columns
    valid_repress = zeros(length(mirs_training), size(gene_training, 1), 3);
    all_indices = zeros(length(mirs_training), size(gene_training, 1), 3); 
        
    utr5 = table2array(gene_training(:,2));
    orfs = table2array(gene_training(:, 3));
    utr3 = table2array(gene_training(:,4));
    
    
    for i = 1:length(mirs_training)                 
        
        waitbar(i/length(mirs_training), f, "Looping through miRNAs...")
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        mirna_seq = char(mirs_training(1, i));          
        seed = mirna_seq(2:8);                          
        mer_site_7 = seqrcomplement(seed);              
        mer_site_8 = strcat(mer_site_7,'A');            
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        for j = 1:size(gene_training, 1)            % j = 1:3947
            
            str_of_utr5 = dna2rna(string(utr5{j}));
            str_of_orf = dna2rna(string(orfs{j}));          
            str_of_utr3 = dna2rna(string(utr3{j}));
            
            temp_utr5 = regexp(str_of_utr5, mer_site_8);
            temp_orf = regexp(str_of_orf, mer_site_8);            
            temp_utr3 = regexp(str_of_utr3, mer_site_8);
            

            all_indices(i, j, 1) = length(temp_utr5);
            all_indices(i, j, 2) = length(temp_orf);
            all_indices(i, j, 3) = length(temp_utr3);

            [first_indices(i, j, 1), valid_repress(i, j, 1)] = valid_combination(temp_utr5, temp_orf, temp_utr3, repress(i, j));
            [first_indices(i, j, 2), valid_repress(i, j, 2)] = valid_combination(temp_orf, temp_utr5, temp_utr3, repress(i, j));
            [first_indices(i, j, 3), valid_repress(i, j, 3)] = valid_combination(temp_utr3, temp_utr5, temp_orf, repress(i, j));
           
        
        end

    end

    reshaped_repress = reshape_nico(valid_repress, "num");
    reshaped_indices = reshape_nico(first_indices, "num");

    save(strcat(path, 'reshaped_repress.mat'), 'reshaped_repress');
    save(strcat(path, 'reshaped_indices.mat'), 'reshaped_indices');
    
    
    close(f)

end
