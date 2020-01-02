%% Will return the first occurance of particular miRNA seed complement in each gene
%{
    This function will return a matrix with the first occurance of a
    bidning side in each of the 3 segments of code for the gene, all in a
    3D array. Additionally, this function will return the number of
    occurances of binding sites in each region. This infomation may prove
    useful as the number of bidning sites may increase ease of binding and
    increase repression.

    NEEDS TESTING: new modifications include keeping count of the number of
    binding sides in each region as well as generating 2 new dimesions to
    the saved matrix that look for the first binding side in the UTRs.
%}


function binding_indices_validation(mirs_training, gene_training, repress, path, method)

    f = waitbar(0, "Calculating Window Energies...");

    repress_truth = table2array(repress(:, 2:end))';
    repress_truth(~isnan(repress_truth) & ~isempty(repress_truth)) = 1;
    repress_truth(isnan(repress_truth) | isempty(repress_truth)) = 0;
 

% first_indices will keep track of the index of the first binding site in
% EACH of the 3 sequence regions. 

% all_indices will count how many occurances of binding incides occur in
% each of the 3 sequence regions.

% both are 74 rows, 3947 columns, 3 dimensions, where each of the dimasions
% corresponds to UTR5', ORF, UTR3'

    first_indices = zeros(length(mirs_training), size(gene_training, 1), 3); %74 rows, 3947 columns
    all_indices = zeros(length(mirs_training), size(gene_training, 1), 3); 
        
    utr5 = table2array(gene_training(:,2));
    orfs = table2array(gene_training(:, 3));
    utr3 = table2array(gene_training(:,4));
   
    length(mirs_training)
                     
        
    i = 1;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mirna_seq = char(mirs_training(1, i));          
    seed = mirna_seq(2:8);                          
    mer_site_7 = (seed);             
    mer_site_8 = strcat(mer_site_7,'A');             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for j = 1:size(gene_training, 1)            % j = 1:3947

        str_of_utr5 = dna2rna(string(utr5{j}));
        str_of_orf = dna2rna(string(orfs{j}));          
        str_of_utr3 = dna2rna(string(utr3{j}));

        temp_utr5 = regexp(str_of_utr5, mer_site_8);
        temp_orf = regexp(str_of_orf, mer_site_8);            %finding indices in each segment
        temp_utr3 = regexp(str_of_utr3, mer_site_8);



        all_indices(i, j, 1) = length(temp_utr5);
        all_indices(i, j, 2) = length(temp_orf);
        all_indices(i, j, 3) = length(temp_utr3);

        % first_index has a value 0 if there is no binding index, num
        % if there is a singe index and NaN if there are multiple
        % indices

        if isempty(temp_utr5)
            first_indices(i, j, 1) = 0;
        elseif length(temp_utr5) > 1
            first_indices(i, j, 1) = NaN;
        elseif length(temp_utr5) == 1
            first_indices(i, j, 1) = temp_utr5(1);
        else 
            disp("Error!!!")
        end


        if isempty(temp_orf)
            first_indices(i, j, 2) = 0;
        elseif length(temp_orf) > 1
            first_indices(i, j, 2) = NaN;
        elseif length(temp_orf) == 1
            first_indices(i, j, 2) = temp_orf(1); 
        else
            disp("Error!")
        end


        if isempty(temp_utr3)
            first_indices(i, j, 3) = 0;
        elseif length(temp_utr3) > 1
            first_indices(i, j, 3) = NaN;
        elseif length(temp_utr3) == 1
            first_indices(i, j, 3) = temp_utr3(1);
        else
            disp("Error!!!")
        end
            
    end

  
    
    size(all_indices)
    
    index_truths = all_indices;
    index_truths(index_truths ~= 1) = 0;
    
    

    usability(:, :, 1) = index_truths(:, :, 1) + repress_truth;
    usability(:, :, 2) = index_truths(:, :, 2) + repress_truth;
    usability(:, :, 3) = index_truths(:, :, 3) + repress_truth;

    usability(usability ~= 2) = 0;
    usability(usability == 2) = 1;

    % usability is an array that tells you which elements in miRNA x gene
    % tables you can actually use
    
    true_indices = first_indices;
    true_indices(usability ~= 1) = NaN;
    usable_repress(:,:,1) = table2array(repress(:, 2:end))';
    usable_repress(:,:,2) = table2array(repress(:, 2:end))';
    usable_repress(:,:,3) = table2array(repress(:, 2:end))';
    
    usable_repress(usability(1, :, :) ~= 1) = NaN;
    
    
    save(strcat(path, 'true_indices.mat'), 'true_indices') %usable 
    save(strcat(path, 'binary_truth.mat'), 'usability')    %binary table 
    save(strcat(path, 'all_indices.mat'), 'all_indices')   %num of binding sites
    save(strcat(path, 'good_repress.mat'), 'usable_repress')   %usable repress values
    
    
    
    reshaped_indices = reshape_nico(true_indices, "num");
    
    save(strcat(path, 'reshaped_indices.mat'), 'reshaped_indices');
    
    
    close(f)

end

