clear, clc
%% Genetic Supression Predictor                             
%   Goal: To predict mRNA degradation and supression as a result of miRNA
%   interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%% Training/Validating Initiation                           
clear, clc

addpath nico_functions
addpath lotem_functions
addpath michal_functions
challenge_path = 'data_sets/challenge_data/';
STR = input("Are you training(1) or validating (2)?:");

if STR == 1

    data_path = 'data_sets/feature_data/';
    codon_weights = load(strcat(challenge_path, 'codon_weights.mat')); 
    path = 'data_sets/feature_data/';
    method = "training";
    
    codon_CAI(1,:) = keys(codon_weights.CAI_weights);
    codon_CAI(2,:) = values(codon_weights.CAI_weights);
    codon_tAI(1,:) = keys(codon_weights.tAI_weights);
    codon_tAI(2,:) = values(codon_weights.tAI_weights);

    clearvars codon_weights

    gene_training = load('data_sets/challenge_data/genes_training.mat');
    gene_training = gene_training.genes;


    miRs_training = load('data_sets/challenge_data/miRs_training.mat');
    mirs_training(1,:) = keys(miRs_training.miRs);
    mirs_training(2,:) = values(miRs_training.miRs);

    clearvars miRs_training

    temp = load('data_sets/challenge_data/repress.mat');
    repress = temp.repress;
    clearvars temp
    repress_use = table2array(repress(:, 2:end))';

    save(strcat(challenge_path, 'codon_CAI.mat'), 'codon_CAI');
    save(strcat(challenge_path, 'codon_tAI.mat'), 'codon_tAI');
    save(strcat(challenge_path, 'gene_training_use.mat'), 'gene_training');
    save(strcat(path, 'repress_use.mat'), 'repress_use');
    save(strcat(challenge_path, 'mirs_training_use.mat'), 'mirs_training');
   
    
elseif STR == 2
        
    data_path = 'data_sets/validation_data/';
    codon_weights = load(strcat(challenge_path, 'codon_weights.mat')); 

    path = 'data_sets/validation_data/';
    method = "validation";
    
    codon_CAI(1,:) = keys(codon_weights.CAI_weights);
    codon_CAI(2,:) = values(codon_weights.CAI_weights);
    codon_tAI(1,:) = keys(codon_weights.tAI_weights);
    codon_tAI(2,:) = values(codon_weights.tAI_weights);

    clearvars codon_weights

    gene_training = load('data_sets/challenge_data/genes_validation.mat');
    gene_training = gene_training.genes;

    miRs_training = load('data_sets/challenge_data/miR_validation.mat');
    mirs_training(1,:) = keys(miRs_training.miRs);
    mirs_training(2,:) = values(miRs_training.miRs);

    clearvars miRs_training

   
    repress_use = ones(1, 5805);

    save(strcat(challenge_path, 'codon_CAI.mat'), 'codon_CAI');
    save(strcat(challenge_path, 'codon_tAI.mat'), 'codon_tAI');
    save(strcat(challenge_path, 'gene_training_use.mat'), 'gene_training');
    save(strcat(path, 'repress_use.mat'), 'repress_use');
    save(strcat(challenge_path, 'mirs_training_use.mat'), 'mirs_training');
    
end

%% INDICED OF BINDING                                       

run_initiation = input("Do you want to recalculate the miRNA-mRNA binding "  +  ...
"indices? This action will take approximatelly 2 minutes... \n([Y] = 1, [N] = 0):  ");

if run_initiation 
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training(2, :), gene_training, repress_use, path)
end
clearvars run_initiation

%% WINDOWS                                                  
run_windows = input("Do you want to recalculate the binding windows? " + ...
    "\n([Y] = 1, [N] = 0):  ");

if run_windows
    load(strcat(path, 'true_indices.mat'));
    fprintf("\nThis might take a minute....\n\n");
    get_gene_windows(gene_training, true_indices, 78, method); %by default, set to 78 (70 + 8)
end

clearvars run_windows

%% Load Data -- RUN THIS IF YOU ARE NOT INITIATING          

addpath nico_functions
addpath lotem_functions
addpath michal_functions
challenge_path = 'data_sets/challenge_data/';
begin = input("Training(1) or Validing(2)\n");
if begin == 1
    method = "training";
elseif begin == 2
    method = "validation";
end

if method == "training"
    path = 'data_sets/feature_data/';
elseif method == "validation"
    path = 'data_sets/validation_data/';
end

load(strcat(path, 'reshaped_repress.mat'))
load(strcat(path, "windows_reshaped.mat"))
load(strcat(path, "reshaped_indices.mat"))

%% Feature: Total Length of Sequence                        
if method == "training"
    load(strcat(path, 'whole_reshaped.mat'))
    load(strcat(path, 'reshaped_repress.mat'))

    total_lengths = cell(1, 3);

    for i = 2:3
        seqs = whole_reshaped{i};
        lengths = zeros(1, length(seqs));
        for j = 1:length(seqs)
            lengths(j) = strlength(seqs(j));
        end
        total_lengths{i} = lengths;
        if method == "training"
            slope_of_sequence_length{i} = data_pipeline(total_lengths{i},reshaped_repress{i}); 
        end
    end
end

%% Feature: Thermodynamics                                  
clc
load(strcat(path, "windows_reshaped.mat"))

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>");
if calc_folding_e == 1
    folding_energies = find_folding_energies(windows_reshaped, method); 
end
clearvars calc_folding_e

load(strcat(path, 'folding_energies.mat'));
load(strcat(path, 'reshaped_repress.mat'));

    if method == "training"

        for dim = 2:3

            if dim == 1
                sequence_dec = "UTR 5'";
            elseif dim == 2
                sequence_dec = "ORF";
            else
                sequence_dec = "UTR 3'";
            end

            fprintf("\n" + strcat("Feature: Folding Energy of Binding Window in ", sequence_dec))
            data_pipeline(folding_energies{1, dim}, reshaped_repress{1, dim});

        end
    end
    
clearvars folding_energies reshaped_repress

%% Feature: Conservation                                    
clc

load(strcat(path, 'conservation.mat'))
load(strcat(path, 'reshaped_repress.mat'))
load(strcat(path, 'whole_conservations_reshaped.mat'))

conservation_ratios = cell(1, 3);
conservation_ratios{1, 1} = conservation{1, 1} ./ whole_conservations_reshaped{1, 1};
conservation_ratios{1, 2} = conservation{1, 2} ./ whole_conservations_reshaped{1, 2};
conservation_ratios{1, 3} = conservation{1, 3} ./ whole_conservations_reshaped{1, 3};

save(strcat(path, 'conservation_ratios.mat'), 'conservation_ratios')

if method == "training"

    titles = ["UTR5", "ORF", "UTR3"];

    fprintf("\nFeature: Average Gene Binding Window Conservation\n")

    for i = 1:length(conservation)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(conservation{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Average Whole Gene Conservation\n")

    for i = 1:length(whole_conservations_reshaped)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(whole_conservations_reshaped{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Ratio of Conservation\n")

    for i = 1:length(conservation_ratios)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(conservation_ratios{1, i}, reshaped_repress{1, i});
    end
end

%% Feature: Distance to terminus                            

clc

load(strcat(path, 'reshaped_repress.mat'));
load(strcat(path, 'reshaped_indices.mat'));
load(strcat(path, 'lengths_reshaped.mat'));

[terminus_distance_one, terminus_distance_two] = distance_edge(reshaped_indices, lengths_reshaped, "training");

titles = ["Distance from End, UTR5", "Distance from End, ORF", "Distance from End, UTR3" ...
    "Distance from Either, UTR5", "Distance from Either, ORF", "Distance from Either, UTR3"];

distance_ratio_one = cell(1, 3);
distance_ratio_two = cell(1, 3);
distance_ratio_three = cell(1, 3);
for i = 1:3
    distance_ratio_three{1, i} = reshaped_indices{1, i}./ lengths_reshaped{1, i};
    distance_ratio_one{1, i} = terminus_distance_one{1, i}./lengths_reshaped{1, i};
    distance_ratio_two{1, i} = terminus_distance_two{1, i}./lengths_reshaped{1, i};
end

if method == "training"
    fprintf("\nFeature: Distance to Begining Terminus\n")
    for i = 2:length(reshaped_indices)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(reshaped_indices{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Distance to Closest Terminus\n")
    for i = 2:length(terminus_distance_one)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(terminus_distance_one{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Distance to Closest Terminus\n")

    for i = 2:length(terminus_distance_two)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(terminus_distance_two{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Ratio of Distance to Begining Terminus\n")
    for i = 2:length(reshaped_indices)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(distance_ratio_three{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: Ratio of Distance End Terminus\n")

    for i = 2:length(terminus_distance_one)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(distance_ratio_one{1, i}, reshaped_repress{1, i});
    end
    fprintf("\nFeature: Ratio of Distance to Closest Terminus\n")

    for i = 2:length(terminus_distance_two)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(distance_ratio_two{1, i}, reshaped_repress{1, i});
    end
end


clearvars i ans titles
save(strcat(path, 'terminus_distance_one.mat'), 'terminus_distance_one')
save(strcat(path, 'terminus_distance_two.mat'), 'terminus_distance_two')
save(strcat(path, 'distance_ratio_one.mat'), 'distance_ratio_one')
save(strcat(path, 'distance_ratio_two.mat'), 'distance_ratio_two')
save(strcat(path, 'distance_ratio_three.mat'), 'distance_ratio_three')

%% Feature: CAI                                             

load(strcat(path, 'reshaped_repress.mat'));
load(strcat(challenge_path, 'codon_CAI.mat'));
load(strcat(path, 'corresponding_orf.mat'));

cai_reshaped = cell(1, 3);
cai_whole_reshaped = cell(1,3);

Sequences_ORF = corresponding_orf{1,2};
cai_reshaped{1,2} = CAI_generator(Sequences_ORF,codon_CAI);


Sequences_UTR5 = corresponding_orf{1,1};
cai_reshaped{1, 1} = CAI_generator(Sequences_UTR5,codon_CAI);


Sequences_UTR3 = corresponding_orf{1,3};
cai_reshaped{1, 3} = CAI_generator(Sequences_UTR3,codon_CAI);

if method == "training"
    titles = ["UTR5", "ORF", "UTR3"];

    fprintf("\nFeature: CAI Score of Binding Window \n")

    for i = 2:length(cai_reshaped)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(cai_reshaped{1, i}, reshaped_repress{1, i});
    end
end

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save(strcat(path, 'cai_reshaped.mat'), 'cai_reshaped');

%% Feature: tAI                                             

clc
load(strcat(path, 'windows_reshaped.mat'));
load(strcat(path, 'reshaped_repress.mat'));
load(strcat(challenge_path, 'codon_tAI.mat'));
load(strcat(path, 'whole_reshaped.mat'));
load(strcat(path, 'corresponding_orf.mat'));

tai_reshaped = cell(1, 3);
tai_whole_reshaped = cell(1,3);

Sequences_ORF = corresponding_orf{1,2};
tai_reshaped{1,2} = CAI_generator(Sequences_ORF,codon_tAI);

Sequences_UTR5 = corresponding_orf{1,1};
tai_reshaped{1, 1} = CAI_generator(Sequences_UTR5,codon_tAI);

Sequences_UTR3 = corresponding_orf{1,3};
tai_reshaped{1, 3} = CAI_generator(Sequences_UTR3,codon_tAI);

if method == "training"

    titles = ["UTR5", "ORF", "UTR3"];

    fprintf("\nFeature: tAI Score of Binding Window \n")

    for i = 2:length(tai_reshaped)
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(tai_reshaped{1, i}, reshaped_repress{1, i});
    end
end
clearvars ans tAI_ORF tAI_UTR3 tAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_tAI 
    save(strcat(path, 'tai_reshaped.mat'), 'tai_reshaped'); 

%% Feature: Codon Count of Windows                          
clc

load(strcat(path,'windows_reshaped.mat'));
load(strcat(path,'reshaped_repress.mat'));

codon_counts = codon_count(windows_reshaped);

if method == "atraining"
    for i = 1
        X_var = cell2mat(codon_counts{i, 2});
        y = reshaped_repress{1, 2};
        data_pipeline(X_var, y);
    end
end

save(strcat(path, 'codon_counts.mat'), 'codon_counts');

%% Feature: Codon Counts of ORF                             

clc

load(strcat(path,'corresponding_orf.mat'));
load(strcat(path,'reshaped_repress.mat'));

codon_counts_orfs = codon_count(corresponding_orf);

if method == "atraining"
    for i = 1
        X_var = cell2mat(codon_counts{i, 2});
        y = reshaped_repress{1, 2};
        data_pipeline(X_var, y);
    end
end

save(strcat(path, 'codon_counts_orf.mat'), 'codon_counts_orfs');

%% Feature: GC content                                      

clc

load(strcat(path, 'windows_reshaped.mat'));
load(strcat(path, 'reshaped_repress.mat'));
load(strcat(path, 'whole_reshaped.mat'));
load(strcat(path, 'corresponding_orf.mat'));
load(strcat(path, 'corresponding_utr5.mat'));
load(strcat(path, 'corresponding_utr3.mat'));

gc_orf = cell(1, 3);
gc_utr3 = cell(1, 3);
gc_utr5 = cell(1, 3);


for i = 1:3
    gc_orf{1, i} = GC_content_generator(corresponding_orf{1,i});
    gc_utr5{1, i} = GC_content_generator(corresponding_utr5{1,i});
    gc_utr3{1, i} =  GC_content_generator(corresponding_utr3{1,i});
end

if method == "training"
    titles = ["UTR5", "ORF", "UTR3"];

    fprintf("\nFeature: GC UTR5 Score \n")

    for i = 1:3
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(gc_utr5{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: GC ORF Score \n")

    for i = 1:3
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(gc_orf{1, i}, reshaped_repress{1, i});
    end

    fprintf("\nFeature: GC UTR3 Score  \n")

    for i = 1:3
        fprintf("\nCurrent Sequence Type: %s", titles(i))
        data_pipeline(gc_utr3{1, i}, reshaped_repress{1, i});
    end
end

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save(strcat(path, 'gc_orf.mat'), 'gc_orf')
save(strcat(path, 'gc_utr3.mat'), 'gc_utr3')
save(strcat(path, 'gc_utr5.mat'), 'gc_utr5')

%% MODELLING STEP                                           

type = input("training (1) or validation?(2)\n");
if type == 1
    method = "training";
    path = 'data_sets/feature_data/';
elseif type == 2
    method = "validation";
    path = 'data_sets/validation_data/';
end

clc



% Repression
load(strcat(path, 'reshaped_repress.mat'))

% Codon Counts
load(strcat(path, 'codon_counts.mat'))
load(strcat(path, 'codon_counts_orf.mat'))

% CAI Scores
load(strcat(path, 'cai_reshaped.mat'))


% GC content
load(strcat(path, 'gc_orf.mat'))
load(strcat(path, 'gc_utr5.mat'))
load(strcat(path, 'gc_utr3.mat'))

% Terminus distances
load(strcat(path, 'terminus_distance_one.mat'))
load(strcat(path, 'terminus_distance_two.mat'))
load(strcat(path, 'distance_ratio_one.mat'))
load(strcat(path, 'distance_ratio_two.mat'))
load(strcat(path, 'distance_ratio_three.mat'))
load(strcat(path, 'lengths_reshaped.mat'))
load(strcat(path, 'reshaped_indices.mat'))

% Energies
load(strcat(path, 'folding_energies.mat'))

%   tAI
load(strcat(path, 'tai_reshaped.mat'))

% Conservation
load(strcat(path, 'conservation.mat'))
load(strcat(path, 'conservation_ratios.mat'))
load(strcat(path, 'whole_conservations_reshaped.mat'))

% mer7
load(strcat(path, 'reshaped_mer7.mat'))

% clengths

load(strcat(path, 'cutr5_lengths.mat'))
load(strcat(path, 'cutr3_lengths.mat'))
load(strcat(path, 'corf_lengths.mat'))

X_all = cell(1, 3);
X_recommended = cell(1,3);
y = cell(1, 3);

for i = 2:3
    
    X_all{i} = [codon_counts_orfs{i}', codon_counts{i}', cai_reshaped{i}', conservation{i}', ...
         conservation_ratios{i}', corf_lengths{i}', ...
         cutr3_lengths{i}', cutr5_lengths{i}', ...
         distance_ratio_one{i}', distance_ratio_three{i}', ...
         distance_ratio_two{i}', folding_energies{i}', gc_orf{i}', ...
         gc_utr3{i}', gc_utr5{i}', lengths_reshaped{i}',...
         reshaped_mer7{i}', reshaped_indices{i}', ...
         tai_reshaped{i}', terminus_distance_one{i}', terminus_distance_two{i}', ...
         whole_conservations_reshaped{i}'];
     
    X_recommended{i} = [folding_energies{i}', corf_lengths{i}', cutr3_lengths{i}',...
        cutr5_lengths{i}', conservation{i}', gc_utr5{i}', gc_orf{i}', gc_utr3{i}',...
        reshaped_mer7{i}', terminus_distance_two{i}', cai_reshaped{i}', tai_reshaped{i}'];
        
    y{i} = reshaped_repress{i}';
    
end




% ORF
if method == "training"
    
    clc
    regress_type = input("\nStepwise(1) or Lasso(2):   ");
    feature_type = input("\nAll Features(1) or Select Features(2):    ");
    
    if feature_type == 1
        X = X_all{2};
    elseif feature_type == 2
        X = X_recommended{2};
    end
    
  
    
    if regress_type == 2
        [B, info] = lasso(X, y{2}, 'CV', 5, 'Alpha', 0.5);
        coef = B(:,1);
        coef0 = info.Intercept(1);
        orf_model = struct('coef', coef, 'coef0', coef0);
        normalized_X = coef' * X';
        y_pred = normalized_X' + coef0;
        save('regression_models/orf_model.mat', 'orf_model')
    elseif regress_type == 1
        model = stepwiselm(X, y{2});
        y_pred = predict(model, X);
    end
    
    correlation = corr(y{2}, y_pred) * 100;
    fprintf("\nPearson Lasso Correlation: %.2f%%\n", correlation)
    correlation = corr(y{2}, y_pred, 'type', 'spearman') * 100;
    fprintf("\nSpearman Lasso Correlation: %.2f%%\n", correlation)


    
    % UTR3
   
    if feature_type == 1
        X = X_all{3};
    elseif feature_type == 2
        X = X_recommended{3};
    end
    
    if regress_type == 2
        [B, info] = lasso(X,y{3}, 'CV', 10, 'Alpha', 0.5);
        coef = B(:,1);
        coef0 = info.Intercept(1);
        utr3_model = struct('coef', coef, 'coef0', coef0);
        normalized_X = coef' * X';
        y_pred = normalized_X' + coef0;
    elseif regress_type == 1
        model = stepwiselm(X, y{3});
        y_pred = predict(model, X);
    end  

    correlation = corr(y_pred, y{3}) * 100;
    fprintf("\nPearson Lasso Correlation: %.2f%%\n", correlation)
    correlation = corr(y_pred, y{3}, 'type', 'spearman') * 100;
    fprintf("\nSpearman Lasso Correlation: %.2f%%\n", correlation)
    save('regression_models/utr3_model.mat', 'utr3_model')

elseif method == "validation"
    load('regression_models/orf_model.mat')
     
    i = 2;
    y = reshaped_repress{i}';
     X_all = [codon_counts_orfs{i}', codon_counts{i}',cai_reshaped{i}', conservation{i}', ...
         conservation_ratios{i}', corf_lengths{i}', ...
         corf_lengths{i}', cutr3_lengths{i}', cutr5_lengths{i}', ...
         distance_ratio_one{i}', distance_ratio_three{i}', ...
         distance_ratio_two{i}', folding_energies{i}', gc_orf{i}', ...
         gc_utr3{i}', gc_utr5{i}', lengths_reshaped{i}',...
         reshaped_mer7{i}', reshaped_indices{i}', ...
         tai_reshaped{i}', terminus_distance_one{i}', terminus_distance_two{i}', whole_conservations_reshaped{i}'];
    
     X = X_all;

    coef = orf_model.coef;
    coef0 = orf_model.coef0;
    validation_orf = coef' * X' + coef0; 
 
    load('data_sets/validation_data/reconstruct_index_two.mat')
    pred_orf = reconstruct_data(validation_orf, reconstruct_index_two);
    save('validation_predictions/pred_orf.mat', 'pred_orf');
    
    clearvars coef coef0 X y X_all
    
    load('regression_models/utr3_model.mat')
    i = 3;
    X_all = [codon_counts_orfs{i}', codon_counts{i}',cai_reshaped{i}', conservation{i}', ...
         conservation_ratios{i}', corf_lengths{i}', ...
         corf_lengths{i}', cutr3_lengths{i}', cutr5_lengths{i}', ...
         distance_ratio_one{i}', distance_ratio_three{i}', ...
         distance_ratio_two{i}', folding_energies{i}', gc_orf{i}', ...
         gc_utr3{i}', gc_utr5{i}', lengths_reshaped{i}',...
         reshaped_mer7{i}', reshaped_indices{i}', ...
         tai_reshaped{i}', terminus_distance_one{i}', terminus_distance_two{i}', whole_conservations_reshaped{i}'];
   
     X = X_all;
     
    coef = utr3_model.coef;
    coef0 = utr3_model.coef0;
    validation_utr3 = coef' * X' + coef0; 
    
    load('data_sets/validation_data/reconstruct_index_three.mat')
    pred_utr3 = reconstruct_data(validation_utr3, reconstruct_index_three);
    save('validation_predictions/pred_utr3.mat', 'pred_utr3');
    
    predictions = pred_orf + pred_utr3;
    save('validation_predictions/predictions.mat', 'predictions');
    
end




%% EXTRA>>>>>
% %% Extra
% 
% % Feature: Length of miRNA and repression (find average repression levels across each of 74 miRNAs)
% clear
% 
% 
% load(strcat(path, 'repress_use.mat'))
% load(strcat(path, 'mirs_training_use.mat'))
% 
% mean_repress_miRNA = nanmean(repress_use, 2)';
% mir_length = zeros(1, length(mirs_training));
% for i = 1:length(mirs_training(2, :))
%     mir_length(i) = strlength(mirs_training(2,i));
% end
% 
% clearvars i mirs_training repress_use
% 
% fprintf("\nFeature: Length of miRNA vs Average Repression in all Genes")
% data_pipeline(mir_length, mean_repress_miRNA);
% save(strcat(path, 'mir_length.mat', 'mir_length'))
% save(strcat(path, 'mean_repress_miRNA.mat', 'mean_repress_miRNA'));
% 

%% Feature: Length of ORF and repression 

clear, clc
 
load('data_sets/challenge_data/repress_use.mat')
load('data_sets/challenge_data/gene_training_use.mat')

mean_repress_gene = nanmean(repress_use);
sequences = table2array(gene_training(:, 2:4))';

seq_lengths = zeros(size(sequences));
titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: Gene Length vs Average Repression Across all miRNAs")
for i = 1:size(seq_lengths, 1)
    for j = 1:size(seq_lengths, 2)
        seq_lengths(i, j) = strlength(sequences(i, j));
    end
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    m_average_length{i} = data_pipeline(seq_lengths(i, :), mean_repress_gene);
end

save('data_sets/feature_data/seq_lengths.mat', 'seq_lengths')
save('data_sets/feature_data/mean_repress_gene.mat', 'mean_repress_gene')
save('regression_models/m_average_length.mat', 'm_average_length')
clearvars gene_training i j