clear, clc
%% Genetic Supression Predictor -- TRAINING -- RUN ME
%   Goal: To predict mRNA degradation and supression as a result of miRNA
%   interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%% Training Initiation
clear, clc

addpath nico_functions
addpath lotem_functions
addpath michal_functions

STR = input("Are you training(1) or validating (2)?:");

if STR == 1

    data_path = 'data_sets/feature_data/';
    challenge_path = 'data_sets/challenge_data/';
    codon_weights = load('data_sets/challenge_data/codon_weights.mat'); 

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
    save(strcat(challenge_path, 'repress_use.mat'), 'repress_use');
    save(strcat(challenge_path, 'mirs_training_use.mat'), 'mirs_training');
    
    path = 'data_sets/feature_data/';
    method = "training";
    
elseif STR == 2
    
    
end

%% Find the first instance of miRNA mRNA binding for each combination (Nico)


run_initiation = input("Do you want to recalculate the miRNA-mRNA binding "  +  ...
"indices? This action will take approximatelly 2 minutes... \n([Y] = 1, [N] = 0):  ");

if run_initiation 
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training(2, :), gene_training, repress_use, path)
end
clearvars run_initiation

%% Obtain windows of specified length for all indices found previously
run_windows = input("Do you want to recalculate the binding windows? " + ...
    "\n([Y] = 1, [N] = 0):  ");

if run_windows
    load('data_sets/feature_data/true_indices.mat');
    fprintf("\nThis might take a minute....\n\n");
    get_gene_windows(gene_training, true_indices, 78, method); %by default, set to 78 (70 + 8)
end

clearvars run_windows

%% Load Data -- RUN THIS IF YOU ARE NOT INITIATING

addpath nico_functions
addpath lotem_functions
addpath michal_functions

if method == "training"
    path = 'data_sets/feature_data/';
elseif method == "validation"
    path = 'data_sets/validation_data/';
end

load(strcat(path, 'reshaped_repress.mat'))
load(strcat(path, "windows_reshaped.mat"))
load(strcat(path, "reshaped_indices.mat"))

%% Feature: Total Length of Sequence
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
    slope_of_sequence_length{i} = data_pipeline(total_lengths{i},reshaped_repress{i}); 

end

%% Feature: Thermodynamics
clc
load(strcat(path, "windows_reshaped.mat"))

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>");
if calc_folding_e == 1
    folding_energies = find_folding_energies(windows_reshaped, "training"); 
end
clearvars calc_folding_e

load(strcat(path, 'folding_energies.mat'));
load(strcat(path, 'reshaped_repress.mat'));

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

save(strcat(path, 'conservation_ratios.mat', 'conservation_ratios'))

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
    distance_ratio_three{1, i} = reshaped_indices{1, i}./ lengths_reshaped{1, i};
    data_pipeline(distance_ratio_three{1, i}, reshaped_repress{1, i});
end

fprintf("\nFeature: Ratio of Distance End Terminus\n")

for i = 2:length(terminus_distance_one)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    distance_ratio_one{1, i} = terminus_distance_one{1, i}./lengths_reshaped{1, i};
    data_pipeline(distance_ratio_one{1, i}, reshaped_repress{1, i});
end


fprintf("\nFeature: Ratio of Distance to Closest Terminus\n")

for i = 2:length(terminus_distance_two)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    distance_ratio_two{1, i} = terminus_distance_two{1, i}./lengths_reshaped{1, i};
    data_pipeline(distance_ratio_two{1, i}, reshaped_repress{1, i});
end

clearvars i ans titles
save(strcat(path, 'terminus_distance_one.mat', 'terminus_distance_one'))
save(strcat(path, 'terminus_distance_two.mat', 'terminus_distance_two'))
save(strcat(path, 'distance_ratio_one.mat', 'distance_ratio_one'))
save(strcat(path, 'distance_ratio_two.mat', 'distance_ratio_two'))
save(strcat(path, 'distance_ratio_three.mat', 'distance_ratio_three'))

%% Feature: CAI 

clear, clc
%load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
load('data_sets/challenge_data/codon_CAI.mat')
load('data_sets/feature_data/corresponding_orf.mat')

cai_reshaped = cell(1, 3);
cai_whole_reshaped = cell(1,3);

Sequences_ORF = corresponding_orf{1,2};
%Sequence_ORF_whole = whole_reshaped{1, 2};
cai_reshaped{1,2} = CAI_generator(Sequences_ORF,codon_CAI);
%cai_whole_reshaped{1, 2} = CAI_generator(Sequence_ORF_whole, codon_CAI);
%cai_ratio{1, 2} = cai_reshaped{1, 2}./cai_whole_reshaped{1, 2};

Sequences_UTR5 = corresponding_orf{1,1};
%Sequence_UTR5_whole = whole_reshaped{1, 1};
cai_reshaped{1, 1} = CAI_generator(Sequences_UTR5,codon_CAI);
% cai_whole_reshaped{1, 1} = CAI_generator(Sequence_UTR5_whole, codon_CAI);
% cai_ratio{1, 1} = cai_reshaped{1, 1}./cai_whole_reshaped{1, 1};

Sequences_UTR3 = corresponding_orf{1,3};
%Sequence_UTR3_whole = whole_reshaped{1, 3};
cai_reshaped{1, 3} = CAI_generator(Sequences_UTR3,codon_CAI);
% cai_whole_reshaped{1, 3} = CAI_generator(Sequence_UTR3_whole, codon_CAI);
% cai_ratio{1, 3} = cai_reshaped{1, 3}./cai_whole_reshaped{1, 3};


titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: CAI Score of Binding Window \n")

for i = 2:length(cai_reshaped)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(cai_reshaped{1, i}, reshaped_repress{1, i});
end

%fprintf("\nFeature: CAI Score of WHole Sequence \n")
% 
% for i = 2:length(cai_whole_reshaped)
%     fprintf("\nCurrent Sequence Type: %s", titles(i))
%     data_pipeline(cai_whole_reshaped{1, i}, reshaped_repress{1, i});
% end
% 
% fprintf("\nFeature: CAI Score Ratio \n")
% 
% for i = 2:length(cai_ratio)
%     fprintf("\nCurrent Sequence Type: %s", titles(i))
%     data_pipeline(cai_ratio{1, i}, reshaped_repress{1, i});
% end

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save('data_sets/feature_data/cai_reshaped.mat', 'cai_reshaped')
% save('data_sets/feature_data/cai_whole_reshaped.mat', 'cai_whole_reshaped')
% save('data_sets/feature_data/cai_ratio.mat', 'cai_ratio')

%% Feature: tAI 

clear, clc
load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
load('data_sets/challenge_data/codon_tAI.mat')
load('data_sets/feature_data/whole_sequence.mat')
load('data_sets/feature_data/corresponding_orf.mat')

tai_reshaped = cell(1, 3);
tai_whole_reshaped = cell(1,3);

Sequences_ORF = corresponding_orf{1,2};
%Sequence_ORF_whole = whole_reshaped{1, 2};
tai_reshaped{1,2} = CAI_generator(Sequences_ORF,codon_tAI);
%tai_whole_reshaped{1, 2} = CAI_generator(Sequence_ORF_whole, codon_tAI);
%tai_ratio{1, 2} = tai_reshaped{1, 2}./tai_whole_reshaped{1, 2};

Sequences_UTR5 = corresponding_orf{1,1};
%Sequence_UTR5_whole = whole_reshaped{1, 1};
tai_reshaped{1, 1} = CAI_generator(Sequences_UTR5,codon_tAI);
%tai_whole_reshaped{1, 1} = CAI_generator(Sequence_UTR5_whole, codon_tAI);
%tai_ratio{1, 1} = tai_reshaped{1, 1}./tai_whole_reshaped{1, 1};

Sequences_UTR3 = corresponding_orf{1,3};
%Sequence_UTR3_whole = whole_reshaped{1, 3};
tai_reshaped{1, 3} = CAI_generator(Sequences_UTR3,codon_tAI);
%tai_whole_reshaped{1, 3} = CAI_generator(Sequence_UTR3_whole, codon_tAI);
%tai_ratio{1, 3} = tai_reshaped{1, 3}./tai_whole_reshaped{1, 3};


titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: tAI Score of Binding Window \n")

for i = 2:length(tai_reshaped)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(tai_reshaped{1, i}, reshaped_repress{1, i});
end

% fprintf("\nFeature: tAI Score of WHole Sequence \n")
% 
% for i = 2:length(tai_whole_reshaped)
%     fprintf("\nCurrent Sequence Type: %s", titles(i))
%     data_pipeline(tai_whole_reshaped{1, i}, reshaped_repress{1, i});
% end
% 
% fprintf("\nFeature: tAI Score Ratio \n")
% 
% for i = 2:length(tai_ratio)
%     fprintf("\nCurrent Sequence Type: %s", titles(i))
%     data_pipeline(tai_ratio{1, i}, reshaped_repress{1, i});
% end



clearvars ans tAI_ORF tAI_UTR3 tAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_tAI 
save('data_sets/feature_data/tai_reshaped.mat', 'tai_reshaped')
% save('data_sets/feature_data/tai_whole_reshaped.mat', 'tai_whole_reshaped')
% save('data_sets/feature_data/tai_ratio.mat', 'tai_ratio')

%% GC content (Michal)

clear, clc

load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
load('data_sets/feature_data/whole_sequence.mat')
load('data_sets/feature_data/corresponding_orf.mat')
load('data_sets/feature_data/corresponding_utr5.mat')
load('data_sets/feature_data/corresponding_utr3.mat')

gc_orf = cell(1, 3);
gc_utr3 = cell(1, 3);
gc_utr5 = cell(1, 3);

%reshaped_nt_windows.mat is windows_reshaped
for i = 1:3
    gc_orf{1, i} = GC_content_generator(corresponding_orf{1,i});
    gc_utr5{1, i} = GC_content_generator(corresponding_utr5{1,i});
    gc_utr3{1, i} =  GC_content_generator(corresponding_utr3{1,i});
end
% 
% Sequences_UTR5 = corresponding_orf{1,1};
% gc_reshaped{1, 1} = GC_content_generator(Sequences_UTR5);
% % gc_whole_reshaped{1, 1} = GC_content_generator(whole_reshaped{1,1});
% % gc_ratio{1, 1} = gc_reshaped{1, 1} ./ gc_whole_reshaped{1, 1};
% 
% Sequences_UTR3 = corresponding_orf{1,3};
% gc_reshaped{1, 3} = GC_content_generator(Sequences_UTR3);
% % gc_whole_reshaped{1, 3} = GC_content_generator(whole_reshaped{1,3});
% % gc_ratio{1, 3} = gc_reshaped{1, 3} ./ gc_whole_reshaped{1, 3};



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

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save('data_sets/feature_data/gc_orf.mat', 'gc_orf')
save('data_sets/feature_data/gc_utr3.mat', 'gc_utr3')
save('data_sets/feature_data/gc_utr5.mat', 'gc_utr5')

%% MODELS!!!!!!!!!
%{
    It can be seen that the features with the best data are 

%}

clear, clc

% Repression
load('data_sets/feature_data/reshaped_repress.mat')

% CAI Scores
load('data_sets/feature_data/cai_reshaped.mat')
load('data_sets/feature_data/cai_whole_reshaped.mat')
load('data_sets/feature_data/cai_ratio.mat')

% GC content
load('data_sets/feature_data/gc_reshaped.mat')
load('data_sets/feature_data/gc_whole_reshaped.mat')
load('data_sets/feature_data/gc_ratio.mat')

% Terminus distances
load('data_sets/feature_data/terminus_distance_one.mat')
load('data_sets/feature_data/terminus_distance_two.mat')
load('data_sets/feature_data/distance_ratio_one.mat')
load('data_sets/feature_data/distance_ratio_two.mat')
load('data_sets/feature_data/distance_ratio_three.mat')
load('data_sets/feature_data/total_lengths.mat')
load('data_sets/feature_data/reshaped_indices.mat')

% Energies
load('data_sets/feature_data/folding_energies.mat')

%   tAI
load('data_sets/feature_data/tai_reshaped.mat')
load('data_sets/feature_data/tai_whole_reshaped.mat')
load('data_sets/feature_data/tai_ratio.mat')

% Conservation
load('data_sets/feature_data/conservations.mat')
load('data_sets/feature_data/conservation_ratios.mat')
load('data_sets/feature_data/whole_conservations.mat')%%

%% ORF
clc
clearvars info B nonan Xnonan MPGnonan

for i = 2
  y = reshaped_repress{i}';
    
    X_standard = [cai_reshaped{i}', conservation{i}', ...
         folding_energies{i}', gc_reshaped{i}', ...
         terminus_distance_two{i}', terminus_distance_one{i}', ...
         tai_reshaped{i}', reshaped_indices{i}'];
    
    
    X_whole = [cai_whole_reshaped{i}', ...
         folding_energies{i}', ...
         lengths_reshaped{i}', ...
         gc_whole_reshaped{i}',...
         tai_whole_reshaped{i}',...
         whole_conservations_reshaped{i}'];
    
   
    X_ratios = [cai_ratio{i}', conservation_ratios{i}', ...
         gc_ratio{i}', ...
         distance_ratio_two{i}', ...
         folding_energies{i}', ...
         tai_ratio{i}'];

    X_all = [cai_reshaped{i}', cai_ratio{i}', cai_whole_reshaped{i}', conservation{i}', ...
         folding_energies{i}', gc_reshaped{i}', conservation_ratios{i}', distance_ratio_one{i}', distance_ratio_two{i}', ...
         terminus_distance_two{i}', terminus_distance_one{i}', distance_ratio_three{i}', ...
         gc_ratio{i}', gc_whole_reshaped{i}', lengths_reshaped{i}',...
         tai_reshaped{i}', reshaped_indices{i}', tai_ratio{i}', tai_whole_reshaped{i}', whole_conservations_reshaped{i}'];
    
     X = X_standard;
     
    [B, info] = lasso(X, y, 'CV', 5, 'Alpha', 0.5);
    
    nonan = ~any(isnan([X y]),2);
    Xnonan = X(nonan,:);
    MPGnonan = y(nonan,:);
    corr(Xnonan)
    
    pnames = {'x1','x2','x3','x4', 'x5', 'x6', 'x7','x8','x9','x10', 'x11', 'x12', 'x13','x14','x15','x16', 'x17', 'x18', 'x19'};
    lassoPlot(B, info,'PlotType','Lambda','XScale','log',...
    'PredictorNames',pnames);
    

    coef = B(:,1);
    coef0 = info.Intercept(1);
    orf_model = struct('coef', coef, 'coef0', coef0);
    
    normalized_X = coef' * X';
    X_regularized = normalized_X' + coef0;
    
    correlation = corr(y, X_regularized) * 100;
    fprintf("\nPearson Lasso Correlation: %.2f%%\n", correlation)
    correlation = corr(y, X_regularized, 'type', 'spearman') * 100;
    fprintf("\nSpearman Lasso Correlation: %.2f%%\n", correlation)
    
%     orf_model = stepwiselm(X, y, 'linear');
%     y_pred = predict(orf_model, X);
%     
%     correlation = corr(y_pred, y) * 100;
%     fprintf("\nPearson Step Wise Regression Correlation: %.2f%%\n", correlation)
%     correlation = corr(y_pred, y, 'type', 'spearman') * 100;
%     fprintf("\nSpearman Step Wise Regression Correlation: %.2f%%\n", correlation)

end

save('regression_models/orf_model.mat', 'orf_model')

%% UTR3
clc
stepwise_model = cell(1);
lasso_model = cell(1);

for i = 3
   
    y = reshaped_repress{i}';
    
    X_standard = [cai_reshaped{i}', conservation{i}', ...
         folding_energies{i}', gc_reshaped{i}', ...
         terminus_distance_two{i}', terminus_distance_one{i}', ...
         tai_reshaped{i}', reshaped_indices{i}'];
    
    
    X_whole = [cai_whole_reshaped{i}', ...
         folding_energies{i}', ...
         lengths_reshaped{i}', ...
         gc_whole_reshaped{i}',...
         tai_whole_reshaped{i}',...
         whole_conservations_reshaped{i}'];
    
   
    X_ratios = [cai_ratio{i}', conservation_ratios{i}', ...
         gc_ratio{i}', ...
         distance_ratio_two{i}', ...
         folding_energies{i}', ...
         tai_ratio{i}'];

    X_all = [cai_reshaped{i}', cai_ratio{i}', cai_whole_reshaped{i}', conservation{i}', ...
         folding_energies{i}', gc_reshaped{i}', conservation_ratios{i}', distance_ratio_one{i}', distance_ratio_two{i}', ...
         terminus_distance_two{i}', terminus_distance_one{i}', distance_ratio_three{i}', ...
         gc_ratio{i}', gc_whole_reshaped{i}', lengths_reshaped{i}',...
         tai_reshaped{i}', reshaped_indices{i}', tai_ratio{i}', tai_whole_reshaped{i}', whole_conservations_reshaped{i}'];
    
     X = X_standard;
     
    [B, info] = lasso(X,y, 'CV', 10, 'Alpha', 0.5);
    
    coef = B(:,1);
    coef0 = info.Intercept(1);
    utr3_model = struct('coef', coef, 'coef0', coef0);
    
    normalized_X = coef' * X';
    y_pred = normalized_X' + coef0;
    
    correlation = corr(y_pred, y) * 100;
    fprintf("\nPearson Lasso Correlation: %.2f%%\n", correlation)
    correlation = corr(y_pred, y, 'type', 'spearman') * 100;
    fprintf("\nSpearman Lasso Correlation: %.2f%%\n", correlation)
    
%     utr3_model = stepwiselm(X, y, 'linear');
%     y_pred = predict(utr3_model, X);
%     
%     correlation = corr(y_pred, y) * 100;
%     fprintf("\nPearson Step Wise Regression Correlation: %.2f%%\n", correlation)
%     correlation = corr(y_pred, y, 'type', 'spearman') * 100;
%     fprintf("\nSpearman Step Wise Regression Correlation: %.2f%%\n", correlation)
% 
%     m = regress(y, X);
%     y_pred = X*m;
%     correlation = corr(y_pred, y) * 100;
%     fprintf("\nPearson Simple Regression Correlation: %.2f%%\n", correlation)
%     correlation = corr(y_pred, y, 'type', 'spearman') * 100;
%     fprintf("\nSpearman Simple Regression Correlation: %.2f%%\n", correlation)

end

save('regression_models/utr3_model.mat', 'utr3_model')






%% Extra

% Feature: Length of miRNA and repression (find average repression levels across each of 74 miRNAs)
clear


load(strcat(path, 'repress_use.mat'))
load(strcat(path, 'mirs_training_use.mat'))

mean_repress_miRNA = nanmean(repress_use, 2)';
mir_length = zeros(1, length(mirs_training));
for i = 1:length(mirs_training(2, :))
    mir_length(i) = strlength(mirs_training(2,i));
end

clearvars i mirs_training repress_use

fprintf("\nFeature: Length of miRNA vs Average Repression in all Genes")
data_pipeline(mir_length, mean_repress_miRNA);
save(strcat(path, 'mir_length.mat', 'mir_length'))
save(strcat(path, 'mean_repress_miRNA.mat', 'mean_repress_miRNA'));

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