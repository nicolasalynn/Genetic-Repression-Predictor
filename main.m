clear, clc
%%  Genetic Supression Predictor -- TRAINING -- RUN ME
%   Goal: To predict mRNA degradation and supression as a result of miRNA
%   interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\

%% Training Initiation
addpath nico_functions
addpath lotem_functions
addpath michal_functions

%  Pull Data  --- THIS DOES NOT HAVE TO BE TOUCHED

clear, clc 
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

%% Find the first instance of miRNA mRNA binding for each combination (Nico)

run_initiation = input("Do you want to recalculate the miRNA-mRNA binding "  +  ...
"indices? This action will take approximatelly 2 minutes... \n([Y] = 1, [N] = 0):  ");

if run_initiation 
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training(2, :), gene_training, repress, 'data_sets/feature_data/')
end
clearvars run_initiation

%% Obtain windows of specified length for all indices found previously
run_windows = input("Do you want to recalculate the binding windows? " + ...
    "\n([Y] = 1, [N] = 0):  ");

if run_windows
    load('data_sets/feature_data/true_indices.mat');
    fprintf("\nThis might take a minute....\n\n");
    get_gene_windows(gene_training, true_indices, 'nt_windows', 74, "training"); %by default, set to 74
end
clearvars run_windows

%% Load Data -- RUN THIS IF YOU ARE NOT INITIATING

clear, clc
addpath nico_functions
addpath lotem_functions
addpath michal_functions

load("data_sets/feature_data/reshaped_repress.mat")
load("data_sets/feature_data/reshaped_nt_windows.mat")
load("data_sets/feature_data/reshaped_indices.mat")


%load("data_sets/feature_data/true_indices.mat")
%load("data_sets/feature_data/nt_windows.mat")
%load("data_sets/feature_data/all_indices.mat")
%load("data_sets/feature_data/good_repress.mat")
%load("data_sets/feature_data/binary_truth.mat")

%% Feature: Number of Binding Sites Across all regions (Nico)
clear, clear, clc

load("data_sets/feature_data/all_indices.mat")
load("data_sets/challenge_data/repress_use.mat")

combined_indices = all_indices(:, :, 1) + all_indices(:, :, 2) + all_indices(:, :, 3); % number of occurances accross all three sequences
repress = repress_use;
clearvars all_indices repress_use

fprintf("\nFeature: Number of binding sites across UTR5, ORF, UTR3")
data_pipeline(combined_indices, repress);

%% Feature: Thermodynamics
clear, clc
load("data_sets/feature_data/reshaped_nt_windows.mat")

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>");
if calc_folding_e == 1
    
    %tic
    dim = 0; % change this value depending of sequence region target
    folding_energies = find_folding_energies(windows_reshaped, "training");
    %fold_energy_time = toc;  
end
clearvars calc_folding_e
  
load('data_sets/feature_data/folding_energies.mat');
load('data_sets/feature_data/reshaped_repress.mat');

for dim = 1:size(folding_energies, 2)
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

%% Feature: Average Repression in presence and absence of binding site
clear, clc

load("data_sets/feature_data/all_indices.mat")
load("data_sets/challenge_data/repress_use.mat")

binding_or_no = all_indices;
binding_or_no(binding_or_no > 0) = 1;
binding_or_no(binding_or_no ~= 1) = 0;
clearvars all_indices 
fprintf("\nFeature: Mean Observations Between Presence of Binding Site and None")
data_pipeline(binding_or_no(:,:,1), repress_use);

%% Feature: Length of miRNA and repression (find average repression levels across each of 74 miRNAs)
clear, clc

load('data_sets/challenge_data/repress_use.mat')
load('data_sets/challenge_data/mirs_training_use.mat')

mean_repress_miRNA = nanmean(repress_use, 2)';
mir_length = zeros(1, length(mirs_training));
for i = 1:length(mirs_training(2, :))
    mir_length(i) = strlength(mirs_training(2,i));
end

clearvars i mirs_training repress_use

fprintf("\nFeature: Length of miRNA vs Average Repression in all Genes")
data_pipeline(mir_length, mean_repress_miRNA);
save('data_sets/feature_data/mir_length.mat', 'mir_length')
save('data_sets/feature_data/mean_repress_miRNA.mat', 'mean_repress_miRNA');

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
    data_pipeline(seq_lengths(i, :), mean_repress_gene);
end

save('data_sets/feature_data/seq_lengths.mat', 'seq_lengths')
save('data_sets/feature_data/mean_repress_gene.mat', 'mean_repress_gene')
clearvars gene_training i j repress_use sequences titles ans

%% Feature: Conservation
clear, clc

load('data_sets/feature_data/conservations.mat')
load('data_sets/feature_data/reshaped_repress.mat')
titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: Average Gene Binding Window Conservation\n")

for i = 1:length(conservation)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(conservation{1, i}, reshaped_repress{1, i});
end

%% Feature: Distance to terminus

clear, clc

load('data_sets/feature_data/reshaped_repress.mat')
load('data_sets/feature_data/reshaped_indices.mat');
load('data_sets/feature_data/total_lengths.mat');

[terminus_distance_one, terminus_distance_two] = distance_edge(reshaped_indices, lengths_reshaped, "training");
reshaped_repress = [reshaped_repress reshaped_repress];
terminus_distance = [terminus_distance_one terminus_distance_two];
clearvars terminus_distance_one terminus_distance_two

titles = ["Distance from End, UTR5", "Distance from End, ORF", "Distance from End, UTR3" ...
    "Distance from Either, UTR5", "Distance from Either, ORF", "Distance from Either, UTR3"];

fprintf("\nFeature: Distance to Closest Terminus\n")
for i = 1:length(terminus_distance)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(terminus_distance{1, i}, reshaped_repress{1, i});
end

clearvars i ans titles
save('data_sets/feature_data/terminus_distance.mat', 'terminus_distance')

%% %%%Feature: MER Site Distance to closest terminus 
% clear, clc
% 
% load('data_sets/feature_data/reshaped_indices.mat');
% load('data_sets/feature_data/total_lengths.mat');
% 
% x = bs_dist_edge();
% dist1 = x{1,1};
% dist2 = x{1,2};
% dist3 = x {1,3};
% 
% repress_dist_utr5 = reshaped_indices{1,1};
% data_pipeline(dist1, repress_dist_utr5);
% repress_dist_orf = reshaped_indices{1,2};
% data_pipeline(dist2, repress_dist_orf);
% repress_dist_utr3 = reshaped_indices{1,3};
% data_pipeline(dist3, repress_dist_utr3);

%% CAI (Michal)

clear, clc
load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
load('data_sets/challenge_data/codon_CAI.mat')

%reshaped_nt_windows.mat is windows_reshaped
Sequences_ORF = windows_reshaped{1,2};
CAI_ORF = CAI_generator(Sequences_ORF,codon_CAI);
Sequences_UTR5 = windows_reshaped{1,1};
CAI_UTR5 = CAI_generator(Sequences_UTR5,codon_CAI);
Sequences_UTR3 = windows_reshaped{1,3};
CAI_UTR3 = CAI_generator(Sequences_UTR3,codon_CAI);

cai_reshaped = cell(1, 3);
cai_reshaped{1, 1} = CAI_UTR5;
cai_reshaped{1, 2} = CAI_ORF;
cai_reshaped{1,3} = CAI_UTR3;

titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: CAI Score of Binding Window \n")

for i = 1:length(cai_reshaped)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(cai_reshaped{1, i}, reshaped_repress{1, i});
end

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save('data_sets/feature_data/cai_reshaped.mat', 'cai_reshaped')

%% tAI (Nico)


clear, clc
load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
load('data_sets/challenge_data/codon_tAI.mat')

Sequences_ORF = windows_reshaped{1,2};
CAI_ORF = CAI_generator(Sequences_ORF,codon_tAI);
Sequences_UTR5 = windows_reshaped{1,1};
CAI_UTR5 = CAI_generator(Sequences_UTR5,codon_tAI);
Sequences_UTR3 = windows_reshaped{1,3};
CAI_UTR3 = CAI_generator(Sequences_UTR3,codon_tAI);

repress_CAI_ORF = reshaped_repress{1,2};
data_pipeline(CAI_ORF, repress_CAI_ORF);        %   1%
repress_CAI_UTR5 = reshaped_repress{1,1};
data_pipeline(CAI_UTR5, repress_CAI_UTR5);      %   9%  winner
repress_CAI_UTR3 = reshaped_repress{1,3};
data_pipeline(CAI_UTR3, repress_CAI_UTR3);      %   .7%

cai_reshaped = cell(1, 3);
cai_reshaped{1, 1} = CAI_UTR5;
cai_reshaped{1, 2} = CAI_ORF;
cai_reshaped{1,3} = CAI_UTR3;

%% GC content (Michal)

clear, clc

load('data_sets/feature_data/reshaped_nt_windows.mat');
load('data_sets/feature_data/reshaped_repress.mat');
%reshaped_nt_windows.mat is windows_reshaped
Sequences_ORF = windows_reshaped{1,2};
GC_content_ORF = GC_content_generator(Sequences_ORF);
Sequences_UTR5 = windows_reshaped{1,1};
GC_content_UTR5 = GC_content_generator(Sequences_UTR5);
Sequences_UTR3 = windows_reshaped{1,3};
GC_content_UTR3 = GC_content_generator(Sequences_UTR3);


gc_reshaped = cell(1, 3);
gc_reshaped{1, 1} = GC_content_UTR5;
gc_reshaped{1, 2} = GC_content_ORF;
gc_reshaped{1,3} = GC_content_UTR3;

titles = ["UTR5", "ORF", "UTR3"];

fprintf("\nFeature: GC Score of Binding Window \n")

for i = 1:length(gc_reshaped)
    fprintf("\nCurrent Sequence Type: %s", titles(i))
    data_pipeline(gc_reshaped{1, i}, reshaped_repress{1, i});
end

clearvars ans CAI_ORF CAI_UTR3 CAI_UTR5 Sequences_ORF Sequences_UTR3 Sequences_UTR5 titles windows_reshaped i codon_CAI 
save('data_sets/feature_data/gc_reshaped.mat', 'gc_reshaped')

%% NOW WE WILL TRY COMBINING DATA TO GET BETTER PREDICTIVE RESULTS
%{
    It can be seen that the features with the best data are 

%}

clear, clc

load('data_sets/feature_data/reshaped_repress.mat')
load('data_sets/challenge_data/gene_training_use.mat')
load('data_sets/feature_data/whole_sequence.mat')
load('data_sets/feature_data/cai_reshaped.mat')
load('data_sets/feature_data/conservations.mat')
load('data_sets/feature_data/gc_reshaped.mat')
load('data_sets/feature_data/terminus_distance.mat')
load('data_sets/feature_data/total_lengths.mat')
load('data_sets/feature_data/reshaped_indices.mat')
load('data_sets/feature_data/folding_energies.mat')

total_lengths = cell(1, 3);
for i = 1:3
    seqs = whole_reshaped{i};
    lengths = zeros(1, length(seqs));
    for j = 1:length(seqs)
        lengths(j) = strlength(seqs(j));
    end
    total_lengths{i} = lengths;
    data_pipeline(total_lengths{i},reshaped_repress{i}); 

end
%%
clc

for i = 1:3
   
    y = reshaped_repress{i}';
    X = [cai_reshaped{i}', conservation{i}', gc_reshaped{i}',reshaped_indices{i}', terminus_distance{i}', folding_energies{i}'];
    
    [B, info] = lasso(X,y);
    
    coef = B(:,1);
    coef0 = info.Intercept(1);
    
    normalized_X = coef' * X';
    y_pred = normalized_X' + coef0;
    
    correlation = corr(y_pred, y) * 100;
    fprintf("Correlation: %.2f%%\n", correlation)
    
end


