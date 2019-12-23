clear, clc
%%  Genetic Supression Predictor
%   Goal: To predict mRNA degradation and supression as a result of miRNA
%   interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
    Technical notes:
                        size(gene_training, 1) = 2947
                        size(gene_training, 2) = 5
                        length(mirtna_training) = 74
%}

addpath nico_functions
addpath lotem_functions
addpath michal_functions

%%  Pull Data  --- THIS DOES NOT HAVE TO BE TOUCHED

clear, clc 

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

%% Step 1: Find the first instance of miRNA mRNA binding for each combination & ...

run_initiation = input("Do you want to recalculate the miRNA-mRNA binding "  +  ...
"indices? This action will take approximatelly 2 minutes... \n([Y] = 1, [N] = 0):  ");

if run_initiation 
    tic;
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training, gene_training, repress)
    load('data_sets/feature_data/binding_indices.mat')
    initiation_time = toc;
end

run_windows = input("Do you want to recalculate the binding windows? " + ...
    "\n([Y] = 1, [N] = 0):  ");

if run_windows
    tic;
    load('data_sets/feature_data/true_indices.mat');
    fprintf("\nTHis will take a minute....\n\n");
    get_gene_windows(gene_training, true_indices);
    genwindows_time = toc;
end

clearvars run_initiation run_windows
%% Load Data

load("data_sets/feature_data/binding_indices.mat")
load("data_sets/feature_data/nt_windows.mat")
load("data_sets/feature_data/all_indices.mat")
load("data_sets/challenge_data/repress.mat")
repress = table2array(repress(:, 2:end))';
%% Feature: Number of Binding Sites Across all regions (Nico)
combined_indices = all_indices(:, :, 1) + all_indices(:, :, 2) + all_indices(:, :, 3); % number of occurances accross all three sequences
[M, I] = max(combined_indices);
previewData(combined_indices, 10);

combined_indices(combined_indices == 0) = NaN;
data = create_usable_data(combined_indices, repress);
m = regress(data(1, :)', data(2, :)');

%% Feature: Thermodynamics

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>");
if calc_folding_e == 1
    tic
    folding_energies = find_folding_energies(nt_windows);
    fold_energy_time = toc;
end

clearvars calc_folding_e

%% Exploring data (what is the averate repression level where there are miRNA binding sites and were there arent)

%non_mirna_binding_repression = binding_average_repress(repress, nt_windows, 'nb');
%mirna_binding_repression = binding_average_repress(repress, nt_windows, 'b');


%% MER Site Distance to closest terminus 


%% tAI (Michal)

tAI = tAI_generator(codon_tAI);

%% GC content (Michal)


%% Cleaning Data
