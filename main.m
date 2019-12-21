clear, clc
%%  Challenge Gene Expression Engineering
%   TAU Biomedical Engineering
%   Nicolas Lynn
%   Goal: To predict mRNA degradation and supression as a result of miRNA
%   interaction.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%  Create Data
clear, clc 

codon_weights = load('../data_sets/challenge_data/codon_weights.mat'); 

codon_CAI(1,:) = keys(codon_weights.CAI_weights);
codon_CAI(2,:) = values(codon_weights.CAI_weights);
codon_tAI(1,:) = keys(codon_weights.tAI_weights);
codon_tAI(2,:) = values(codon_weights.tAI_weights);

clearvars codon_weights

gene_training = load('../data_sets/challenge_data/genes_training.mat');
gene_training = gene_training.genes;
% a table of the training set genes, with the following columns:
% 1: ID ? The unique identifier of the gene
% 2: UTR5 (Sequences of the corresponding regions of the gene)
% 3: ORF (Sequences of the corresponding regions of the gene)
% 4: UTR3 (Sequences of the corresponding regions of the gene)
% 5: Conservation ? A conservation vector of each coordinate in the mRNA 
% (beginning at the first coordinate of the 5'UTR and ending at the last 
% coordinate of the 3'UTR). For additional information, see appendix 1.
% Literally, each nucleic acid in each of the three sections of DNA has a
% conservation value associated to it. 

miRs_training = load('../data_sets/challenge_data/miRs_training.mat');
% a map container with the training set miRNAs
% The keys are the miRNA names, as they appear in repress.mat.
% The values are the corresponding miRNA sequences.
mirs_training(1,:) = keys(miRs_training.miRs);
mirs_training(2,:) = values(miRs_training.miRs);

clearvars miRs_training
temp = load('../data_sets/challenge_data/repress.mat');
repress = temp.repress;
clearvars temp

%% Step 1: Find the first instance of miRNA mRNA binding for each combination & ...
% Generate an array of nucleotide 'windows'.

% find the miRNA strand that binds to mRNA (index 2:8) and convert this
% region to the reverse complement. Use this string as a referance to find
% the start index of its first occurance in each ORF

% 21 codons (63 nucleotides) windows pulled where the seed binds in the
% middle... So if the seed is 7 nts long and the window is 63 nts, the we
% need to grab 28 nucleotides before the first index and 35 nucleotides
% after the first index.

run_initiation = input("Do you want to recalculate the miRNA-mRNA binding "  +  ...
"indices? \n([Y] = 1, [N] = 0):  ");

if run_initiation
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training, gene_training)
    get_gene_windows(gene_training)
end


%% Load Data

load("binding_indices.mat")
load("nt_windows.mat")

%% Feature: Thermodynamics

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>")
if calc_folding_e == 1
    folding_energies = find_folding_energies(nt_windows);
end

%% Folding Energy



%% MER Site Distance to closest terminus



%% Cleaning Data
