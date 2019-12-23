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


%%  Create Data
% The data we have obtained contains data for ~5000 genes and ~200 miRNAs.
% Not all gene miRNA combinations bind. We only wish to look at the
% combinations with binding. 


clear, clc 

codon_weights = load('data_sets/challenge_data/codon_weights.mat'); 

codon_CAI(1,:) = keys(codon_weights.CAI_weights);
codon_CAI(2,:) = values(codon_weights.CAI_weights);
codon_tAI(1,:) = keys(codon_weights.tAI_weights);
codon_tAI(2,:) = values(codon_weights.tAI_weights);

clearvars codon_weights

gene_training = load('data_sets/challenge_data/genes_training.mat');
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

miRs_training = load('data_sets/challenge_data/miRs_training.mat');
% a map container with the training set miRNAs
% The keys are the miRNA names, as they appear in repress.mat.
% The values are the corresponding miRNA sequences.
mirs_training(1,:) = keys(miRs_training.miRs);
mirs_training(2,:) = values(miRs_training.miRs);

clearvars miRs_training
temp = load('data_sets/challenge_data/repress.mat');
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
"indices? This action will take approximatelly 2 minutes... \n([Y] = 1, [N] = 0):  ");

if run_initiation 
    %tic;
    fprintf("\nThis will take a minute...\n\n");
    binding_indices(mirs_training, gene_training)
    load('data_sets/feature_data/binding_indices.mat')
    previewData(first_indices, 10);
    %initiation_time = toc;
end

% site_indices will be a 74 x 3947 x 3  matrix with the index values for
% the orrucance of miRNA in gene (1D: UTR5, 2D: ORF, 3D: UTR3

run_windows = input("Do you want to recalculate the binding windows? " + ...
    "\n([Y] = 1, [N] = 0):  ");

if run_windows
    tic;
    fprintf("\nTHis will take a minute....\n\n");
    get_gene_windows(gene_training)
    genwindows_time = toc;
end

clearvars run_initiation run_windows
%% Load Data

load("binding_indices.mat")
load("nt_windows.mat")

%% Feature: Thermodynamics

calc_folding_e = input("\nWould you like to calculate folding " + ...
    "energies?\nThis will take a few minutes..\n [Y]:1, [N]:0\n>>");
if calc_folding_e == 1
    folding_energies = find_folding_energies(nt_windows);
end

clearvars calc_folding_e

%% Exploring data (what is the averate repression level where there are miRNA binding sites and were there arent)

%non_mirna_binding_repression = binding_average_repress(repress, nt_windows, 'nb');
%mirna_binding_repression = binding_average_repress(repress, nt_windows, 'b');


%% Folding energy of a window around the binding side of the miRNA

%disp(all_indices(:,:,3));

%% Folding Energy (calculate folding energy of the beginning of each ORF and 5'NCR



%% MER Site Distance to closest terminus


%% tAI (Michal)

tAI = tAI_generator(codon_tAI);


%% Cleaning Data
