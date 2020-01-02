%% Create model here
clear, clc

%utr5 folding energies
%utr3 conservation
%distances to terminus
%utr3 cai

load('data_sets/feature_data/conservations.mat')
load('data_sets/feature_data/folding_energies.mat')
load('data_sets/feature_data/cai_reshaped.mat')
load('data_sets/feature_data/terminus_distance.mat')
load('data_sets/feature_data/reshaped_repress.mat')
load('data_sets/feature_data/reshaped_indices.mat')


%% UTR5' Data
X = [conservation{1, 1}' folding_energies{1, 1}' terminus_distance{1, 1}' cai_reshaped{1, 1}' reshaped_indices{1,1}'];
Y = reshaped_repress{1, 1}';
m = regress(Y, X);
tog = [X Y];
y_predict =  m' * X';
mstep = stepwiselm(tog, 'linear');
y_predict_step = predict(mstep, X);


correlation = corr(Y, y_predict')


%% ORF Data

X = [conservation{1, 2}' folding_energies{1, 2}' terminus_distance{1, 2}' cai_reshaped{1, 2}' reshaped_indices{1,2}'];
Y = reshaped_repress{1, 2}';
m = regress(Y, X);
y_predict =  m' * X';

correlation = corr(Y, y_predict')

%% UTR3' Data

X = [conservation{1, 3}' folding_energies{1, 3}' terminus_distance{1, 3}' cai_reshaped{1, 3}' reshaped_indices{1,3}'];
Y = reshaped_repress{1, 3}';
m = regress(Y, X);
y_predict =  m' * X';

correlation = corr(Y, y_predict')

%% Other Data
clear, clc

load('data_sets/feature_data/orf_length.mat')
X = orf_length';
load('data_sets/feature_data/mean_repress_gene.mat')
Y = mean_repress_gene';

m = regress(Y, X);
y_predict =  m' * X';

correlation = corr(Y, y_predict')

%%
clear, clc

load('data_sets/feature_data/mir_length.mat')
load('data_sets/feature_data/mean_repress_miRNA.mat');
X = mir_length';
Y = mean_repress_miRNA';
m = regress(Y, X);
y_predict =  m' * X';

correlation = corr(Y, y_predict')
