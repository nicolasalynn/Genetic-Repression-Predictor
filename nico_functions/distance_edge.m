function [distance_cell_end, disatnce_cell_tot] = distance_edge()

reshaped_indices = load('data_sets/feature_data/reshaped_indices.mat');
total_lengths = load('data_sets/feature_data/total_lengths.mat');

dim = length(total_lengths.lengths_reshaped);

distance_cell_end = cell(1, dim);
disatnce_cell_tot = cell(1, dim);

for i = 1:dim
    length_data = total_lengths.lengths_reshaped{1, i};
    index_data = reshaped_indices.reshaped_indices{1, i};
    distances = zeros(size(index_data));
    distances_two = zeros(size(index_data));
    distances(1, :) = abs(length_data - index_data);
    distances_two(1,:) = min(distances, index_data);
    distance_cell_end{1, i} = distances;
    disatnce_cell_tot{1, i} = distances_two;
end

    
    
    