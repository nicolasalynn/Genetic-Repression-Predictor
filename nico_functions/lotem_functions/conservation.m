function [conservation] = conservation(seq,bs_indices)

load('genes_training.mat');
cons = genes_training.conservation;
L = length(bs_indices);
avg = 0;
conservation = [];

for j=1:1:L
for i=0:1:6
    avg = avg + seq(j+i);
end
conservation = [conservation, avg/7];
avg = 0;
end

conservation = conservation;
end

