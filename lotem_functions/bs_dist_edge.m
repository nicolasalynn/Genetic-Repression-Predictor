function bs_dist_edge[dist] = untitled()
%seq is the mrna from genes_training
%bs_index keeps the location of the miRNA bs - only for those that have 1 bs
% x equals 2:5'UTR  3:ORF  4:3'UTR   depending on bs location

load('data_sets\feature_data\reshaped_indices.mat');
load('data_sets\feature_data\total_lengths.mat');

dist1 = [];
dist2 = [];
dist3 = [];

utr5 = reshaped_indices{1,1};
orf = reshaped_indices{1,2};
utr3 = reshaped_indices{1,3};

l_utr5 = lengths_reshaped{1,1};
l_orf = lengths_reshaped{1,2};
l_utr3 = lengths_reshaped{1,3};

num1 = size(utr5);
num2 = size(orf);
num3 = size(utr3);

    for i = 1:1:num1
      x1 = l_utr5(i);
      x2 = utr5(i);
      dist1 = [dist1,min(x1-x2, x2)];
    end
    
    for i = 1:1:num2
      x1 = l_orf(i);
      x2 = orf(i);
      dist2 = [dist2,min(x1-x2, x2)];
    end
    
    for i = 1:1:num3
      x1 = l_utr3(i);
      x2 = utr3(i);
      dist3 = [dist3,min(x1-x2, x2)];
    end

dist = {dist1,dist2, dist3};
end

