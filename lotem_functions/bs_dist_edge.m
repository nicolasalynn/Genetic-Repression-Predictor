function bs_dist_edge[dist] = untitled(seq,bs_index, x)
%seq is mrna from genes_training
%bs_index keeps the location of the miRNA bs - only for those that have 1 bs
% x equeals 2:5'UTR  3:ORF  4:3'UTR   depending on bs location

load('genes_training.mat');
L = length(seq);

for i=1:1:L
    if(x==2)
        
    elseif(x==3)
    elseif(x==4)
    end
    
end


dist = dist;
end

