function bs_dist_edge[dist] = untitled(seq,bs_index, x)
%seq is the mrna from genes_training
%bs_index keeps the location of the miRNA bs - only for those that have 1 bs
% x equals 2:5'UTR  3:ORF  4:3'UTR   depending on bs location

load('genes_training.mat');

    if(x==2)
      num1 = length(seq.utr5);
      num2 = bs_index
      dist = min(num1-num2, num2);
    elseif(x==3)
      num1 = length(seq.orf);
      num2 = bs_index
      dist = min(num1-num2, num2);
    elseif(x==4)
      num1 = length(seq.utr3);
      num2 = bs_index
      dist = min(num1-num2, num2);
    end

dist = dist;
end
