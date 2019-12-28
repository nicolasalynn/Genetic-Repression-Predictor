function [conservation] = Conservation(seq,bs_index)
% seq is 1 mrna
% bs_index is the 1st coordinate of the binding site

load('genes_training.mat');
%cons = genes_training.conservation;
%L = length(bs_indices);
avg = 0;
conservation = [];


    %a = bs_indices(j); %1st binding site index
    a = bs_index;
    for i=0:1:6
        avg = avg + seq.conservation(i+a); %???? ????? ??????? ?? ????? 
    end
    %conservation = [conservation, avg/7];


conservation = (avg/7);
end

