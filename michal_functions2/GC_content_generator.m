function [GC_content_calculation] = GC_content_generator(nt_windows)
%nt_windows contains the relevant binding sites for each gene
%this function will calculate the GC content of each sequence of the
%binding sites

[miRNA_length,gene_length] = size(nt_windows);
z=0;
j=0;
for z=1:gene_length
    for j=1:miRNA_length
        seq=nt_windows(j,z); 
        A_number=length(strfind(seq,'A'));
        T_number=length(strfind(seq,'T'));
        G_number=length(strfind(seq,'G'));
        C_number=length(strfind(seq,'C'));
        All_nucleotides=A_number+T_number+G_number+C_number;
        GC_nucleotides=G_number+C_number;
        GC_cont(j,z)=(GC_nucleotides)/(All_nucleotides);
    end
end
GC_content_calculation=GC_cont;
end