function [GC_content_calculation] = GC_content_generator(Sequences)
%Sequences contains the relevant binding sites
%this function will calculate the GC content of each sequence of the
%binding sites

[r,Sequences_number]= size(Sequences);
j=0;
    for j=1:Sequences_number
        seq=Sequences(j); 
        A_number=length(strfind(seq,'A'));
        T_number=length(strfind(seq,'T'));
        G_number=length(strfind(seq,'G'));
        C_number=length(strfind(seq,'C'));
        All_nucleotides=A_number+T_number+G_number+C_number;
        GC_nucleotides=G_number+C_number;
        GC_cont(j)=(GC_nucleotides)/(All_nucleotides);
    end
GC_content_calculation=GC_cont;
end