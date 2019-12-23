% function [CAI_calculations] =  CAI_generator (nt_windows,codon_CAI)
% codons_CAI contains the CAI weights per each codon
%nt_windows contains the relevant binding sites for each gene

codon_CAI(1,:) = keys(CAI_weights);
codon_CAI(2,:) = values(CAI_weights);

Gene_length = size(nt_windows, 1);
miRNA_length = size(nt_windows, 2);

z=0;
j=0;

for z=1:Gene_length
    for j=1:miRNA_length
        %We take the binding site sequence for gene z and miRNA j and calcuate the CAI value of this sequence
        seq=nt_windows(z,j); 
        length_seq=strlength(seq);
        codons_length=length_seq/3;
        codon_number=floor(codons_length); %taking a sequence that can be divided by 3
        w=1;
        i=1;
        counter=0;
        while i<=(codon_number-2) 
                counter=counter+1;
                seq_char=char(seq);
                str=seq_char(i:i+2);
                CAI_str_index=find(codon_CAI(1,:)==str);
                CAI_str_value=codon_CAI(2,CAI_str_index(2));
                w(counter)=CAI_str_value;
                i=i+3;
        end
        CAI(j,z)=exp(sum(log(w))/codon_number);
    end
end
CAI_calculations = CAI;
% end