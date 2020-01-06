function [CAI_calculations] =  CAI_generator (Sequences,codon_CAI)
% codons_CAI contains the CAI weights per each codon
%Sequences contains the relevant binding sites
%this function will calculate the CAI value of each sequence of the
%binding sites

[~,Sequences_number]= size(Sequences);
    for j=1:Sequences_number
        %We take the binding site sequence in column j and calcuate the CAI value of the sequence there
        seq=Sequences(j); 
        length_seq=strlength(seq);
        codons_length=length_seq/3;
        codon_number=floor(codons_length); %taking a sequence that can be divided by 3
        w=1;
        i=3;
        counter=0;
       while i<=(codon_number*3 - 2) 
                counter=counter+1;
                seq_char=char(seq);
                str=seq_char(i:i+2);            
                [~,column] = find(strcmp(codon_CAI, str));            
                CAI_str_value=codon_CAI{2,column};
                w(counter)=CAI_str_value;
                i=i+3;
        end
        %CAI(j)=exp(sum(log(w))/codon_number);
        CAI(j) = geomean(w);
    end
CAI_calculations = CAI;
end




