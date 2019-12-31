function [CAI_calculations] =  CAI_generator (Sequences,codon_CAI)
% codons_CAI contains the CAI weights per each codon
%nt_windows contains the relevant binding sites for each gene
%this function will calculate the CAI value of each sequence of the
%binding sites

[r,Sequences_number]= size(Sequences);
j=0;
    for j=1:Sequences_number
        %We take the binding site sequence in column j and calcuate the CAI value of the sequence there
        seq=Sequences(j); 
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
                [row,column] = find(strcmp(codon_CAI, str));            
                CAI_str_value=codon_CAI{2,column};
                w(counter)=CAI_str_value;
                i=i+3;
        end
        CAI(j)=exp(sum(log(w))/codon_number);
    end
CAI_calculations = CAI;
end



% function [CAI_calculations] =  CAI_generator (nt_windows,codon_CAI)
% % codons_CAI contains the CAI weights per each codon
% %nt_windows contains the relevant binding sites for each gene
% %this function will calculate the CAI value of each sequence of the
% %binding sites
% 
% [miRNA_length,gene_length] = size(nt_windows);
% z=0;
% j=0;
% for z=1:gene_length
%     for j=1:miRNA_length
%         %We take the binding site sequence for mRNA binding site z and
%         %miRNA j and calcuate the CAI value of the sequence there
%         seq=nt_windows(j,z); 
%         length_seq=strlength(seq);
%         codons_length=length_seq/3;
%         codon_number=floor(codons_length); %taking a sequence that can be divided by 3
%         w=1;
%         i=1;
%         counter=0;
%         while i<=(codon_number-2) 
%                 counter=counter+1;
%                 seq_char=char(seq);
%                 str=seq_char(i:i+2);            
%                 [row,column] = find(strcmp(codon_CAI, str));            
%                 CAI_str_value=codon_CAI{2,column};
%                 w(counter)=CAI_str_value;
%                 i=i+3;
%         end
%         CAI(j,z)=exp(sum(log(w))/codon_number);
%     end
% end
% CAI_calculations = CAI;
% end