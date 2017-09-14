function [ greatest ] = lexico_order( seq1, seq2 )
%LEXICO_ORDER compares seq 1 and 2 in returns the indice of the greatest
%using the lexicographical order

i = 1;
eq = 1;
greatest = 2;
if (length(seq1)>0)&(length(seq2)>0),
    while (eq)&(i <= min(length(seq1),length(seq2))),
        if seq1(i)<seq2(i),
            greatest =  2;
            eq = 0;
        elseif seq1(i)>seq2(i),
            greatest = 1;
            eq = 0;
        end
        i = i+1;
    end
elseif (length(seq1)==0)&(length(seq2)==0),
    greatest = 2;
elseif length(seq2) == 0,
    greatest = 1;
end
end

