function [ eq ] = eq_seq( seq1, seq2 )
%EQ_SEQ Summary of this function goes here
%   Detailed explanation goes here
eq = 1;

if length(seq1)~=length(seq2),
    eq = 0;
else
    i = 1;
    while (eq)&(i <= length(seq1)),
        if seq1(i) ~= seq2(i),
            eq = 0;
        end
        i = i+1;
    end
end

end

