function [ proba_K6 ] = pval_K6( N, m, ntrials )
%PVAL_Max this function computes the probability of obtaining each values of
%the largest count for sample size N and number of bins m


proba_K6 = zeros(ceil(N/6),1);

for k = 1:ntrials
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = length(find(r>=6));
    proba_K6(itemp+1) = proba_K6(itemp+1) +1;
end
proba_K6 = proba_K6./ntrials;

end

