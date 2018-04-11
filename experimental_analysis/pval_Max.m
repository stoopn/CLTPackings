function [ proba_Max ] = pval_Max( N, m, ntrials )
%PVAL_Max this function computes the probability of obtaining each values of
%the largest count for sample size N and number of bins m


proba_Max = zeros(N+1,1);

for k = 1:ntrials,
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = max(r);
    proba_Max(itemp+1) = proba_Max(itemp+1) +1;
end
proba_Max = proba_Max./ntrials;

end

