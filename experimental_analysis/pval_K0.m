function [ proba_K0 ] = pval_K0( N, m, ntrials )
%PVAL_Max this function computes the probability of obtaining each values of
%the largest count for sample size N and number of bins m


proba_K0 = zeros(N+1,1);

for k = 1:ntrials,
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = length(find(r==0));
    proba_K0(itemp+1) = proba_K0(itemp+1) +1;
end
proba_K0 = proba_K0./ntrials;

end

