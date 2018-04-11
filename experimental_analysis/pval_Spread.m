function [ proba_Spread ] = pval_Spread( N, m, ntrials )
%PVAL_Max this function computes the probability of obtaining each values of
%the largest count for sample size N and number of bins m


proba_Spread = zeros(N+1,1);

for k = 1:ntrials,
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = max(r)-min(r);
    proba_Spread(itemp+1) = proba_Spread(itemp+1) +1;
end
proba_Spread = proba_Spread./ntrials;

end

