function [ proba_K5 ] = pval_K5( N, m, ntrials )
%PVAL_Max this function computes the probability of obtaining each values of
%the largest count for sample size N and number of bins m


proba_K5 = zeros(ceil(N/5),1);

for k = 1:ntrials,
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = length(find(r>=5));
    proba_K5(itemp+1) = proba_K5(itemp+1) +1;
end
proba_K5 = proba_K5./ntrials;

end

