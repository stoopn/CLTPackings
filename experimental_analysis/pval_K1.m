function [ proba_K1 ] = pval_K1( N, m, ntrials )
%PVAL_K1 this function computes the probability of obtaining each values of
%K1 for sample size N and number of bins m


proba_K1 = zeros(N+1,1);

for k = 1:ntrials,
    %fprintf([num2str(k),'\n']);
    r = mnrnd(N,1/m*ones(1,m));
    itemp = length(find(r==1));
    proba_K1(itemp+1) = proba_K1(itemp+1) +1;
end
proba_K1 = proba_K1./ntrials;

end

