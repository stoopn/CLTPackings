%% Probability of single occupancy K1 with uniform probability distribution

ndraws = 10^5;
fprintf('Computing probabilities \n')
fprintf(['N = ',num2str(neggchambers),' ... \n']);
t = cputime;
pK1_neggchambers = pval_K1(neggchambers,72,ndraws);
e = cputime - t;

fprintf(['OK in ',num2str(e),' s \n']);
fprintf(['N = ',num2str(nflatcysts),' ... \n'])
t = cputime;
pK1_flatcysts = pval_K1(nflatcysts,72,ndraws);
e = cputime - t;

fprintf(['OK in ',num2str(e),' s \n']);

%%
figure, 
subplot(1,2,2)  
bar(0:neggchambers,pK1_neggchambers);xlim([0 neggchambers+1])
title(['K1 probability distribution for N = ',num2str(neggchambers),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(14,pK1_neggchambers(15),'r')

subplot(1,2,1)
bar(0:nflatcysts,pK1_flatcysts);xlim([0 nflatcysts+1])
title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(22,pK1_neggchambers(23),'r')
