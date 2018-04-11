clear all;
include_flat = true;
cd '/Users/norbert/Dropbox (MIT)/StasPacking/Data/experiment'
figure_dir = '/Users/norbert/Dropbox (MIT)/StasPacking/figure_data/'

rounded = load('ExpStates_rounded_all.mat');
% Uncomment to include flat states:
if include_flat == true
flat = load('ExpStates_flat_data1.mat');
flat_prefix = 'r+f';
else
    flat_prefix = '';
end
%%
ndraws = 10^5;
totalcounts = zeros(2, 72);
% Rounded samples
[totalcounts_r,EDGES] = histcounts(rounded.ExpStates,0.5:72.5);
totalcounts(1,:)=totalcounts_r;
if include_flat
% Flat samples
[totalcounts_f,EDGES] = histcounts(flat.ExpStates,0.5:72.5);
totalcounts(2,:)=totalcounts_f;
end
% Sum them up
totalcounts = sum(totalcounts);

nsingle = size(find(totalcounts==1),2);
nK2 = size(find(totalcounts==2),2);
nK0 = size(find(totalcounts==0),2);
nK4 = size(find(totalcounts==4),2);
nK5 = size(find(totalcounts>=5),2);
nK6 = size(find(totalcounts>=6),2);

nmax=max(totalcounts);
nSpread = max(totalcounts)-min(totalcounts);
ntotal = sum(totalcounts); % = all data counts

%%
fprintf('Computing K0 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK0_ntotal = pval_K0(ntotal,72,ndraws);
toc;

fprintf('Computing K1 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK1_ntotal = pval_K1(ntotal,72,ndraws);
toc;

fprintf('Computing max probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pMax_ntotal = pval_Max(ntotal,72,ndraws);
toc;

fprintf('Computing Spread probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pSpread_ntotal = pval_Spread(ntotal,72,ndraws);
toc;

fprintf('Computing K2 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK2_ntotal = pval_K2(ntotal,72,ndraws);
toc;

fprintf('Computing K4 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK4_ntotal = pval_K4(ntotal,72,ndraws);
toc;

fprintf('Computing K5 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK5_ntotal = pval_K5(ntotal,72,ndraws);
toc;

fprintf('Computing K6 probabilities \n')
fprintf(['N = ',num2str(ntotal),' ... \n']);
tic;
pK6_ntotal = pval_K6(ntotal,72,ndraws);
toc;
%%
close all;
figure(1);
hold off;
histogram('BinEdges', -0.5+[0:ntotal+1],'BinCounts', pK1_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0 55]);
ylim([0,0.107]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',24)
set(gcf, 'color', 'white');
hold on,
bar(nsingle,pK1_ntotal(nsingle+1),1,'r');
axis square;
set(gca,'linewidth', 2);
export_fig([figure_dir, 'probK1_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');


figure(2);
histogram('BinEdges', -0.5+[0:ntotal+1],'BinCounts', pMax_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0 10]);
if include_flat
    xlim([0,15]);
end
ylim auto;
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(nmax,pMax_ntotal(nmax+1),1,'r');
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);
axis square;
export_fig([figure_dir, 'probMax_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');
%print(gcf, '-dpdf', 'probK1_both.pdf');


figure(3);
histogram('BinEdges', -0.5+[0:ceil(ntotal/2)],'BinCounts', pK2_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0 30]);
%ylim([0,0.13]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(nK2,pK2_ntotal(nK2+1),1,'r');
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);
axis square;
export_fig([figure_dir, 'probK2_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');


figure(4);
histogram('BinEdges', -0.5+[0:ntotal+1],'BinCounts', pK0_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0 45]);
ylim([0,0.17]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(nK0,pK0_ntotal(nK0+1),1,'r');
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);
axis square;
export_fig([figure_dir, 'probK0_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');


figure(5);
histogram('BinEdges', -0.5+[0:ntotal+1],'BinCounts', pSpread_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0, 10]);
ylim([0,0.65]);
if include_flat
    xlim([0,13]);
    ylim([0,0.55]);
end
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(nSpread,pSpread_ntotal(nSpread+1),1, 'r');
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);
axis square;
export_fig([figure_dir, 'probSpread_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');

close all;
figure(6);
hold off;
histogram('BinEdges', -0.5+[0:ceil(ntotal/4)],'BinCounts', pK4_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([0 15]);
ylim([0,0.25]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',24)
set(gcf, 'color', 'white');
hold on,
bar(nK4,pK4_ntotal(nK4+1),1,'r');
axis square;
set(gca,'linewidth', 2);
export_fig([figure_dir, 'probK4_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');

figure(7);
hold off;
histogram('BinEdges', -0.5+[0:ceil(ntotal/5)],'BinCounts', pK5_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([-0.5 10]);
ylim([0,0.35]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',24)
set(gcf, 'color', 'white');
hold on,
bar(nK5,pK5_ntotal(nK5+1),1,'r');
axis square;
set(gca,'linewidth', 2);
export_fig([figure_dir, 'probK5_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');

figure(8);
hold off;
histogram('BinEdges', -0.5+[0:ceil(ntotal/6)],'BinCounts', pK6_ntotal, 'FaceColor',[127,127,127]/256, 'FaceAlpha', 1);
xlim([-0.5 10]);
ylim([0,0.35]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',24)
set(gcf, 'color', 'white');
hold on,
bar(nK6,pK6_ntotal(nK6+1),1,'r');
axis square;
set(gca,'linewidth', 2);
export_fig([figure_dir, 'probK6_all_n=',num2str(ntotal), flat_prefix],'-m2', '-pdf');


%% OLD
figure(1)
hold off;
histogram('BinEdges', -0.5+[0:neggchambers+1],'BinCounts', pK1_neggchambers,'FaceColor' , [102,102,255]/255, 'FaceAlpha', 1);
xlim([0 neggchambers+1])
ylim([0,0.14])
%title(['K1 probability distribution for N = ',num2str(neggchambers),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);

%export_fig('probK1_rounded','-m2', '-pdf');


%figure(2)
%hold off;
histogram('BinEdges', -0.5+[0:nflatcysts+1],'BinCounts', pK1_flatcysts, 'FaceColor', 'r', 'FaceAlpha', 0.5);
xlim([0 nflatcysts+1]);
ylim([0,0.13]);
%title(['K1 probability distribution for N = ',num2str(nflatcysts),' and m = 72']);
set(gca,'FontSize',14)
set(findall(gcf,'type','text'),'FontSize',14)
hold on,
bar(14,pK1_neggchambers(15),'k','LineStyle','none');

bar(22,pK1_flatcysts(23),'k','LineStyle','none');
axis square;
set(gcf, 'color', 'white');
set(gca, 'FontSize', 28);
set(gca,'linewidth', 2);
axis square;
export_fig('probK1_flat','-m2', '-pdf');
print(gcf, '-dpdf', 'probK1_both.pdf');
%%
% Cramer van-Mises criterium
% Experiment:
N = ntotal;
oi = totalcounts;
% Simulations:
data_dir_paul='/Users/norbert/Dropbox (MIT)/StasPacking/Data/eqSpheres/processed_data_paul/';
data_dir_norbert='/Users/norbert/Dropbox (MIT)/StasPacking/Data/eqSpheres/processed_data_norbert/';
data = load([data_dir_paul, 'adjac_counts_prob.mat']);
descrcounts = data.Expression1(1,:);
adjaccounts = data.Expression1(2,:);
data = load([data_dir_norbert, 'adjac_counts_prob.mat']);
tmp1 = data.Expression1(1,:);
tmp2 = data.Expression1(2,:);
descrcounts = descrcounts + tmp1;
adjaccounts = adjaccounts + tmp2;
N = sum(descrcounts);
oi = descrcounts;


pi = 1/72*ones(1,72);
ei = N*pi;
Sj = cumsum(oi);
Tj = cumsum(ei);
Hj = Tj/N;
Zj = Sj - Tj;
Zb = sum(Zj.*pi);

W2 = 1/N*sum(Zj.^2.*pi);
U2 = 1/N*sum((Zj - Zb).^2.*pi);
tmp = Zj.^2.*pi./(Hj.*(1-Hj));
tmp(end) = 0;
A2 = 1/N*sum(tmp);

W2
U2
A2
