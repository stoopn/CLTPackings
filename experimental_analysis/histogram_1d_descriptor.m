% for better visualization settings
presets

%%
% list of all 1d descriptors

names_all = [
    [10,12,16,14,11,15,9,13];
    [10,12,16,14,11,15,13,9];
    [10,12,16,14,13,11,15,9];
    [10,12,16,14,13,15,11,9];
    [10,12,16,14,15,11,9,13];
    [10,12,16,14,15,11,13,9];
    %
    [10,14,12,16,11,15,9,13];
    [10,14,12,16,11,15,13,9];
    [10,14,12,16,13,11,15,9];
    [10,14,12,16,13,15,11,9];
    [10,14,12,16,15,11,9,13];
    [10,14,12,16,15,11,13,9];
    %
    [10,14,16,12,11,15,9,13];
    [10,14,16,12,11,15,13,9];
    [10,14,16,12,13,11,15,9];
    [10,14,16,12,13,15,11,9];
    [10,14,16,12,15,11,9,13];
    [10,14,16,12,15,11,13,9];
    %
    [10,16,12,14,11,15,9,13];
    [10,16,12,14,11,15,13,9];
    [10,16,12,14,13,11,15,9];
    [10,16,12,14,13,15,11,9];
    [10,16,12,14,15,11,9,13];
    [10,16,12,14,15,11,13,9];
    %
    [12,16,10,14,11,15,9,13];
    [12,16,10,14,11,15,13,9];
    [12,16,10,14,13,11,15,9];
    [12,16,10,14,13,15,11,9];
    [12,16,10,14,15,11,9,13];
    [12,16,10,14,15,11,13,9];
    %
    [12,16,14,10,11,15,9,13];
    [12,16,14,10,11,15,13,9];
    [12,16,14,10,13,11,15,9];
    [12,16,14,10,13,15,11,9];
    [12,16,14,10,15,11,9,13];
    [12,16,14,10,15,11,13,9];
    %
    [14,10,12,16,11,15,9,13];
    [14,10,12,16,11,15,13,9];
    [14,10,12,16,13,11,15,9];
    [14,10,12,16,13,15,11,9];
    [14,10,12,16,15,11,9,13];
    [14,10,12,16,15,11,13,9];
    %
    [14,10,16,12,11,15,9,13];
    [14,10,16,12,11,15,13,9];
    [14,10,16,12,13,11,15,9];
    [14,10,16,12,13,15,11,9];
    [14,10,16,12,15,11,9,13];
    [14,10,16,12,15,11,13,9];
    %
    [14,12,16,10,11,15,9,13];
    [14,12,16,10,11,15,13,9];
    [14,12,16,10,13,11,15,9];
    [14,12,16,10,13,15,11,9];
    [14,12,16,10,15,11,9,13];
    [14,12,16,10,15,11,13,9];
    %
    [14,16,12,10,11,15,9,13];
    [14,16,12,10,11,15,13,9];
    [14,16,12,10,13,11,15,9];
    [14,16,12,10,13,15,11,9];
    [14,16,12,10,15,11,9,13];
    [14,16,12,10,15,11,13,9];
    %
    [16,12,10,14,11,15,9,13];
    [16,12,10,14,11,15,13,9];
    [16,12,10,14,13,11,15,9];
    [16,12,10,14,13,15,11,9];
    [16,12,10,14,15,11,9,13];
    [16,12,10,14,15,11,13,9];
    %
    [16,12,14,10,11,15,9,13];
    [16,12,14,10,11,15,13,9];
    [16,12,14,10,13,11,15,9];
    [16,12,14,10,13,15,11,9];
    [16,12,14,10,15,11,9,13];
    [16,12,14,10,15,11,13,9];
    ];


%% hist flat cysts
ind_names_flat = [];
res_flat = [1];
res_flat_names = [names_flat{orderedflatcysts(1)}];


intemp = find(ismember(names_all,names_flat{orderedflatcysts(1)},'rows'));

if ~isempty(intemp),
    ind_names_flat = intemp;
else
    ind_names_flat = 0;
end

for i = 2:nflatcysts,
    if eq_seq(names_flat{orderedflatcysts(i-1)},names_flat{orderedflatcysts(i)}) == 1,
        
        res_flat(end) = res_flat(end)+1;
    else
        % indice of the 1d descriptor in the list of all names
        intemp = find(ismember(names_all,names_flat{orderedflatcysts(i)},'rows'));
        
        if ~isempty(intemp),
            ind_names_flat = [ind_names_flat,intemp];
        else
            ind_names_flat = [ind_names_flat,0];
            orderedflatcysts(i)
        end
        
        res_flat = [res_flat,1];
        res_flat_names = [res_flat_names;names_flat{orderedflatcysts(i)}];
    end
end
res_flat;


%% hist rounded cysts

ind_names = [];
intemp = find(ismember(names_all,names{orderedeggchambers(1)},'rows'));

if ~isempty(intemp),
    ind_names = intemp;
else
    ind_names = 0;
end



res = [1];
res_names = [[10,12,16,14,13,11,15,9]];
for i = 2:neggchambers,
    if eq_seq(names{orderedeggchambers(i-1)},names{orderedeggchambers(i)}) == 1,
        i,
        res(end) = res(end)+1;
    else
        res = [res,1];
        res_names = [res_names;names{orderedeggchambers(i)}];
        % indice of the 1d descriptor in the list of all names
        intemp = find(ismember(names_all,names{orderedeggchambers(i)},'rows'));
        if ~isempty(intemp),
            ind_names = [ind_names,intemp];
        else
            ind_names = [ind_names,0];
            orderedeggchambers(i)
        end
    end
end

%% histogram of 1d descriptor

figure(5),

subplot(2,1,1)
pos_flat = [];

bar(ind_names_flat,res_flat(1:end))

xlim([0 73]);
ylim([0 4]);
ax = gca;
ax.YTick = [0 1 2 3 4];
ax.XTick = [1:72];

%set(findall(gcf,'type','text'),'FontSize',14)
title([num2str(nflatcysts),' Flat Cysts - E(K1) = ',num2str(nflatcysts*((71/72)^(nflatcysts-1))),' - K1 =', num2str(length(find(res_flat == 1))),' - T = ',num2str(nflatcysts*((71/72)^(nflatcysts-1)) - length(find(res_flat == 1))), ' - T_{\alpha} = ',num2str(nflatcysts*sqrt(log(5))/(2*sqrt(72)))]);
% formula for expectation of K1 is from Liam Paninski 2008

set(gca,'fontsize', 16);


%late stage cysts
subplot(2,1,2)
bar(ind_names,res(1:end)),

ylim([0 4]);
xlim([0 73]);
ax = gca;

ax.YTick = [0 1 2 3 4]
ax.XTick = [1:72]


title([num2str(neggchambers),' Round Cysts - E(K1) = ',num2str(neggchambers*((71/72)^(neggchambers-1))),' - K1 =', num2str(length(find(res == 1))),' - T = ',num2str(neggchambers*((71/72)^(neggchambers-1)) - length(find(res == 1))), ' - T_{\alpha} = ',num2str(neggchambers*sqrt(log(5))/(2*sqrt(72)))]);
% formula for expectation of K1 is from Liam Paninski 2008
set(gca,'fontsize', 16);
%early flat

set(gca, 'ydir', 'reverse');

