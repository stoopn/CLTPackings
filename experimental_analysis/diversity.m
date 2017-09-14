presets;

%% get adjacencies
neggchambers = 35;
ncells = 16;
n = 37  ;
ringCanal = [[1,2];[1,3];[1,5];[1,9];...
    [2,4];[2,6];[2,10]; ...
    [3,7];[3,11];...
    [4,8];[4,12];...
    [5,13];...
    [6,14];...
    [7,15];...
    [8,16];...
    ];
% matrixRingCanal = zeros(16); for i = 1:size(ringCanal,1),matrixRingCanal(ringCanal(i,1),ringCanal(i,2)) = 1;matrixRingCanal(ringCanal(i,2),ringCanal(i,1)) = 1; end
% [V,D] = eig(matrixRingCanal);

matrixfreq = zeros(16);
matrix = {};
for j = 1:neggchambers,
    if j<26,
        [~, text, ~] = xlsread('7-7-15 Ring Canal Data.xlsx', n, [char('A'+j),num2str(3)])
    else
        [~, text, ~] = xlsread('7-7-15 Ring Canal Data.xlsx', n, ['A',char('A'+(j-26)),num2str(3)])
    end
    matrixtemp = zeros(16);
    for i = 5:19,
        if j < 26,
            [~, text, ~] = xlsread('7-7-15 Ring Canal Data.xlsx', n, [char('A'+j),num2str(i)]);
        else
            [~, text, ~] = xlsread('7-7-15 Ring Canal Data.xlsx', n, ['A',char('A'+(j-26)),num2str(i)]);
        end
        num = regexprep(text, '[', '');
        num = regexprep(num, ']', '');
        num = regexprep(num, '?', '');
        num = sscanf(num{1}, '%d')';
        %disp(i)
        
        for k = 1:length(num)
            if (num(k) ~= i-3),
                matrixtemp(i-3, num(k)) = 1;
                matrixtemp(num(k), i-3) = 1;
            end
        end
    end
    
    for k = 1:length(ringCanal),
        matrixtemp(ringCanal(k,1),ringCanal(k,2)) = 1;
        matrixtemp(ringCanal(k,2),ringCanal(k,1)) = 1;
    end
    matrix{j} = matrixtemp;
    matrixfreq = matrixfreq + matrix{j};
    
end
matrixfreq = matrixfreq./(neggchambers);


figure, imagesc(matrixfreq); colorbar
clear text;

%%

diff = zeros(neggchambers);
diff_all = cell(neggchambers);
diffs = [];
indices = [];
figure,
for j = 1:neggchambers-1,
    imagesc(matrix{j});
    pause(0.01);
    for k = j+1:neggchambers,
        diff(j,k) = sum(abs(matrix{j}(:) - matrix{k}(:)));
        diff_all{j,k} = matrix{j}(:) - matrix{k}(:);
        diffs = [diffs;sum(abs(matrix{j}(:) - matrix{k}(:)))];
        if sum(abs(matrix{j}(:) - matrix{k}(:))) == 0,
            indices = [indices;[j,k]];
        end
    end
end

diff = diff + diff';

%% histogram of differences between egg chambers adjacencies


figure,

hist(diffs);
xlabel('Pairwise differences')
ylabel('Number of pairs');

