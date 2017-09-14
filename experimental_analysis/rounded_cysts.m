

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


%% what is the 1d descriptor of the network of cellular contacts

clc,
names = {};
for i = 1:neggchambers,
    c = 1;
    % cells that are one edge away from the oocyte
    cellsrank1 = [3,5,9];
    if length(find(matrix{i}(2,cellsrank1))) == 2,
        indrank1 = find(matrix{i}(2,cellsrank1));
        rank1 = [];% ordered cells that are one edge away from the oocyte
        rank1 = [2,cellsrank1(indrank1(1)),cellsrank1(find([1,2,3]~=indrank1(1) & [1,2,3]~=indrank1(2))),cellsrank1(indrank1(2))];
        rank1,
    end
    cellsrank2 = [4,6,10];
    rank2 = []; % ordered cells descendants of cell 2
    if length(find(matrix{i}(rank1(4),cellsrank2))) == 1,
        if length(find(matrix{i}(rank1(2),cellsrank2))) == 1,
            idx1 = find(matrix{i}(rank1(4),cellsrank2));
            idx3 = find(matrix{i}(rank1(2),cellsrank2));
            rank2 = [cellsrank2(idx1),cellsrank2(find([1,2,3] ~= idx1 & [1,2,3] ~= idx3)),cellsrank2(idx3)];
            rank2,
        else
            c = 0;
            i,
        end
    else
        c = 0;
        i,
    end
    cellsrank3 = [8,12];
    rank3 = [];
    idx1 = find(rank2 == 4);
    if idx1 == 1,
        if length(find(matrix{i}(rank1(4), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank1(4), cellsrank3));
            rank3 = [cellsrank3(idx2), cellsrank3(find([1,2] ~= idx2))],
        elseif length(find(matrix{i}(rank2(2), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank2(2), cellsrank3));
            rank3 = [cellsrank3(find([1,2] ~= idx2)), cellsrank3(idx2)],
        else
            c = 0;
            i
        end
    elseif idx1 == 2,
        if length(find(matrix{i}(rank2(1), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank2(1), cellsrank3));
            rank3 = [cellsrank3(idx2), cellsrank3(find([1,2] ~= idx2))],
        elseif length(find(matrix{i}(rank2(3), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank2(3), cellsrank3));
            rank3 = [cellsrank3(find([1,2] ~= idx2)), cellsrank3(idx2)],
        else
            c = 0;
            i,
        end
    elseif idx1 == 3,
        if length(find(matrix{i}(rank1(2), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank1(2), cellsrank3));
            rank3 = [cellsrank3(find([1,2] ~= idx2)), cellsrank3(idx2)],
        elseif length(find(matrix{i}(rank2(2), cellsrank3))) == 1,
            idx2 = find(matrix{i}(rank2(2), cellsrank3));
            rank3 = [cellsrank3(idx2), cellsrank3(find([1,2] ~= idx2))],
        else
            c = 0;
            i,
        end
    end
    
    cellsrank4 = [7,11];
    rank4 = [];
    idx1 = find(rank1 == 3);
    if length(find(matrix{i}(rank1(idx1-1),cellsrank4))) == 1,
        idx2 = find(matrix{i}(rank1(idx1-1),cellsrank4));
        rank4 = [cellsrank4(idx2), cellsrank4(find([1,2] ~= idx2))],
        
    elseif length(find(matrix{i}(rank1(idx1+1),cellsrank4))) == 1,
        idx2 = find(matrix{i}(rank1(idx1+1),cellsrank4));
        rank4 = [cellsrank4(find([1,2] ~= idx2)),cellsrank4(idx2)],
        
    else
        c = 0;
        i;
    end
    i,c;
    if c,
        names{i} = [];
        for j = 1:length(rank2),
            if rank2(j) == 4,
                for k = 1:length(rank3),
                    if rank3(k) == 8,
                        names{i} = [names{i},16];
                    else
                        names{i} = [names{i},rank3(k)];
                    end
                end
            elseif rank2(j) == 6,
                names{i} = [names{i},14];
            else
                names{i} = [names{i},rank2(j)];
            end
        end
        for j = 2:(length(rank1)),
            if rank1(j) == 3,
                for k = 1:length(rank4),
                    if rank4(k) == 7,
                        names{i} = [names{i},15];
                    else
                        names{i} = [names{i},rank4(k)];
                    end
                end
            elseif rank1(j) == 5,
                names{i} = [names{i},13];
            elseif rank1(j) == 9,
                names{i} = [names{i},9];
            end
        end
    end
    
    if i == 7 & length(names) < i,
        names{i} = [16,12,10,14,13,15,11,9];
    elseif i == 15 & length(names) < i,
        names{i} = [14,12,16,10,11,15,13,9];
    elseif i == 17 & length(names) < i,
        names{i} = [14,16,12,10,15,11,13,9];
    elseif i == 19 & length(names) < i,
        names{i} = [14,10,12,16,11,15,13,9];
    elseif i == 22 & length(names) < i,
        names{i} = [12,16,10,14,13,15,11,9];
    elseif i == 23 ,
        names{i} = [12, 16, 10, 14, 11, 15, 9, 13];
    elseif i == 24 & length(names) < i,
        names{i} = [14,16,12,10,15,11,13,9];
    elseif i == 25 & length(names) < i,
        names{i} = [10,14,16,12,11,15,13,9];
    elseif i == 28 & length(names) < i,
        names{i} = [12,16,10,14,13,15,11,9];
    elseif i == 30 & length(names) < i,
        names{i} = [14,10,12,16,11,15,9,13];
    elseif i == 32 & length(names) < i,
        names{i} = [14, 12, 16, 10, 11, 15, 9, 13];
    end
    
end


%% ordering egg chamber with lexicographical order on 1d descriptor

orderedeggchambers = 1:neggchambers;

swapped = 1;
while swapped
    swapped = 0;
    for i = 2:neggchambers,
        if lexico_order(names{orderedeggchambers(i-1)},names{orderedeggchambers(i)}) == 1,
            idxtemp = orderedeggchambers(i-1);
            orderedeggchambers(i-1) = orderedeggchambers(i);
            orderedeggchambers(i) = idxtemp;
            swapped = 1;
        end
    end
end

for i = 1:neggchambers,
    names{orderedeggchambers(i)},
end



