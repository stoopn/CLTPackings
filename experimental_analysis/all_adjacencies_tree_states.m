%% Reads all adjacency data from the different sets of rounded cyst experiments
clear all;
cd('/Users/norbert/Dropbox (MIT)/StasPacking/Data/experiment');

ncells = 16;
ringCanal = [[1,2];[1,3];[1,5];[1,9];...
    [2,4];[2,6];[2,10]; ...
    [3,7];[3,11];...
    [4,8];[4,12];...
    [5,13];...
    [6,14];...
    [7,15];...
    [8,16];...
    ];

%% get adjacencies of old rounded data set
load('ADJ_experiment_rounded_data1.mat');

%% Get adjacencies from extra data set 1: Run the part above first because here we append the data!
neggchambers = 45;
n = 1;
ctr = size(matrix,2);

for j = 1:neggchambers
    if j<26
        xlrange = sprintf('%c1:%c19', char('A'+j), char('A'+j));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra1', n, xlrange);
    else
        xlrange = sprintf('A%c1:A%c19', char('A'+(j-26)), char('A'+(j-26)));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra1', n, xlrange);
    end
    
    matrixtemp = zeros(16);
    for i = 5:19
        text=raw{i};
        
        num = regexprep(text, '[', '');
        num = regexprep(num, ']', '');
        num = regexprep(num, '?', '');
        num = sscanf(num, '%d')';
        %disp(i)
        
        for k = 1:length(num)
            if (num(k) ~= i-3)
                matrixtemp(i-3, num(k)) = 1;
                matrixtemp(num(k), i-3) = 1;
            end
        end
    end
    
    for k = 1:length(ringCanal),
        matrixtemp(ringCanal(k,1),ringCanal(k,2)) = 1;
        matrixtemp(ringCanal(k,2),ringCanal(k,1)) = 1;
    end
    matrix{ctr+j} = matrixtemp;
end

%% Do not use - just for n=80 partial statistics
% Rebuild matrixfreq:
matrixfreq = zeros(16);
for i=1:size(matrix,2)
    matrixfreq = matrixfreq + matrix{i};
end
% Normalize matrixfreq:
matrixfreq = matrixfreq/size(matrix,2);

neggchambers = size(matrix,2);
save(['ADJ_experiments_rounded_n=', num2str(neggchambers), '.mat'],'matrix', 'matrixfreq','neggchambers');
%% Now extra data 2:
%  Run the part above first because here we append the data!
neggchambers = 21;
n = 1;
ctr = size(matrix,2);

for j = 1:neggchambers
    if j<26
        xlrange = sprintf('%c1:%c19', char('A'+j), char('A'+j));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra2', n, xlrange);
    else
        xlrange = sprintf('A%c1:A%c19', char('A'+(j-26)), char('A'+(j-26)));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra2', n, xlrange);
    end
    
    matrixtemp = zeros(16);
    for i = 5:19
        text=raw{i};
        
        num = regexprep(text, '[', '');
        num = regexprep(num, ']', '');
        num = regexprep(num, '?', '');
        num = sscanf(num, '%d')';
        %disp(i)
        
        for k = 1:length(num)
            if (num(k) ~= i-3)
                matrixtemp(i-3, num(k)) = 1;
                matrixtemp(num(k), i-3) = 1;
            end
        end
    end
    
    for k = 1:length(ringCanal)
        matrixtemp(ringCanal(k,1),ringCanal(k,2)) = 1;
        matrixtemp(ringCanal(k,2),ringCanal(k,1)) = 1;
    end
    matrix{ctr+j} = matrixtemp;
end

%% Now extra data 3:
%  Run the part above first because here we append the data!
neggchambers = 20;
n = 1;
ctr = size(matrix,2);

for j = 1:neggchambers
    if j<26
        xlrange = sprintf('%c1:%c19', char('A'+j), char('A'+j));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra3', n, xlrange);
    else
        xlrange = sprintf('A%c1:A%c19', char('A'+(j-26)), char('A'+(j-26)));
        [~, column, raw] = xlsread('rounded_cysts_Data_extra3', n, xlrange);
    end
    
    matrixtemp = zeros(16);
    for i = 5:19
        text=raw{i};
        
        num = regexprep(text, '[', '');
        num = regexprep(num, ']', '');
        num = regexprep(num, '?', '');
        num = sscanf(num, '%d')';
        %disp(i)
        
        for k = 1:length(num)
            if (num(k) ~= i-3)
                matrixtemp(i-3, num(k)) = 1;
                matrixtemp(num(k), i-3) = 1;
            end
        end
    end
    
    for k = 1:length(ringCanal)
        matrixtemp(ringCanal(k,1),ringCanal(k,2)) = 1;
        matrixtemp(ringCanal(k,2),ringCanal(k,1)) = 1;
    end
    matrix{ctr+j} = matrixtemp;
end
%%
% Rebuild matrixfreq:
matrixfreq = zeros(16);
for i=1:size(matrix,2)
    matrixfreq = matrixfreq + matrix{i};
end
% Normalize matrixfreq:
matrixfreq = matrixfreq/size(matrix,2);

%%
neggchambers = size(matrix,2);
save('ADJ_experiments_rounded_all.mat','matrix', 'matrixfreq','neggchambers');

%% Now do all the 1d descriptors:
% For data1 (the originzl rounded cysts), we automatically generate 1d
% descriptors in rounded_cysts_data1.m

% Load all 72 tree state names and save them in ordered form for easier identification: 
load('../names_all.mat');
names_all_ordered = [];
for i=1:size(names_all,1)
    [val,idx]=find(names_all(i,:)==16)
    names_all_ordered(i,:) = circshift(names_all(i,:),-(idx-1));
end
save('../names_all_ordered.mat','names_all_ordered');

% First and second set of tree state data:
load('names_rounded_data1.mat');
ExpNames=names;
ExpNames_ordered=[];
ctr=size(names,2);

extra1=csvread('1D_descriptors_extra1.csv');
for i=1:size(extra1,1)
    ExpNames{ctr+i} = extra1(i,:); 
end

% Extra data 2:
ctr=size(ExpNames,2);
extra1=csvread('1D_descriptors_extra2.csv');
for i=1:size(extra1,1)
    ExpNames{ctr+i} = extra1(i,:); 
end

%Extra data 3:
ctr=size(ExpNames,2);
extra1=csvread('1D_descriptors_extra3.csv');
for i=1:size(extra1,1)
    ExpNames{ctr+i} = extra1(i,:); 
end

% Convert to state IDs
for i=1:size(ExpNames,2)
    [val,idx]=find(ExpNames{i}==16);
    ExpNames_ordered(i,:) = circshift(ExpNames{i},-(idx-1));
    intemp = find(ismember(names_all_ordered,ExpNames_ordered(i,:),'rows'));
    ExpStates(i) = intemp;
end
save('ExpStates_rounded_all.mat','ExpStates', 'ExpNames_ordered');
