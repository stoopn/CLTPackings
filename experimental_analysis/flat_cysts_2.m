% opening 1d descriptor obtained by hand

nflatcysts = 46;
% 
names_hand = {};
opposite_dir = [4,7,8,9,10,11,12,13,15,16,17,22,23,24,27,30,32,34,36,37,38,40,41,43,45];% correspond to flat cysts with 1d descriptor written in the "counterclockwise sense"

for j = 1:nflatcysts,
    i = j+1;
    [NUM,TXT,RAW] = xlsread('FLAT_cysts_descriptors4.xlsx',['C',num2str(i),':J',num2str(i)]);
    if length(find(ismember(opposite_dir,j))),
        % read name on the opposite direction
        names_flat{j} = [NUM(4:-1:1),NUM(end:-1:5)];
    else
        names_flat{j} = NUM;
    end
end


%% ordering flat cysts 1d descriptor using lexicographical order

orderedflatcysts = 1:nflatcysts;

swapped = 1;
while swapped,
    swapped = 0;
    for i = 2:nflatcysts,
        if lexico_order(names_flat{orderedflatcysts(i-1)},names_flat{orderedflatcysts(i)}) == 1,
            idxtemp = orderedflatcysts(i-1);
            orderedflatcysts(i-1) = orderedflatcysts(i);
            orderedflatcysts(i) = idxtemp;
            swapped = 1;
        end
    end
end

for i = 1:nflatcysts,
    names_flat{orderedflatcysts(i)},
end

