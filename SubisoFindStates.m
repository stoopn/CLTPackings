% Here we import the Mathematica coordinates of all subisomoprhism of the
% RC tree into the J90 solid in order to calculate the corresponding 1d
% descriptor and hence the tree state. We save the states at the end for
% later analysis.
clear all;
cd('/Users/norbert/Dropbox (MIT)/StasPacking/Data/Thomson/')
data = load('T0SubisoCoords.mat');
RCdata = load('../RCAdj.mat');
load('../names_all.mat');

sphere_center=[0,0,0];
sphere_radius=1;

expnames = fieldnames(data);

states = [];
names={};

for m=1:numel(expnames)
    if mod(m, 100)==0
        disp(m)
    end
   X =  getfield(data, expnames{m});
   [ rank1,rank2,rank3,rank4, name ] = extract_1d_descriptor( X, RCdata.C, sphere_center );
    names{m}=name; 
    state=find(ismember(names_all,name,'rows'));
    if isempty(state)
       % Reverse name and check again:
       nameR = fliplr(name);
       idx = min(find(nameR==14 | nameR==10 | nameR==12 | nameR==16));
       nameT = nameR(mod([idx-1:idx+6],8)+1);
       state=find(ismember(names_all,nameT,'rows'));
       names{m}=names_all(state,:);
    end
    states(m) = state;
end
save('T0SubisoStates.mat', 'states');
%%
figure(1)
histogram(states,[0.5:72.5]);

%% Repeat the same for the T1 Thomson graph
data = load('T1SubisoCoords.mat');
RCdata = load('../RCAdj.mat');
load('../names_all.mat');

sphere_center=[0,0,0];
sphere_radius=1;

expnames = fieldnames(data);

states = [];
names={};

for m=1:numel(expnames)
    if mod(m, 100)==0
        disp(m)
    end
   X =  getfield(data, expnames{m});
   [ rank1,rank2,rank3,rank4, name ] = extract_1d_descriptor( X, RCdata.C, sphere_center );
    names{m}=name; 
    state=find(ismember(names_all,name,'rows'));
    if isempty(state)
       % Reverse name and check again:
       nameR = fliplr(name);
       idx = min(find(nameR==14 | nameR==10 | nameR==12 | nameR==16));
       nameT = nameR(mod([idx-1:idx+6],8)+1);
       state=find(ismember(names_all,nameT,'rows'));
       names{m}=names_all(state,:);
    end
    states(m) = state;
end
save('T1SubisoStates.mat', 'states');
%%
figure(2)
histogram(states,[0.5:72.5]);
