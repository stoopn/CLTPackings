clear all;
%On cluster:
cd('/home/stoopn/Dropbox (MIT)/StasPacking/Data/matchedRadii/processed_data');

%cd('/Users/norbert/Dropbox (MIT)/StasPacking/Data/matchedRadii/processed_data');
%
input_data = load('processed_data.mat');
mycolors=[133, 142, 141; 55, 144, 166; 194, 45, 45; 0, 107, 72; 252, 183, 21]/256;

%%
% debugging:
  contrad=[10,5,3,2.5,2,1.7,1.5,1.4];
  dr=2*mean(Params.radii)
sphere_rad=dr*contrad(8);
        sphere_center=[0,0,(sphere_rad-Params.radii(1))+0.5*Params.radii(1)];
        SCMat = ones(16,3);
        SCMat(:,1)=sphere_center(1);
        SCMat(:,2)=sphere_center(2);
        SCMat(:,3)=sphere_center(3);
        Params.SCMat = SCMat;
        Params.sphere_rad=sphere_rad;
        Params.sphere_center = sphere_center;
%%
EquilData=struct();
EquilData.Params = input_data.Params;
EquilData.Params.Orig_sphere_rad = input_data.Params.sphere_rad;

Data=struct();
IData=input_data.Data;
Params = input_data.Params;
% Fix last column of SCMat:
SCMat = Params.SCMat;
SCMat = [SCMat; SCMat(1,:)];
EquilData.Params.SCMat = SCMat;
Params.SCMat = SCMat;
CU=triu(Params.C);

 Nsteps=3000;
 dt=0.05;
 visualize=false;
parfor sample=1:input_data.Params.NSamples
 X=IData(sample).X;
 X0=X;
 disp(sample);
 for step=1:Nsteps
     [force, Etot] = PackingForceMatchedRadii(X,Params);
     Xnew = X+dt*force;
     X= Xnew;
     if mod(step,100)==0 && visualize
     figure(1);
     hold off;
     draw_configuration(X,CU,0.5*Params.radii,Params.GD,Params.sphere_rad, Params.sphere_center);
     drawnow;
     end
 end
 % Calculate final energy:
 [force, Etot] = PackingForceMatchedRadii(X,Params, true);
 
 Data(sample).filename = IData(sample).filename;
 Data(sample).X = X;
 Data(sample).OldX = X0;
 Data(sample).Energy = Etot;
 Data(sample).OldEnergy = IData(sample).MinEnergy;
 
 figure(2);
      draw_configuration(X0,CU,0.5*Params.radii,Params.GD,Params.sphere_rad, Params.sphere_center);
end
EquilData.Data = Data;
save('equil_data.mat', 'EquilData');

%%
% Extract 1d descriptor state and convex hull adjacency matrix:
%load('equil_data.mat');
load('equil_data.mat');

load('names_all.mat');

CHAdjMatrix={};
AdjMatrix={};

states = [];
names={};

Data=EquilData.Data;
Params=EquilData.Params;
tol=0.1;
parfor sample=1:Params.NSamples
    X=Data(sample).X;
    % Construct convex hull adjacency matrix:
    K=convhull(X(:,1),X(:,2),X(:,3),'simplify', true);
    TR = triangulation(K, X(:,1), X(:,2), X(:,3));
    TADJ=zeros(16,16);
    idx=TR.edges;
    for i=1:size(idx,1)
       j=idx(i,2);
        i=idx(i,1);
      TADJ(i,j)=1;
      TADJ(j,i)=1;
    end
    CHAdjMatrix{sample}=TADJ;
    % Construct contact adjacency matrices:
    D=Params.S*X;
    D2=sum(D.^2,2);
    idx = find(D2-(1+tol)*(Params.Srad).^2<0);
    ADJ = zeros(16,16);
    for ii=1:numel(idx)
        idd = find(Params.S(idx(ii),:));
        ADJ(idd(1),idd(2))=1;
    end
    ADJ=ADJ+ADJ';
    ADJ=bitor(ADJ,Params.C);
    AdjMatrix{sample}=ADJ;
    
    % Extract 1d descriptor:
    [ rank1,rank2,rank3,rank4, name ] = extract_1d_descriptor( X, Params.C, Params.sphere_center );
    names{sample}=name; 
    state=find(ismember(names_all,name,'rows'));
    if isempty(state)
       % Reverse name and check again:
       nameR = fliplr(name);
       idx = min(find(nameR==14 | nameR==10 | nameR==12 | nameR==16));
       nameT = nameR(mod([idx-1:idx+6],8)+1);
       state=find(ismember(names_all,nameT,'rows'));
       names{sample}=names_all(state,:);
    end
    states(sample) = state;
end

save('Equil_CHAdjMatrix.mat', 'CHAdjMatrix');
save(['Equil_AdjMatrix_tol=',num2str(tol),'.mat'], 'AdjMatrix');
save('Equil_states.mat','states');

%%
load('equil_data.mat');
Data=EquilData.Data;
Params=EquilData.Params;
Energy=[];
for sample=1:Params.NSamples
    Energy(sample) = Data(sample).Energy;
end
save('Equil_Energy.mat', 'Energy');

%%
% Extract convex hull adjacency matrix by projecting:
load('equil_data.mat');

PrCHAdjMatrix={};
Data=EquilData.Data;
Params=EquilData.Params;
parfor sample=1:Params.NSamples
    X=Data(sample).X;
    Xs=X-repmat(Params.sphere_center,16,1);
    Xs=Xs./sqrt(sum(Xs.^2,2));
    % Construct convex hull adjacency matrix:
    K=convhull(Xs(:,1),Xs(:,2),Xs(:,3),'simplify', true);
    TR = triangulation(K, Xs(:,1), Xs(:,2), Xs(:,3));
    TADJ=zeros(16,16);
    idx=TR.edges;
    for i=1:size(idx,1)
       j=idx(i,2);
        i=idx(i,1);
      TADJ(i,j)=1;
      TADJ(j,i)=1;
    end
    PrCHAdjMatrix{sample}=TADJ;
end
save('Equil_PrCHAdjMatrix.mat', 'PrCHAdjMatrix');

%%
% Save so mathematica can read the coordinates:
load('equil_data.mat');
Xs = {EquilData.Data.X};
save('equil_coords_mathematica.mat', '-v6', 'Xs')