%% Packing using MC Metropolis updates --- Equal sphere packings with tree constraints
% Particle positions:
clear all;
cd '/home/norbert/StasPacking/eqSpheres/'
%cd '/Users/norbert/Projects/StasPacking/eqSpheres/'

radii=[];
for i=1:16
    radii(i) = 1;
end
dr=2*mean(radii);
NRUNS=4000;
% Cooling steps:
Nstepss=[700000,1000000,1000000,700000,500000,500000,5000000];
ddisps=[dr/5,dr/6,dr/8,dr/16,dr/20,dr/20,dr/30];
betas=[1,2,5,10,20,40,80];
%testing:
%Nstepss=[10,10,10,10,10,10,10];
%NRUNS=12;

visualize=false;

kcont = 1.5; %spherical container energy constant
k=5.0;
ks=1.0;

%Graph distance to oocyte:
GD(1)=0;
GD(2)=1; GD(3)=1; GD(5)=1; GD(9)=1;
GD(4)=2; GD(6)=2; GD(10)=2; GD(7)=2; GD(11)=2; GD(13)=2;
GD(8)=3; GD(12)=3; GD(14)=3; GD(15)=3;
GD(16)=4;
% Adjacency matrix:
C = zeros(16,16);
C(1,[2,3,5,9])=1;
C(2,[1,4,6,10])=1;
C(3,[1,7,11])=1;
C(4,[2,8,12])=1;
C(5,[1,13])=1;
C(6,[2,14])=1;
C(7,[3,15])=1;
C(8,[4,16])=1;
C(9,1)=1;
C(10,2)=1;
C(11,3)=1;
C(12,4)=1;
C(13,5)=1;
C(14,6)=1;
C(15,7)=1;
C(16,8)=1;

% Number of links:
Nlinks = nnz(C)/2;
%Spring connection matrix:
L=[];
CU=triu(C);
for i=1:16
   [a,idx]=find(CU(i,:));
   rows = zeros(numel(idx),16);
   rows(:,i)=1;
   for ii=1:numel(idx)
       rows(ii,idx(ii))=-1;
   end
   L=[L;rows];
end
for i=1:size(L,1)
    idx1 = find(L(i,:)==1);
    idx2 = find(L(i,:)==-1);
    Drad(i) = radii(idx1)+radii(idx2);
end
Drad=Drad';
% Hard-sphere repulsion matrix:
S=[];
CSU = triu(ones(16,16)-diag(ones(1,16)));
for i=1:16
   [a,idx]=find(CSU(i,:));
   rows = zeros(numel(idx),16);
   rows(:,i)=1;
   for ii=1:numel(idx)
       rows(ii,idx(ii))=-1;
   end
   S=[S;rows];
end
for i=1:size(S,1)
    idx1 = find(S(i,:)==1);
    idx2 = find(S(i,:)==-1);
    Srad(i) = radii(idx1)+radii(idx2);
end
Srad=Srad';
SradSq = Srad.^2;

sphere_rad=dr;
sphere_center=[0,0,0];
SCMat = ones(16,3);
SCMat(:,1)=sphere_center(1);
SCMat(:,2)=sphere_center(2);
SCMat(:,3)=sphere_center(3);

dd=dir('data');
if isempty(dd)
    disp('data directory does not exist! Exiting...')
end
data_dir='data/';

%%
parfor RUN=1:NRUNS
% Initial conditions:
X0 = 8*max(radii).*(rand(16,3)-0.5); % X0 is now a 16x3 matrix with the position of each sphere center
X0(:,3)=0;
X0(1,2)=0;
X0(1,1)=0;

X=X0;
SAMPNO = char(java.util.UUID.randomUUID);
%
% Metropolis MC algorithm with cooling
for rstep=1:numel(Nstepss)
    ddisp=ddisps(rstep);    
    Nsteps = Nstepss(rstep);
    Id=ones(1,3);
    beta=betas(rstep);
    
    D=L*X;
    D2=D.^2;
    Elinks = k/2.*sum((sqrt(sum(D2,2))-Drad).^2);
    D=S*X;
    D2=sum(D.^2,2);
    idx = find(D2-SradSq<0);
    Erepuls = ks/2.*sum((sqrt(D2(idx))-Srad(idx)).^2);

    %Spherical constraint:
    Econt=kcont/2.*sum((sqrt(sum((X-SCMat).^2,2))-sphere_rad*ones(16,1)).^2);
    Etot = Elinks + Erepuls + Econt;
    Ntrials=0;
    Nacc = 0;
    
    Etots=zeros(1,Nsteps);
    if rstep==1
        Emin=Etot;
        if visualize
          figure(1);
          draw_configuration(X,CU,0.5*radii,GD,sphere_rad, sphere_center);
        end
    end
    Configurations = struct();
    for i=1:1000
        Configurations(i).Energy = Emin;
        Configurations(i).X = X;
    end
    EnIdxMat = zeros(1000,2);
    EnIdxMat(:,2)=Emin;
    EnIdxMat(:,1)=[1:1000];
    SEnIdx = sortrows(EnIdxMat,2);
    ConfigEmax = Emin;
    for step=1:Nsteps
        id=randi(16); % Don't pick the mother cell, this is our point of reference...
        Xtrial=X;
        %Move sphere id to a random position inside a cube of side length dr:
        drv = normrnd(0,ddisp,1,3);
        if id==1
            drv(1)=0; drv(2)=0; % keep the central sphere at x=y=0
        end
        Xtrial(id,:) = X(id,:)+drv;
        %Calculate energy of new state:
        D=L*Xtrial;
        D2=D.^2;
        Elinks = k/2.*sum((sqrt(sum(D2,2))-Drad).^2);
        D=S*Xtrial;
        D2=sum(D.^2,2);
        idx = find(D2-SradSq<0);
        Erepuls = ks/2.*sum((sqrt(D2(idx))-Srad(idx)).^2);
              
        spdist = (sqrt(sum((Xtrial-SCMat).^2,2))-sphere_rad*ones(16,1)).^2;
        Econt=kcont/2*sum(spdist);
        acc=false;
        dE = Elinks + Erepuls + Econt - Etot;
        betadE = beta*dE;
        if betadE<75
            if betadE<0
                Etot = Elinks + Erepuls + Econt;
                X = Xtrial;
                Nacc = Nacc+1;
                acc=true;
            elseif exp(-betadE) > rand()
                Etot = Elinks + Erepuls + Econt;
                X = Xtrial;
                Nacc = Nacc+1;
                acc=true;
            end
        end
        Ntrials = Ntrials+1;
        draw=false;
        if Etot<Emin
            XEmin = X;
            Emin = Etot;
            draw=true;
        end;
        if Etot<ConfigEmax && acc==true
            ConfigEmaxIdx = SEnIdx(end,1);
            Configurations(ConfigEmaxIdx).Energy = Etot;
            Configurations(ConfigEmaxIdx).X = X;
            EnIdxMat(ConfigEmaxIdx,2) = Etot;
            SEnIdx = sortrows(EnIdxMat,2);
            ConfigEmaxIdx = SEnIdx(end,1);
            ConfigEmax = Configurations(ConfigEmaxIdx).Energy;
        end
        
        Etots(step)=Etot;
        if visualize
            drawnow limitrate;
        end
        if visualize && (mod(step,500)==0 || draw==true )
            Ntrials
            Nacc
            %figure(1);
            %scatter3(X(:,1),X(:,2),X(:,3));
            %drawnow;
            if draw
                figure(1);
                hold off;
                draw_configuration(XEmin,CU,0.5*radii,GD,sphere_rad, sphere_center);
                
                drawnow limitrate;
            end
        end
    end
end
    State=struct();
    State.XEmin = XEmin;
    State.Emin = Emin;
    State.Etots = Etots;
    State.Configurations = Configurations;
    MCparsave([data_dir, 'Configuration_',SAMPNO,'.mat'],State);
end
