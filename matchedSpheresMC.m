% Packing using MC Metropolis updates
clear all;
cd '/home/norbert/StasPacking/SphericalConstraint2/'
load('meanradii.mat');
NRUNS=3000;

%Graph distance to oocyte:
GD(1)=0;
GD(2)=1; GD(3)=1; GD(5)=1; GD(9)=1;
GD(4)=2; GD(6)=2; GD(10)=2; GD(7)=2; GD(11)=2; GD(13)=2;
GD(8)=3; GD(12)=3; GD(14)=3; GD(15)=3;
GD(16)=4;

gdrad=[1.0,1,1,1,1];
gdrad=[1.0,1,1,1,1];
radii=[];
for i=1:16
    %radii(i) = gdrad(GD(i)+1);
    radii(i) = meanradii(i);
end

dr=2*mean(radii);

kcont = 1.5; %spherical container energy constant
%k=5.0; % Link force constant
k=5.0;
%ks=1.0; % sphere repulsion potential
ks=1.0;
%kv=0.0000; %volume energy constant
%ka=0.005; % area energy constant
kv=0;
ka=0;

Vspheres= 4*pi/3*sum(radii.^3); %The total volume of all cells
Vtarget = 0.8*Vspheres/0.74; % The target volume is the volume corrresponding to the densest packing of spheres.
Vtarget = 0;
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
%

visualize=false;
parfor RUN=1:NRUNS
    % Initial conditions:
    X0 = 8*max(radii).*(rand(16,3)-0.5); % X0 is now a 16x3 matrix with the position of each sphere center
    X0(:,3)=5*rand(16,1);
    X0(1,:)=[0,0,0];
    angles1=rand(4,1);
    
    X=X0;
    SAMPNO = char(java.util.UUID.randomUUID);
    %
    % Metropolis MC algorithm:
    % 1. Pick random point and random direction
    % 2. Calculate energy of new configuration
    % 3. Accept/reject move with detailed balance condition
    % 4. Iterate...
    
    contrad=[10,5,3,2.5,2,1.7,1.5,1.4];
    Nstepss=[700000,1000000,1000000,700000,500000,500000,500000,5000000];
    betas=[5,30,70,70,80,80,80,90,90];
    for rstep=1:numel(contrad)
        
        sphere_rad=dr*contrad(rstep);
        sphere_center=[0,0,(sphere_rad-radii(1))+0.5*radii(1)];
        SCMat = ones(15,3);
        SCMat(:,1)=sphere_center(1);
        SCMat(:,2)=sphere_center(2);
        SCMat(:,3)=sphere_center(3);
        
        Nsteps = Nstepss(rstep);
        %Nsteps = 1;
        MeanRad=mean(radii);
        Id=ones(1,3);
        beta=betas(rstep);
        
        D=L*X;
        D2=D.^2;
        Elinks = k/2.*sum((sqrt(sum(D2,2))-Drad).^2);
        D=S*X;
        D2=sum(D.^2,2);
        idx = find(D2-SradSq<0);
        Erepuls = ks/2.*sum((sqrt(D2(idx))-Srad(idx)).^2);
        % Zaxis:
        %Z=X(:,3);
        %zn=Z(Z<0);
        Ez=0;
        %Ez = k/2.*sum((zn).^2);
        
        %Volume:
        %[Vol, Area] = area3d(X(1:16,:));
        % The volume is augmented to include the finite thickness of the
        % points:
        %Vol = Vol + 0.05*Area*MeanRad;
        %Evol = kv/2.*(Vol-Vtarget)^2;
        Evol=0;
        %Earea = ka/2.*(Area)^2;
        Earea=0;
        %Spherical constraint:
        Econt=kcont/2.*sum((sqrt(sum((X(2:16,:)-SCMat).^2,2))-(sphere_rad*ones(15,1)-radii(2:16)')).^2);
        Etot = Elinks + Erepuls + Econt;
        Ntrials=0;
        Nacc = 0;
        if visualize
         figure(1);    
            draw_configuration(X,CU,radii,GD,sphere_rad, sphere_center);
        end
        Etots=zeros(1,Nsteps);
        Emin=Etot;
        Configurations=struct();
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
            id=randi(15)+1; % Don't pick the mother cell, this is our point of reference...
            Xtrial=X;
            %Move sphere id to a random position inside a cube of side length dr:
            drv = normrnd(0,dr/2,1,3);
            %drv = (2.0*rand(1,3)-Id)*dr;
            Xtrial(id,:) = X(id,:)+drv;
            %Calculate energy of new state:
            D=L*Xtrial;
            D2=D.^2;
            Elinks = k/2.*sum((sqrt(sum(D2,2))-Drad).^2);
            D=S*Xtrial;
            D2=sum(D.^2,2);
            idx = find(D2-SradSq<0);
            Erepuls = ks/2.*sum((sqrt(D2(idx))-Srad(idx)).^2);
            %Z=Xtrial(:,3);
            %zn=Z(Z<0);
            %Ez = k/2.*sum((zn).^2);
            Ez = 0;
            
            % Volume:
            %[Vol, Area] = area3d(Xtrial(1:16,:));
            % The volume is augmented to include the finite thickness of the
            % points:
            %Vol = Vol + 0.05*Area*MeanRad;
            %Evol = kv/2.*(Vol-Vtarget)^2;
            Evol=0;
            %Earea = ka/2.*(Area)^2;
            Earea=0;
            spdist = (sqrt(sum((Xtrial(2:16,:)-SCMat).^2,2))-(sphere_rad*ones(15,1)-radii(2:16)')).^2;
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
                    X= Xtrial;
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
                    figure(3);
                    hold off;
                    draw_configuration(XEmin,CU,radii,GD,sphere_rad, sphere_center);
                    xlim([-7,7]);
                    ylim([-7,7]);
                    zlim([-1,13]);
                    drawnow limitrate;
                end
            end
        end
    end
    State=struct();
    State.Configurations=Configurations;
    State.XEmin=XEmin;
    State.Emin = Emin;
    State.radii=radii;
    MCparsave(['Configuration_',SAMPNO,'.mat'],State);
end
