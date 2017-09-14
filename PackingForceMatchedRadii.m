function [ FRES, Etot ] = PackingForceMatchedRadii( X, Params, calc_energy )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
if nargin<3
    calc_energy=false;
    Etot=-1;
end

DD=(Params.L*X);
LINKLEN=sqrt(sum(DD.^2,2));
LINKFORCE = - Params.k*(LINKLEN-Params.Drad).*DD./LINKLEN;

VFORCES=zeros(16,3);
VFORCES(:,1)=Params.L'*LINKFORCE(:,1);
VFORCES(:,2)=Params.L'*LINKFORCE(:,2);
VFORCES(:,3)=Params.L'*LINKFORCE(:,3);

DD=Params.S*X;
% D2=sum(DD.^2,2);
DISTANCES = sqrt(sum(DD.^2,2));
idx = find(DISTANCES-Params.Srad>=0);
REPULFORCE = - Params.ks*(DISTANCES-Params.Srad).*DD./DISTANCES;
REPULFORCE(idx,1)=0; REPULFORCE(idx,2)=0; REPULFORCE(idx,3)=0; % Set those pair forces to zero that are too far apart
RFORCES=zeros(16,3);
RFORCES(:,1)=Params.S'*REPULFORCE(:,1);
RFORCES(:,2)=Params.S'*REPULFORCE(:,2);
RFORCES(:,3)=Params.S'*REPULFORCE(:,3);

VECS = X-Params.SCMat;
DISTANCES = sqrt(sum(VECS.^2,2));
CONTFORCE = - Params.kcont*(DISTANCES - (Params.sphere_rad*ones(16,1)-Params.radii')).*VECS./DISTANCES;

FRES=RFORCES+VFORCES+CONTFORCE;
%FRES(1,:)=0; % Set force on cell 1 to zero;

if calc_energy
    D=Params.L*X;
    D2=D.^2;
    Elinks = Params.k/2.*sum((sqrt(sum(D2,2))-Params.Drad).^2);
    D=Params.S*X;
    D2=sum(D.^2,2);
    idx = find(D2-(Params.Srad).^2<0);
    Erepuls = Params.ks/2.*sum((sqrt(D2(idx))-Params.Srad(idx)).^2);

    %Spherical constraint (assume first sphere does not count, as in MC simulations:)
    Econt=Params.kcont/2.*sum((sqrt(sum((X(2:16,:)-Params.SCMat(2:16,:)).^2,2))-(Params.sphere_rad*ones(15,1)-Params.radii(2:16)')).^2);

    Etot = Elinks + Erepuls + Econt; 
    disp(Elinks);
    disp(Erepuls);
    disp(Econt);
end

end

