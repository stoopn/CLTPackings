function [ rank1,rank2,rank3,rank4, name ] = extract_1d_descriptor( X, C, sphere_center )
% Well-defined computation of the 1d-descriptor from the actual 3d configuration by
%repeated tangent-plane projections. 
%Outline: Start with cell 0, construct best fitting plane through cell 0
%and its RC neighbors. Project neighbouring cells on this plane, extract
%orientation of neighbours => gives permutation of the first layer.
% Then do the same for all other cells until done.

%First layer - cells around the oocyte:
c=1;
nb=find(C(c,:)); % neighbour indices of cell c
points = [X(c,:); X(nb,:)]; % Position of neighbours and cell c
angles=tangent_projection(points,sphere_center);
[~,I]=sort(angles);
cycle=[nb(I), nb(I)];
idx=min(find(cycle==2));
rank1=cycle(idx:idx+3);

c=2;
nb=find(C(c,:)); % neighbour indices of cell c
points = [X(c,:); X(nb,:)]; % Position of neighbours and cell c
angles=tangent_projection(points,sphere_center);
[~,I]=sort(angles);
cycle=[nb(I), nb(I)];
idx=min(find(cycle==1));
rank2=cycle(idx+1:idx+3);

c=4;
nb=find(C(c,:)); % neighbour indices of cell c
points = [X(c,:); X(nb,:)]; % Position of neighbours and cell c
angles=tangent_projection(points,sphere_center);
[~,I]=sort(angles);
cycle=[nb(I), nb(I)];
idx=min(find(cycle==2));
rank3=cycle(idx+1:idx+2);

c=3;
nb=find(C(c,:)); % neighbour indices of cell c
points = [X(c,:); X(nb,:)]; % Position of neighbours and cell c
angles=tangent_projection(points,sphere_center);
[~,I]=sort(angles);
cycle=[nb(I), nb(I)];
idx=min(find(cycle==1));
rank4=cycle(idx+1:idx+2);

% Construct 1d descriptor name from ranks (from Paul's code):
name = [];
for j = 1:length(rank2),
    if rank2(j) == 4,
        for k = 1:length(rank3),
            if rank3(k) == 8,
                name = [name,16];
            else
                name = [name,rank3(k)];
            end
        end
    elseif rank2(j) == 6,
        name = [name,14];
    else
        name = [name,rank2(j)];
    end
end
for j = 2:(length(rank1)),
    if rank1(j) == 3,
        for k = 1:length(rank4),
            if rank4(k) == 7,
                name = [name,15];
            else
                name = [name,rank4(k)];
            end
        end
    elseif rank1(j) == 5,
        name = [name,13];
    elseif rank1(j) == 9,
        name = [name,9];
    end
end

end


function angles = tangent_projection(points,sphere_center) 
% Helper function. Returns the polar angles of the tangent plane projection of points,
% centered around points(1), i.e. the first point in the list of points.
[n,V,p] = affine_fit(points); % V is the basis of the plane.
% Construct orthogonal plane basis vectors:
a=V(:,1)/norm(V(:,1));
n=n/norm(n);
%Check orientation of normal:
if dot(n,points(1,:)-sphere_center)>0
   % disp('flipping normal');
    n=-n;
end
b=cross(n,a)/norm(cross(n,a));
V=[a,b];
T=inv(V'*V); %Not really necessary since we have already an orthogonal basis...
P= V*T*V';
projp=P*points'; %Projected point coordinates
% Obtain coordinates in terms of plane basis vectors
abprojp=[a';b']*projp;
%Subtract coordinate of central cell:
cc=repmat(abprojp(:,1),1,size(abprojp,2));
coord = abprojp - cc;
%Angles of neighbours of central cell:
angles=mod(atan2(coord(2,2:end),coord(1,2:end)),2*pi);
end
