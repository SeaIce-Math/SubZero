function [Floe] = corners(Floe,Nb,Floe0,c2_boundary)
%Observations of older floe fields show a tendency to form rounder shapes through 
%repeated interactions with other floes. The corner grinding mechanism uses the 
%contact overlap areas to determine whether a floe could have its corner fractured
id ='MATLAB:polyshape:repairedBySimplify';
warning('off',id)
id3 ='MATLAB:polyshape:boundary3Points';
warning('off',id3)

Lx= max(c2_boundary(1,:));
Ly= max(c2_boundary(2,:));

%Create ghost floes
ghostFloeX=[];
ghostFloeY=[];

x=cat(1,Floe0.Xi);
y=cat(1,Floe0.Yi);
alive=cat(1,Floe0.alive);

for i=1:length(Floe0)
    poly = polyshape(Floe0(i).c_alpha'+[x(i) y(i)]);
    if alive(i) && (max(abs(poly.Vertices(:,1)))>Lx)
        
        ghostFloeX=[ghostFloeX  Floe0(i)];
        ghostFloeX(end).Xi=Floe0(i).Xi-2*Lx*sign(x(i));
        
    end   
    
end

Floe0=[Floe0 ghostFloeX];

x=cat(1,Floe0.Xi);
y=cat(1,Floe0.Yi);
alive=cat(1,Floe0.alive);

for i=1:length(Floe0)
    poly = polyshape(Floe0(i).c_alpha'+[x(i) y(i)]);
    if alive(i) && (max(abs(poly.Vertices(:,2)))>Ly)
        
        ghostFloeY=[ghostFloeY  Floe0(i)];
        ghostFloeY(end).Yi=Floe0(i).Yi-2*Ly*sign(y(i));
        
    end
    
end

Floe0=[Floe0 ghostFloeY];
    

N0 = length(Floe0);
floenew = [];
KeepF = ones(length(Floe),1);
for ii = 1+Nb:length(Floe)
    if ~isempty(Floe(ii).interactions)
        %Identify which vertices are interacting with other floes to see if
        %it is eligible to be fractured off
        inter = Floe(ii).interactions(:,1);
        Xi = Floe(ii).interactions(inter<=N0,4);
        Yi = Floe(ii).interactions(inter<=N0,5);
        
        poly = polyshape(Floe(ii).c_alpha');
        polytrue = polyshape(Floe(ii).c_alpha'+[Floe(ii).Xi, Floe(ii).Yi]);
        if poly.NumRegions == 1 && poly.NumHoles == 0 
            if norm(poly.Vertices(1,:)-poly.Vertices(end,:)) == 0
                poly.Vertices(end,:) = [];
            end
            angles = polyangles(poly.Vertices(:,1),poly.Vertices(:,2));
            Anorm = 180-360/length(angles);%180;%
            break1=(rand(length(angles),1)>angles/Anorm);
            da = zeros(length(angles),1);
            [break2,~] = dsearchn(polytrue.Vertices,[Xi Yi]);
            in = zeros(length(polytrue.Vertices),1);
            bnd = max(isinf(inter));
            inter2 = inter(~isinf(inter));
            inter2(inter2>N0) = [];
            for jj = 1:length(inter2)
                c0 = Floe0(inter2(jj)).c_alpha'+[Floe0(inter2(jj)).Xi Floe0(inter2(jj)).Yi];
                [inn,~] = inpolygon(polytrue.Vertices(:,1),polytrue.Vertices(:,2),c0(:,1),c0(:,2));
                in = in+inn;
            end
            if bnd
                [inn,~] = inpolygon(polytrue.Vertices(:,1),polytrue.Vertices(:,2),c2_boundary(1,:)',c2_boundary(2,:)');
                in = in+~inn;
            end
            da(in>0) = 1;    
            da(break2) = 1;
            
            %Identify which corners to break
            grind = logical(break1 + da==2);
            
            %fracture off those corners
            if sum(grind)>1
                KeepF(ii) = 0;
                fracturedFloes = frac_corner(Floe(ii),grind,poly);
                floenew=[floenew fracturedFloes];
            end
        end
    end
end
if isfield(floenew,'poly')
    floenew=rmfield(floenew,{'poly'});
end
if isfield(Floe,'poly')
    Floe=rmfield(Floe,{'poly'});
end
if isfield(floenew,'potentialInteractions')
    floenew=rmfield(floenew,{'potentialInteractions'});
end
if isfield(Floe,'potentialInteractions')
    Floe=rmfield(Floe,{'potentialInteractions'});
end

Floe = [Floe(logical(KeepF)) floenew];

warning('on',id)
warning('on',id3)

end

